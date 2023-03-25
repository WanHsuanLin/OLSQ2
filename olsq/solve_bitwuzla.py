import math
import datetime
from pybitwuzla import *
from olsq.input import input_qasm
from olsq.output import output_qasm
from olsq.device import qcdevice
from olsq.run_h_compiler import run_sabre
import pkgutil
from enum import Enum
import timeit

# pairwise    = 0
# seqcounter  = 1
# sortnetwrk  = 2
# cardnetwrk  = 3
# bitwise     = 4
# ladder      = 5
# totalizer   = 6
# mtotalizer  = 7
# kmtotalizer = 8
# native      = 9

class Mode(Enum):
    transition = 1
    normal = 2


def collision_extracting(list_gate_qubits):
    """Extract collision relations between the gates,
    If two gates g_1 and g_2 both acts on a qubit (at different time),
    we say that g_1 and g_2 collide on that qubit, which means that
    (1,2) will be in collision list.

    Args:
        list_gate_qubits: a list of gates in OLSQ IR
    
    Returns:
        list_collision: a list of collisions between the gates
    """

    list_collision = list()
    # We sweep through all the gates.  For each gate, we sweep through all the
    # gates after it, if they both act on some qubit, append them in the list.
    for g in range(len(list_gate_qubits)):
        for gg in range(g + 1, len(list_gate_qubits)):
            
            if list_gate_qubits[g][0] == list_gate_qubits[gg][0]:
                    list_collision.append((g, gg))
                
            if len(list_gate_qubits[gg]) == 2:
                if list_gate_qubits[g][0] == list_gate_qubits[gg][1]:
                    list_collision.append((g, gg))
            
            if len(list_gate_qubits[g]) == 2:
                if list_gate_qubits[g][1] == list_gate_qubits[gg][0]:
                    list_collision.append((g, gg))
                if len(list_gate_qubits[gg]) == 2:
                    if list_gate_qubits[g][1] == list_gate_qubits[gg][1]:
                        list_collision.append((g, gg))
    # print(list_collision)
    return tuple(list_collision)

def dependency_extracting(list_gate_qubits, count_program_qubit: int):
    """Extract dependency relations between the gates.
    If two gates g_1 and g_2 both acts on a qubit *and there is no gate
    between g_1 and g_2 that act on this qubit*, we then say that
    g2 depends on g1, which means that (1,2) will be in dependency list.

    Args:
        list_gate_qubits: a list of gates in OLSQ IR
        count_program_qubit: the number of logical/program qubit
    
    Returns:
        list_dependency: a list of dependency between the gates
    """

    list_dependency = []
    list_last_gate = [-1 for i in range(count_program_qubit)]
    # list_last_gate records the latest gate that acts on each qubit.
    # When we sweep through all the gates, this list is updated and the
    # dependencies induced by the update is noted.
    for i, qubits in enumerate(list_gate_qubits):
        
        if list_last_gate[qubits[0]] >= 0:
            list_dependency.append((list_last_gate[qubits[0]], i))
        list_last_gate[qubits[0]] = i

        if len(qubits) == 2:
            if list_last_gate[qubits[1]] >= 0:
                list_dependency.append((list_last_gate[qubits[1]], i))
            list_last_gate[qubits[1]] = i
    print("list_dependency")
    print(list_dependency)
    return tuple(list_dependency)


class OLSQ:
    def __init__(self, obj_is_swap, mode, encoding, swap_up_bound = -1, thread = 1):
        """Set the objective of OLSQ, and whether it is transition-based

        Args:
            objective_name: can be "depth", "swap", or "fidelity"
            mode: can be "normal" or "transition" (TB-OLSQ in the paper)       
        """

        if mode == "transition":
            self.mode = Mode.transition
        elif mode == "normal":
            self.mode = Mode.normal
        else:
            raise ValueError("Invalid Choice of Transition-Based or Not")

        self.obj_is_swap = False
        self.obj_is_swap = obj_is_swap

        # These values should be updated in setdevice(...)
        self.device = None
        self.count_physical_qubit = 0
        self.list_qubit_edge = []
        self.swap_duration = 0
        self.dict_gate_duration = dict()
        # self.list_gate_duration = []

        # These values should be updated in setprogram(...)
        self.list_gate_qubits = []
        self.count_program_qubit = 0
        self.list_gate_name = []
        
        # bound_depth is a hyperparameter
        self.bound_depth = 0

        self.input_dependency = False
        self.list_gate_dependency = []
        self.card_encoding = encoding
        self.swap_up_bound = swap_up_bound
        self.thread = thread
        self.start = 0
        self.swap_sabre = 0
        # self.ancillary_var_counter = 0

    def setdevice(self, device: qcdevice):
        """Pass in parameters from the given device.  If in TB mode,
           swap_duration is set to 1 without modifying the device.

        Args:
            device: a qcdevice object for OLSQ
        """

        self.device = device
        self.count_physical_qubit = device.count_physical_qubit
        self.list_qubit_edge = device.list_qubit_edge
        self.swap_duration = device.swap_duration
        # if self.if_transition_based:
        if self.mode == Mode.transition:
            self.swap_duration = 1

    def setprogram(self, program, input_mode: str = None, gate_duration: dict = None):
        """Translate input program to OLSQ IR, and set initial depth
        An example of the intermediate representation is shown below.
        It contains three things: 1) the number of qubit in the program,
        2) a list of tuples representing qubit(s) acted on by a gate,
        the tuple has one index if it is a single-qubit gate,
        two indices if it is a two-qubit gate, and 3) a list of
        type/name of each gate, which is not important to OLSQ,
        and only needed when generating output.
        If in TB mode, initial depth=1; in normal mode, we perform ASAP
        scheduling without consideration of SWAP to calculate depth.

        Args:
            program: a qasm string, or a list of the three things in IR.
            input_mode: (optional) can be "IR" if the input has ben
                translated to OLSQ IR; can be "benchmark" to use one of
                the benchmarks.  Default mode assumes qasm input.

        Example:
            For the following circuit
                q_0: ───────────────────■───
                                        │  
                q_1: ───────■───────────┼───
                     ┌───┐┌─┴─┐┌─────┐┌─┴─┐
                q_2: ┤ H ├┤ X ├┤ TDG ├┤ X ├─
                     └───┘└───┘└─────┘└───┘ 
            count_program_qubit = 3
            gates = ((2,), (1,2), (2,), (0,1))
            gate_spec = ("h", "cx", "tdg", "cx")
        """
        
        if input_mode == "IR":
            self.count_program_qubit = program[0]
            self.list_gate_qubits = program[1]
            self.list_gate_name = program[2]
        elif input_mode == "benchmark":
            f = pkgutil.get_data(__name__, "benchmarks/" + program + ".qasm")
            program = input_qasm(f.decode("utf-8"))
            self.count_program_qubit = program[0]
            self.list_gate_qubits = program[1]
            self.list_gate_name = program[2]
        else:
            program = input_qasm(program)
            self.count_program_qubit = program[0]
            self.list_gate_qubits = program[1]
            self.list_gate_name = program[2]

        # create a dict to store the gate duration. Gate name to duration => given as spec
        self.dict_gate_duration = gate_duration
        # create a list to remember the gate duration. gate => duration => construct in when setting program and use for construct constraint.
        # if gate_duration != None:
        #     # self.list_gate_duration
        #     for gate_name in self.list_gate_name:
        #         self.list_gate_duration.append(self.dict_gate_duration[gate_name])
        # else:
        #     self.list_gate_duration = [1]*len(self.list_gate_qubits)

        # calculate the initial depth
        if self.mode == Mode.transition:
        # if self.if_transition_based:
            self.bound_depth = 1
        else:
            push_forward_depth = [0 for i in range(self.count_program_qubit)]
            for qubits in self.list_gate_qubits:
                if len(qubits) == 1:
                    push_forward_depth[qubits[0]] += 1
                else:
                    tmp_depth = push_forward_depth[qubits[0]]
                    if tmp_depth < push_forward_depth[qubits[1]]:
                        tmp_depth = push_forward_depth[qubits[1]]
                    push_forward_depth[qubits[1]] = tmp_depth + 1
                    push_forward_depth[qubits[0]] = tmp_depth + 1
            self.bound_depth = max(push_forward_depth)
        
        count_gate = len(self.list_gate_qubits)
        self.list_gate_two = []
        self.list_gate_single = []
        self.list_span_edge = None
        for l in range(count_gate):
            if len(self.list_gate_qubits[l]) == 1:
                self.list_gate_single.append(l)
            else:
                self.list_gate_two.append(l)

    def setdependency(self, dependency: list):
        """Specify dependency (non-commutation)

        Args:
            dependency: a list of gate index pairs
        
        Example:
            For the following circuit
                q_0: ───────────────────■───
                                        │  
                q_1: ───────■───────────┼───
                     ┌───┐┌─┴─┐┌─────┐┌─┴─┐
                q_2: ┤ H ├┤ X ├┤ TDG ├┤ X ├─
                     └───┘└───┘└─────┘└───┘ 
                gate   0    1     2     3
            dependency = [(0,1), (1,2), (2,3)]

            However, for this QAOA subcircuit (ZZ gates may have phase
            parameters, but we neglect them for simplicity here)
                         ┌──┐ ┌──┐
                q_0: ────┤ZZ├─┤  ├─
                     ┌──┐└┬─┘ │ZZ│  
                q_1: ┤  ├─┼───┤  ├─
                     │ZZ│┌┴─┐ └──┘
                q_2: ┤  ├┤ZZ├──────
                     └──┘└──┘ 
                gate   0   1   2
            dependency = []    # since ZZ gates are commutable
        """
        self.list_gate_dependency = dependency
        self.input_dependency = True

    def dump(self, folder: str = None):
        """
        dump constraints for OLSQ mode
        """
            
        return


    def solve(self, use_sabre, output_mode: str = None, output_file_name: str = None):
        """Formulate an SMT, pass it to z3 solver, and output results.
        CORE OF OLSQ, EDIT WITH CARE.

        Returns:
            a pair of int that specifies the upper and lower bound of SWAP gate counts
            a list of results depending on output_mode
            "IR": 
            | list_scheduled_gate_name: name/type of each gate
            | list_scheduled_gate_qubits: qubit(s) each gate acts on
            | initial_mapping: logical qubit |-> physical qubit 
            | final_mapping: logical qubit |-> physical qubit in the end 
            | objective_value: depth/#swap/fidelity depending on setting
            None:
              a qasm string
              final_mapping
              objective_value
        """
        if self.count_physical_qubit < self.count_program_qubit:
            raise ValueError("[ERROR] number of physical qubits is less than number of program qubits")
        self._preprocessing()
        if self.mode == Mode.transition:
            print("Using transition based mode...")
            results = self._solve(use_sabre, output_mode, output_file_name, memory_max_size, verbose)
        elif self.mode == Mode.normal:
            print("Using normal mode...")
            results = self._solve(use_sabre, output_mode, output_file_name, memory_max_size, verbose)
        else:
            raise ValueError( ("Wrong type") )
        return results 

    def _preprocessing(self):
        if not self.input_dependency:
            self.list_gate_dependency = collision_extracting(self.list_gate_qubits)
        # list_adjacency_qubit takes in a physical qubit index _p_, and
        # returns the list of indices of physical qubits adjacent to _p_
        list_adjacent_qubit = list()
        # list_span_edge takes in a physical qubit index _p_,
        # and returns the list of edges spanned from _p_
        list_qubit_edge = self.list_qubit_edge
        count_qubit_edge = len(list_qubit_edge)
        self.list_span_edge = list()
        for n in range(self.count_physical_qubit):
            list_adjacent_qubit.append(list())
            self.list_span_edge.append(list())
        for k in range(count_qubit_edge):
            list_adjacent_qubit[list_qubit_edge[k][0]].append(
                                                        list_qubit_edge[k][1])
            list_adjacent_qubit[list_qubit_edge[k][1]].append(
                                                        list_qubit_edge[k][0])
            self.list_span_edge[list_qubit_edge[k][0]].append(k)
            self.list_span_edge[list_qubit_edge[k][1]].append(k)
        

    def _solve(self, use_sabre = False, output_mode: str = None, output_file_name: str = None):
        

        # pre-processing
        count_gate = len(self.list_gate_qubits)
        

        not_solved = True
        model = None
        iteration = 1
        # tight_bound_depth: use for depth constraint

        # bound_depth: generate constraints until t = bound_depth
        self.start = timeit.default_timer()
        if self.mode == Mode.transition:
            # bound_depth = self.get_swap_upper_bound() + 1
            tight_bound_depth = 0
            bound_depth = 3
        else:
            tight_bound_depth = self.bound_depth
            bound_depth = int(1.5 * self.bound_depth)
        # bound_depth = 21
        self.swap_sabre, max_depth = self.get_swap_upper_bound()


        while not_solved:
            print("start adding constraints...")
            # variable setting 
            lsqc = Bitwuzla()
            lsqc.set_option(Option.INCREMENTAL, True)
            pi, time, sigma = self._construct_variable(lsqc, bound_depth)
            
            # constraint setting
            self._add_injective_mapping_constraints(bound_depth, pi, lsqc)

            # Consistency between Mapping and Space Coordinates.
            self._add_consistency_gate_constraints(bound_depth, pi, time, lsqc)
            
            # Avoiding Collisions and Respecting Dependencies. 
            self._add_dependency_constraints(lsqc, time, bound_depth)

            # # No swap for t<s
            # # swap gates can not overlap with swap
            if self.mode == Mode.transition:
                self._add_swap_constraints(bound_depth, sigma, lsqc)
            else:
                self._add_swap_constraints(bound_depth, sigma, lsqc, True, time, pi)
            # Mapping Not Transformations by SWAP Gates.
            # Mapping Transformations by SWAP Gates.
            self._add_transformation_constraints(bound_depth, lsqc, sigma, pi)
            if self.mode == Mode.transition:
                not_solved = self._optimize_circuit_tran(tight_bound_depth, lsqc , time, sigma, count_gate, bound_depth, use_sabre)
            else:
                not_solved = self._optimize_circuit_normal(tight_bound_depth, lsqc , time, sigma, count_gate, bound_depth, use_sabre)
            if not_solved:
                tight_bound_depth = bound_depth
                bound_depth = max(2*bound_depth, max_depth)
            else:
                result = self._extract_results(lsqc, time, pi, sigma, output_mode, output_file_name)
                
        print(f"Total compilation time = {timeit.default_timer() - self.start}.")
        return result

    def _construct_variable(self, lsqc, bound_depth):
        # at cycle t, logical qubit q is mapped to pi[q][t]
        
        count_qubit_edge = len(self.list_qubit_edge)
        length = int(math.log2(self.count_physical_qubit))+1
        # length = int(math.log2(max(count_qubit_edge, self.count_physical_qubit)))+1
        pi_sort = lsqc.mk_bv_sort(length)
        pi = [[lsqc.mk_const(pi_sort, "map_q{}_t{}".format(i, j)) for j in range(bound_depth)]
                for i in range(self.count_program_qubit)]

        length = int(math.log2(bound_depth))+1
        # time coordinate for gate l is time[l]
        time_sort = lsqc.mk_bv_sort(length)
        time = [lsqc.mk_const(time_sort, "time_{}".format(i)) for i in range(len(self.list_gate_qubits))]

        # if at cycle t, a SWAP finishing on edge k, then sigma[k][t]=1
        bool_sort = lsqc.mk_bool_sort()
        sigma = [[lsqc.mk_const(bool_sort, "ifswap_e{}_t{}".format(i, j))
            for j in range(bound_depth)] for i in range(count_qubit_edge)]

        return pi, time, sigma


    def _add_injective_mapping_constraints(self, bound_depth, pi, model):
        # Injective Mapping
        for t in range(bound_depth):
            model.assert_formula(model.mk_term(Kind.DISTINCT, [pi[m][t] for m in range(self.count_program_qubit)]))
            for m in range(self.count_program_qubit):
                model.assert_formula(model.mk_term(Kind.BV_UGE, [pi[m][t], 0]))
                model.assert_formula(model.mk_term(Kind.BV_ULT, [pi[m][t], self.count_physical_qubit]))

    def _add_consistency_gate_constraints(self, bound_depth, pi, time, model):
        # Consistency between Mapping and Space Coordinates.
        list_gate_qubits = self.list_gate_qubits
        
        for l in self.list_gate_two:
            for t in range(bound_depth):
                term1 = model.mk_term(Kind.NOT, [model.mk_term(Kind.EQUAL, [time[l], t])])
                term2 = model.mk_term(Kind.OR, [model.mk_term(Kind.OR, [model.mk_term(Kind.AND, [model.mk_term(Kind.EQUAL, [pi[list_gate_qubits[l][0]][t], edge[0]]),
                                                                                                  model.mk_term(Kind.EQUAL, [pi[list_gate_qubits[l][1]][t], edge[1]]) ]),
                                                                        model.mk_term(Kind.AND, [model.mk_term(Kind.EQUAL, [pi[list_gate_qubits[l][0]][t], edge[1]]),
                                                                                                 model.mk_term(Kind.EQUAL, [pi[list_gate_qubits[l][1]][t], edge[0]]) ])]) for edge in self.list_qubit_edge ])
                model.assert_formula(model.mk_term(Kind.OR, [term1, term2]))
                
    def _add_dependency_constraints(self, model, time, bound_depth):
        # list_gate_duration = self.list_gate_duration
        list_gate_dependency = self.list_gate_dependency
        count_gate = len(self.list_gate_qubits)
        if self.mode == Mode.transition:
            for d in list_gate_dependency:
                model.assert_formula(model.mk_term(Kind.BV_ULE, [time[d[0]],time[d[1]]]))
        else:
            # length = int(math.log2(bound_depth))+1
            # bit_duration_list = [BitVecVal(list_gate_duration[l], length) for l in range(count_gate)]
            # bit_duration_minus_one_list = [BitVecVal(list_gate_duration[l]-1, length) for l in range(count_gate)]

            for d in list_gate_dependency:
                model.assert_formula(model.mk_term(Kind.BV_ULT, [time[d[0]],time[d[1]]]))
            # add initial condition for gates
            for l in range(count_gate):
                model.add(ULE(0, time[l]))
                model.assert_formula(model.mk_term(Kind.BV_ULE, [0, time[l]]))
                # model.add(ULE(bit_duration_minus_one_list[l], time[l]))

    def _add_swap_constraints(self, bound_depth, sigma, model, normal = False, time = None, pi = None):
        # if_overlap_edge takes in two edge indices _e_ and _e'_,
        # and returns whether or not they overlap
        list_qubit_edge = self.list_qubit_edge
        count_qubit_edge = len(list_qubit_edge)
        list_gate_qubits = self.list_gate_qubits

        # No swap for t<s
        for t in range(min(self.swap_duration - 1, bound_depth)):
            for k in range(count_qubit_edge):
                model.assert_formula(model.mk_term(Kind.NOT, [sigma[k][t]]))
        # swap gates can not overlap with swap
        for t in range(self.swap_duration - 1, bound_depth):
            for k in range(count_qubit_edge):
                for tt in range(t - self.swap_duration + 1, t):
                    model.assert_formula(model.mk_term(Kind.OR, [model.mk_term(Kind.NOT, [sigma[k][t]]), model.mk_term(Kind.NOT, [sigma[k][tt]])]))
                for tt in range(t - self.swap_duration + 1, t + 1):
                    for i in range(2):
                        for kk in self.list_span_edge[list_qubit_edge[k][i]]:
                            if kk != k:
                                model.assert_formula(model.mk_term(Kind.OR, [model.mk_term(Kind.NOT, [sigma[k][t]]), model.mk_term(Kind.NOT, [sigma[kk][tt]])]))

        if normal:
            # swap gates can not ovelap with other gates
            # the time interval should be modified
            count_gate = len(self.list_gate_qubits)
            for t in range(self.swap_duration - 1, bound_depth):
                for k in range(count_qubit_edge):
                    for l in range(count_gate):
                        # for tt in range(t - self.swap_duration + 1, t + self.list_gate_duration[l]):

                        for tt in range(t - self.swap_duration + 1, t + 1):
                            if l in self.list_gate_single:
                                term_or = model.mk_term(Kind.OR, [model.mk_term(Kind.EQUAL, [pi[list_gate_qubits[l][0]][tt], list_qubit_edge[k][0]]), 
                                                                  model.mk_term(Kind.EQUAL, [pi[list_gate_qubits[l][0]][tt], list_qubit_edge[k][1]])])
                                term_and = model.mk_term(Kind.AND, [model.mk_term(Kind.EQUAL, [time[l], tt]), term_or])
                                model.assert_formula(model.mk_term(Kind.IMPLIES, [term_and, model.mk_term(Kind.NOT, [sigma[k][t]])]))
                                
                            elif l in self.list_gate_two:
                                term_or = model.mk_term(Kind.OR, [model.mk_term(Kind.EQUAL, [pi[list_gate_qubits[l][0]][tt], list_qubit_edge[k][0]]), 
                                                                  model.mk_term(Kind.EQUAL, [pi[list_gate_qubits[l][0]][tt], list_qubit_edge[k][1]]),
                                                                  model.mk_term(Kind.EQUAL, [pi[list_gate_qubits[l][1]][tt], list_qubit_edge[k][0]]),
                                                                  model.mk_term(Kind.EQUAL, [pi[list_gate_qubits[l][1]][tt], list_qubit_edge[k][1]])])
                                term_and = model.mk_term(Kind.AND, [model.mk_term(Kind.EQUAL, [time[l], tt]), term_or])
                                model.assert_formula(model.mk_term(Kind.IMPLIES, [term_and, model.mk_term(Kind.NOT, [sigma[k][t]])]))
     
    def _add_transformation_constraints(self, bound_depth, model, sigma, pi):
        list_qubit_edge = self.list_qubit_edge
        count_qubit_edge = len(list_qubit_edge)
        list_span_edge = self.list_span_edge
        # Mapping Not Transformations by SWAP Gates.
        for t in range(bound_depth - 1):
            for n in range(self.count_physical_qubit):
                for m in range(self.count_program_qubit):
                    term_or = model.mk_term(Kind.OR, [sigma[k][t] for k in list_span_edge[n]])
                    term_and = model.mk_term(Kind.AND, [model.mk_term(Kind.NOT, [term_or]), model.mk_term(Kind.EQUAL, [pi[m][t], n])])
                    model.assert_formula(model.mk_term(Kind.IMPLIES, [term_and, model.mk_term(Kind.EQUAL, [pi[m][t+1], n])]))
        
        # Mapping Transformations by SWAP Gates.
        for t in range(bound_depth - 1):
            for k in range(count_qubit_edge):
                for m in range(self.count_program_qubit):
                    term_and = model.mk_term(Kind.AND, [sigma[k][t], model.mk_term(Kind.EQUAL, [pi[m][t], list_qubit_edge[k][0]])])
                    model.assert_formula(model.mk_term(Kind.IMPLIES, [term_and, model.mk_term(Kind.EQUAL, [pi[m][t+1], list_qubit_edge[k][1]])]))
                    term_and = model.mk_term(Kind.AND, [sigma[k][t], model.mk_term(Kind.EQUAL, [pi[m][t], list_qubit_edge[k][1]])])
                    model.assert_formula(model.mk_term(Kind.IMPLIES, [term_and, model.mk_term(Kind.EQUAL, [pi[m][t+1], list_qubit_edge[k][0]])]))

           
    def _add_atmostk_cnf(self, model, sigma, k, tight_bound_depth):
        from pysat.card import CardEnc
        num_sigma = 1+(tight_bound_depth+1)*len(self.list_qubit_edge)
        sigma_list = [i for i in range(1,num_sigma)]
        # print(len(sigma_list))
        print("Using encoding mehtod {}".format(self.card_encoding))
        cnf = CardEnc.atmost(lits = sigma_list, bound = k, encoding = self.card_encoding)
        # and_list = []
        bool_sort = lsqc.mk_bool_sort()
        ancillary = dict()
        for c in cnf.clauses:
            or_list = []
            for l in c:
                var = abs(l)
                if var < num_sigma:
                    sigma_idx1 = (var - 1)//(1+tight_bound_depth)
                    sigma_idx2 = (var - 1) % (1+tight_bound_depth)
                    if l < 0:
                        term = model.mk_term(Kind.NOT, [sigma[sigma_idx1][sigma_idx2]])
                        or_list.append(term)
                    else:
                        or_list.append(sigma[sigma_idx1][sigma_idx2])
                else:
                    if var not in ancillary.keys():
                        ancillary[var] = model.mk_const(bool_sort, "anx_{}".format(var))
                    if l < 0:
                        term = model.mk_term(Kind.NOT, [ancillary[var]])
                        or_list.append(term)
                    else:
                        or_list.append(ancillary[var])
            # print("or_list: ", or_list)
            # input()
            # and_list.append(Or(or_list))
            if len(or_list) > 1:
                model.assert_formula(model.mk_term(Kind.OR, or_list))
            else:
                model.assert_formula(or_list[0])
        # self.ancillary_var_counter = max_idx + 1
        # print("clauses num: {}".format(len(cnf.clauses)))
        # print("anxillary num: {}".format(len(ancillary)))
        # print(and_list)
        # model.add(And(and_list))

    def _count_swap(self, model, sigma, result_depth):
        n_swap = 0
        for k in range(len(self.list_qubit_edge)):
            for t in range(result_depth):
                if model.get_value_str(sigma[k][t]) == '1':
                    n_swap += 1
        return n_swap

    def _optimize_circuit_tran(self, tight_bound_depth, lsqc, time, sigma, count_gate, bound_depth, use_sabre):
        if use_sabre:
            swap_sabre = self.swap_sabre
            if swap_sabre > 0:
                upper_b_swap = swap_sabre - 1
        elif self.swap_up_bound > -1:
            upper_b_swap = self.swap_up_bound
        else:
            upper_b_swap = len(self.list_gate_two)
        print("set initial swap bound {}".format(upper_b_swap))
        lower_b_swap = 0
        bound_swap_num = 0
        find_min_depth = False
        # incremental solving use pop and push

        n_swap = -1
        length = int(math.log2(bound_depth))+1
        bv_sort = lsqc.mk_bv_sort(length)
        bit_tight_bound_depth = None
        term = None
        while not find_min_depth:
            bit_tight_bound_depth = lsqc.mk_bv_value(bv_sort, tight_bound_depth)
            print("Trying maximal depth = {}...".format(tight_bound_depth))
            start_time = datetime.datetime.now()
            # for depth optimization
            term = [model.mk_term(Kind.BV_UGE, [bit_tight_bound_depth, time[l]]) for l in range(count_gate)]
            lsqc.assume_formula(term)
            satisfiable = lsqc.check_sat()
            print("Depth optimization time = {}, time including preprocessing = {}.".format(datetime.datetime.now() - start_time, timeit.default_timer()-self.start))
            if satisfiable == Result.SAT:
                find_min_depth = True
                model = lsqc.model()
                n_swap = self._count_swap(model, sigma, tight_bound_depth)
                upper_b_swap = min(n_swap-1, upper_b_swap)
                bound_swap_num = upper_b_swap
                lower_b_swap = tight_bound_depth
                print("Find minimal depth {} with swap num {}".format(tight_bound_depth, n_swap))
            else:
                # lsqc.pop()
                tight_bound_depth += 1
                if tight_bound_depth >= bound_depth:
                    print("FAIL to find depth witnin {}.".format(bound_depth-1))
                    break
        if not find_min_depth:
            return True
        
        # stats = lsqc.statistics()
        # print(stats)
        lsqc.assert_formula(term)

        # for swap optimization
        if n_swap == 0:
            find_min_swap = True
            not_solved = False
        else:
            find_min_swap = False
        while not find_min_swap:
            print("Bound of Trying min swap = {}...".format(bound_swap_num))
            start_time = datetime.datetime.now()
            lsqc.push()
            # TODO: add atmost-k constraint
            self._add_atmostk_cnf(lsqc, sigma, bound_swap_num, tight_bound_depth)
            satisfiable = lsqc.check_sat()
            print("status:{}, optimization time = {}, time including preprocessing = {}".format(satisfiable, datetime.datetime.now() - start_time, timeit.default_timer()-self.start))
            if satisfiable == Result.SAT:
                cur_swap = self._count_swap(lsqc, sigma, tight_bound_depth)
                # print(cur_swap, lower_b_swap)
                if cur_swap > lower_b_swap:
                    upper_b_swap = cur_swap
                    if use_sabre:
                        bound_swap_num = upper_b_swap - 1
                    else:
                        bound_swap_num = (upper_b_swap + lower_b_swap) // 2
                else: 
                    find_min_swap = True
                    not_solved = False
            else:
                lsqc.pop()
                lower_b_swap = bound_swap_num + 1
                if upper_b_swap <= lower_b_swap:
                    start_time = datetime.datetime.now()
                    print("Trying min swap = {}...".format(upper_b_swap))
                    # TODO: add atmost-k constraint
                    lsqc.push()
                    self._add_atmostk_cnf(lsqc, sigma, upper_b_swap, tight_bound_depth)
                    satisfiable = lsqc.check_sat()
                    # satisfiable = lsqc.check(PbLe([(sigma[k][t],1) for k in range(count_qubit_edge)
                    #      for t in range(bound_depth)], upper_b_swap) ):;q
                    print("status:{}, optimization time = {}, time including preprocessing = {}".format(satisfiable, datetime.datetime.now() - start_time, timeit.default_timer()-self.start))
                    assert(satisfiable == Result.SAT)
                    find_min_swap = True
                    not_solved = False
                elif not use_sabre:
                    bound_swap_num = (upper_b_swap + lower_b_swap)//2
            lsqc.pop()
        return not_solved
    

    def _optimize_circuit_normal(self, tight_bound_depth, lsqc, time, sigma, count_gate, bound_depth, use_sabre):
        if use_sabre:
            swap_sabre = self.swap_sabre
            if swap_sabre > 0:
                upper_b_swap = swap_sabre - 1
        elif self.swap_up_bound > -1:
            print("set initial swap bound {}".format(self.swap_up_bound))
            upper_b_swap = self.swap_up_bound
        else:
            upper_b_swap = len(self.list_gate_two)
        lower_b_swap = 0
        bound_swap_num = 0
        find_min_depth = False
        # incremental solving use pop and push
        length = int(math.log2(bound_depth))+1
        bv_sort = lsqc.mk_bv_sort(length)
        bit_tight_bound_depth = None
        term = None
        n_swap = -1
        step = 5
        while not find_min_depth:
            bit_tight_bound_depth = lsqc.mk_bv_value(bv_sort, tight_bound_depth)
            print("Trying maximal depth = {}...".format(tight_bound_depth))
            term = [model.mk_term(Kind.BV_UGE, [bit_tight_bound_depth, time[l]]) for l in range(count_gate)]
            lsqc.assume_formula(term)
            satisfiable = lsqc.check_sat()
            print("Depth optimization time = {}, time including preprocessing = {}.".format(datetime.datetime.now() - start_time, timeit.default_timer()-self.start))
            if satisfiable == Result.SAT:
                find_min_depth = True
                n_swap = self._count_swap(lsqc, sigma, tight_bound_depth)
                upper_b_swap = min(n_swap-1, upper_b_swap)
                bound_swap_num = upper_b_swap
                print("Find minimal depth {} with swap num {}".format(tight_bound_depth, n_swap))
            else:
                # lsqc.pop()
                if tight_bound_depth > 100: 
                    step = 10
                tight_bound_depth = step + tight_bound_depth
                if tight_bound_depth > bound_depth:
                    print("FAIL to find depth witnin {}.".format(bound_depth))
                    break
        if not find_min_depth:
            return True
        if not self.obj_is_swap:
            return False
        # stats = lsqc.statistics()
        # print(stats)
        lsqc.assert_formula(term)

        # for swap optimization
        if n_swap == 0:
            find_min_swap = True
            not_solved = False
        else:
            find_min_swap = False
        while not find_min_swap:
            print("Bound of Trying min swap = {}...".format(bound_swap_num))
            start_time = datetime.datetime.now()
            lsqc.push()
            # TODO: add atmost-k constraint
            self._add_atmostk_cnf(lsqc, sigma, bound_swap_num, tight_bound_depth)
            # lsqc.add(PbLe([(sigma[k][t],1) for k in range(count_qubit_edge)
            #              for t in range(bound_depth)], bound_swap_num) )
            satisfiable = lsqc.check_sat()
            # satisfiable = lsqc.check(PbLe([(sigma[k][t],1) for k in range(count_qubit_edge)
            #              for t in range(bound_depth)], bound_swap_num) )
            # lsqc.pop()
            print("status:{}, optimization time = {}, time including preprocessing = {}".format(satisfiable, datetime.datetime.now() - start_time, timeit.default_timer()-self.start))
            if satisfiable == Result.SAT:
                cur_swap = self._count_swap(lsqc, sigma, tight_bound_depth)
                # print(cur_swap, lower_b_swap)
                if cur_swap > lower_b_swap:
                    upper_b_swap = cur_swap
                    if use_sabre:
                        bound_swap_num = upper_b_swap - 1
                    else:
                        bound_swap_num = (upper_b_swap + lower_b_swap) // 2
                else: 
                    find_min_swap = True
                    not_solved = False
            else:
                lsqc.pop()
                lower_b_swap = bound_swap_num + 1
                if upper_b_swap <= lower_b_swap:
                    start_time = datetime.datetime.now()
                    if upper_b_swap == swap_sabre - 1:
                        upper_b_swap = swap_sabre
                    print("Trying min swap = {}...".format(upper_b_swap))
                    lsqc.push()
                    self._add_atmostk_cnf(lsqc, sigma, upper_b_swap, tight_bound_depth)
                    satisfiable = lsqc.check_sat()
                    # satisfiable = lsqc.check(PbLe([(sigma[k][t],1) for k in range(count_qubit_edge)
                    #      for t in range(bound_depth)], upper_b_swap) ):;q
                    print("status:{}, optimization time = {}, time including preprocessing = {}".format(satisfiable, datetime.datetime.now() - start_time, timeit.default_timer()-self.start))
                    assert(satisfiable == Result.SAT)
                    find_min_swap = True
                    not_solved = False
                elif not use_sabre:
                    bound_swap_num = (upper_b_swap + lower_b_swap)//2
            lsqc.pop()
        return not_solved

    
    def _extract_results(self, model, time, pi, sigma, output_mode, output_file_name):
        # post-processing
        list_gate_two = self.list_gate_two
        list_gate_single = self.list_gate_single
        list_qubit_edge = self.list_qubit_edge
        list_gate_qubits = self.list_gate_qubits
        count_qubit_edge = len(list_qubit_edge)
        count_gate = len(list_gate_qubits)
        list_gate_name = self.list_gate_name
        count_program_qubit = self.count_program_qubit
        result_time = []
        result_depth = 0
        for l in range(count_gate):
            result_time.append(self._get_model_bv(model, time[l]))
            result_depth = max(result_depth, result_time[-1])
        result_depth += 1
        list_result_swap = []
        # print(m.evaluate(f(x)))
        for k in range(count_qubit_edge):
            for t in range(result_depth):
                if model.get_value_str(sigma[k][t]) == '1':
                    list_result_swap.append((k, t))
                    print(f"SWAP on physical edge ({list_qubit_edge[k][0]},"\
                        + f"{list_qubit_edge[k][1]}) at time {t}")
        for l in range(count_gate):
            if len(list_gate_qubits[l]) == 1:
                qq = list_gate_qubits[l][0]
                tt = result_time[l]
                print(f"Gate {l}: {list_gate_name[l]} {qq} on qubit "\
                    + f"{self._get_model_bv(model, pi[qq][tt])} at time {tt}")
            else:
                qq = list_gate_qubits[l][0]
                qqq = list_gate_qubits[l][1]
                tt = result_time[l]
                print(f"Gate {l}: {list_gate_name[l]} {qq}, {qqq} on qubits "\
                    + f"{self._get_model_bv(model, pi[qq][tt])} and "\
                    + f"{self._get_model_bv(model, pi[qqq][tt])} at time {tt}")

        tran_detph = result_depth

        # transition based
        real_time = None
        if self.mode == Mode.transition:
            real_time = [0] * count_gate
            self.swap_duration = self.device.swap_duration
            list_depth_on_qubit = [-1] * self.count_physical_qubit
            list_real_swap = []
            for block in range(result_depth):
                for tmp_gate in range(count_gate):
                    if result_time[tmp_gate] == block:
                        qubits = list_gate_qubits[tmp_gate]
                        if len(qubits) == 1:
                            p0 = self._get_model_bv(model, pi[qubits[0]][block])
                            real_time[tmp_gate] = \
                                list_depth_on_qubit[p0] + 1
                            list_depth_on_qubit[p0] = \
                                real_time[tmp_gate]
                        else:
                            p0 = self._get_model_bv(model, pi[qubits[0]][block])
                            p1 = self._get_model_bv(model, pi[qubits[1]][block])
                            real_time[tmp_gate] = max(
                                list_depth_on_qubit[p0],
                                list_depth_on_qubit[p1]) + 1
                            list_depth_on_qubit[p0] = \
                                real_time[tmp_gate]
                            list_depth_on_qubit[p1] = \
                                real_time[tmp_gate]
                            # print(f"{tmp_gate} {p0} {p1} real-time={real_time[tmp_gate]}")

                if block < result_depth - 1:
                    for (k, t) in list_result_swap:
                        if t == block:
                            p0 = list_qubit_edge[k][0]
                            p1 = list_qubit_edge[k][1]
                            tmp_time = max(list_depth_on_qubit[p0],
                                list_depth_on_qubit[p1]) \
                                + self.swap_duration
                            list_depth_on_qubit[p0] = tmp_time
                            list_depth_on_qubit[p1] = tmp_time
                            list_real_swap.append((k, tmp_time))
                # print(list_depth_on_qubit)
            real_depth = 0
            for tmp_depth in list_depth_on_qubit:
                if real_depth < tmp_depth + 1:
                    real_depth = tmp_depth + 1
            result_depth = real_depth
            list_result_swap = list_real_swap
            # print(list_result_swap)
        else:
            real_time = result_time

        print(f"result- additional SWAP count = {len(list_result_swap)}.")
        print(f"result- circuit depth = {result_depth}.")

        list_scheduled_gate_qubits = [[] for i in range(result_depth)]
        list_scheduled_gate_name = [[] for i in range(result_depth)]
        for l in range(count_gate):
            t = result_time[l]
            list_scheduled_gate_name[t].append(list_gate_name[l])
            if l in list_gate_single:
                q = self._get_model_bv(model, pi[list_gate_qubits[l][0]][result_time[l]])
                list_scheduled_gate_qubits[t].append((q,))
            elif l in list_gate_two:
                [q0, q1] = list_gate_qubits[l]
                tmp_t = t
                if self.mode == Mode.transition:
                    tmp_t = self._get_model_bv(model, time[l])
                q0 = self._get_model_bv(model, pi[q0][tmp_t])
                q1 = self._get_model_bv(model, pi[q1][tmp_t])
                # q0 = model.evaluate(f_map(q0,tmp_t)).as_long()
                # q1 = model.evaluate(f_map(q1,tmp_t)).as_long()
                list_scheduled_gate_qubits[t].append((q0, q1))
            else:
                raise ValueError("Expect single-qubit or two-qubit gate.")

        tmp_depth = result_depth - 1
        if self.mode == Mode.transition:
            tmp_depth = tran_detph - 1
        final_mapping = [self._get_model_bv(model, pi[m][tmp_depth]) for m in range(count_program_qubit)]
        initial_mapping = [self._get_model_bv(model, pi[m][0]) for m in range(count_program_qubit)]

        # final_mapping = [model.evaluate(f_map(m,tmp_depth)).as_long() for m in range(count_program_qubit)]
        # initial_mapping = [model.evaluate(f_map(m,0)).as_long() for m in range(count_program_qubit)]

        for (k, t) in list_result_swap:
            q0 = list_qubit_edge[k][0]
            q1 = list_qubit_edge[k][1]
            if self.swap_duration == 1:
                list_scheduled_gate_qubits[t].append((q0, q1))
                list_scheduled_gate_name[t].append("SWAP")
            elif self.swap_duration == 3:
                list_scheduled_gate_qubits[t].append((q0, q1))
                list_scheduled_gate_name[t].append("cx")
                list_scheduled_gate_qubits[t - 1].append((q1, q0))
                list_scheduled_gate_name[t - 1].append("cx")
                list_scheduled_gate_qubits[t - 2].append((q0, q1))
                list_scheduled_gate_name[t - 2].append("cx")
            else:
                raise ValueError("Expect SWAP duration one, or three")


        if output_mode == "IR":
            if output_file_name:
                output_file = open(output_file_name, 'w')
                output_file.writelines([list_scheduled_gate_name,
                                        list_scheduled_gate_qubits,
                                        final_mapping])
            return (result_depth,
                    list_scheduled_gate_name,
                    list_scheduled_gate_qubits,
                    final_mapping,
                    initial_mapping)
        else:
            return (output_qasm(self.device, result_depth, list_scheduled_gate_name,
                                list_scheduled_gate_qubits, final_mapping,
                                True, output_file_name),
                    final_mapping,
                    initial_mapping)

    def get_swap_upper_bound(self, heuristic = "sabre"):
        if heuristic == "sabre":
            swap_num, depth, _ = run_sabre("olsq", self.list_gate_qubits, self.list_qubit_edge, self.count_physical_qubit)
            print("Run heuristic compiler sabre to get upper bound for SWAP: {}, depth: {}".format(swap_num, depth))
        else:
            raise TypeError("Only support sabre.")
        return swap_num, depth

    def _get_model_bv(self, model, var):
        r_str = model.get_value_str(var)
        return int(r_str, 2)
