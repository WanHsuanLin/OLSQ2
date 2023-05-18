import math
import datetime
from z3 import Bool, Implies, And, Or, sat, Solver, set_option, BitVec, ULT, ULE, UGE, BitVecVal, Not, Then, AtMost, Goal
from olsq.input import input_qasm
from olsq.output import output_qasm
from olsq.device import qcdevice
from olsq.run_h_compiler import run_sabre
import pkgutil
from enum import Enum
import timeit

TIMEOUT = 90000
# MEMORY_MAX_SIZE = 1000000000 * 58
MEMORY_MAX_SIZE = 0
MAX_TREAD_NUM = 8
VERBOSE = 10
# encoding 0 4 5 can only be applied to bound = 1 or len(lit) - 1
CARD_ENCODING = 1
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
    def __init__(self, obj_is_swap, mode, encoding, swap_up_bound = -1):
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

    def dump(self, folder: str = None, bound_depth = 5, bound_swap = 2):
        """
        dump constraints for OLSQ mode
        """
        # bound_depth: generate constraints until t = bound_depth
        print("start adding constraints...")
        # variable setting 
        self._preprocessing()
        pi, time, sigma = self._construct_variable(bound_depth)

        lsqc = Goal()
        start = timeit.default_timer()
        self._add_injective_mapping_constraints(bound_depth, pi, lsqc)
        self._add_consistency_gate_constraints(bound_depth, pi, time, lsqc)
        self._add_dependency_constraints(lsqc, time, bound_depth)
        self._add_swap_constraints(bound_depth, sigma, lsqc)
        self._add_transformation_constraints(bound_depth, lsqc, sigma, pi)
        lsqc.add([UGE(1, time[l]) for l in range(len(self.list_gate_qubits))])
        print("time to generate constraints: {}".format(timeit.default_timer() - start))
        self._add_atmostk_cnf(lsqc, sigma, bound_swap, bound_depth-1)
        start = timeit.default_timer()
        tactic = Then('simplify','propagate-values','solve-eqs','card2bv','bit-blast', 'tseitin-cnf')
        output_file_name = folder+"/"+str(self.count_physical_qubit)+"_"+str(self.count_program_qubit) + "_21.txt"
        cnf = tactic(lsqc)[0]
        print("time to generate cnf: {}".format(timeit.default_timer() - start))
        with open(output_file_name,"w") as ous:
            ous.write(cnf.dimacs())
        return


    def solve(self, use_sabre, output_mode: str = None, output_file_name: str = None, memory_max_size=MEMORY_MAX_SIZE, verbose = VERBOSE):
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
        

    def _solve(self, use_sabre = False, output_mode: str = None, output_file_name: str = None, memory_max_size=MEMORY_MAX_SIZE, verbose = VERBOSE):
        

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
            tight_bound_depth = 1
            bound_depth = 5
        else:
            tight_bound_depth = self.bound_depth
            bound_depth = int(1.5 * self.bound_depth)
        # bound_depth = 21
        
        self.swap_sabre, max_depth = self.get_swap_upper_bound()
        # return

        while not_solved:
            print("start adding constraints...")
            # variable setting 
            pi, time, sigma = self._construct_variable(bound_depth)
            # lsqc = Solver()
            lsqc = Then('simplify', 
                 'solve-eqs', 'bit-blast', 
                 'sat').solver()
            set_option("memory_max_size", memory_max_size)
            # set_option("timeout", TIMEOUT)

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
                not_solved, model = self._optimize_circuit_tran(tight_bound_depth, lsqc , pi, time, sigma, count_gate, bound_depth, use_sabre)
            else:
                not_solved, model = self._optimize_circuit_normal(tight_bound_depth, lsqc , pi, time, sigma, count_gate, bound_depth, use_sabre)
            if not_solved:
                tight_bound_depth = bound_depth
                bound_depth = max(2*bound_depth, max_depth)
            else:
                result = self._extract_results(model, time, pi, sigma, output_mode, output_file_name)

            stats = lsqc.statistics()
            print(stats)
                
        print(f"Total compilation time = {timeit.default_timer() - self.start}.")
        return result

    def _construct_variable(self, bound_depth):
        # at cycle t, logical qubit q is mapped to pi[q][t]
        
        count_qubit_edge = len(self.list_qubit_edge)
        length = int(math.log2(self.count_physical_qubit))+1
        # length = int(math.log2(max(count_qubit_edge, self.count_physical_qubit)))+1

        pi = [[BitVec(("map_q{}_t{}".format(i, j)), length) for j in range(bound_depth)]
                for i in range(self.count_program_qubit)]

        length = int(math.log2(bound_depth))+1
        # time coordinate for gate l is time[l]
        time = [BitVec("time_{}".format(i), length) for i in range(len(self.list_gate_qubits))]

        # if at cycle t, a SWAP finishing on edge k, then sigma[k][t]=1
        sigma = [[Bool("ifswap_e{}_t{}".format(i, j))
            for j in range(bound_depth)] for i in range(count_qubit_edge)]

        return pi, time, sigma

    def _add_injective_mapping_constraints(self, bound_depth, pi, model):
        # Injective Mapping
        for t in range(0, bound_depth):
            # model.add(Distinct([pi[m][t] for m in range(self.count_program_qubit)]))
            for m in range(self.count_program_qubit):
                model.add(UGE(pi[m][t], 0))
                model.add(ULT(pi[m][t], self.count_physical_qubit))
                for mm in range(m):
                    model.add(pi[m][t] != pi[mm][t])

    def _add_consistency_gate_constraints(self, bound_depth, pi, time, model):
        # Consistency between Mapping and Space Coordinates.
        list_gate_qubits = self.list_gate_qubits
        
        for l in self.list_gate_two:
            for t in range(bound_depth):
                model.add(Or(Not(time[l] == t), Or([Or(And(pi[list_gate_qubits[l][0]][t] == edge[0], pi[list_gate_qubits[l][1]][t] == edge[1]), \
                                And(pi[list_gate_qubits[l][0]][t] == edge[1], pi[list_gate_qubits[l][1]][t] == edge[0])) for edge in self.list_qubit_edge ] )))

    def _add_dependency_constraints(self, model, time, bound_depth):
        # list_gate_duration = self.list_gate_duration
        list_gate_dependency = self.list_gate_dependency
        count_gate = len(self.list_gate_qubits)
        if self.mode == Mode.transition:
            for d in list_gate_dependency:
                # lsqc.add(time[d[0]] <= time[d[1]])
                model.add(ULE(time[d[0]],time[d[1]]))
        else:
            # length = int(math.log2(bound_depth))+1
            # bit_duration_list = [BitVecVal(list_gate_duration[l], length) for l in range(count_gate)]
            # bit_duration_minus_one_list = [BitVecVal(list_gate_duration[l]-1, length) for l in range(count_gate)]

            for d in list_gate_dependency:
                model.add(ULT(time[d[0]], time[d[1]]))
            # add initial condition for gates
            for l in range(count_gate):
                model.add(ULE(0, time[l]))
                # model.add(ULE(bit_duration_minus_one_list[l], time[l]))

    def _add_swap_constraints(self, bound_depth, sigma, model, normal = False, time = None, pi = None):
        # if_overlap_edge takes in two edge indices _e_ and _e'_,
        # and returns whether or not they overlap
        list_qubit_edge = self.list_qubit_edge
        count_qubit_edge = len(list_qubit_edge)
        list_gate_qubits = self.list_gate_qubits
        # if_overlap_edge = [[0] * count_qubit_edge
        #     for k in range(count_qubit_edge)]
        # # list_over_lap_edge takes in an edge index _e_,
        # # and returnsthe list of edges that overlap with _e_
        # list_overlap_edge = list()
        # # list_count_overlap_edge is the list of lengths of
        # # overlap edge lists of all the _e_
        # # list_count_overlap_edge = list()
        # for k in range(count_qubit_edge):
        #     list_overlap_edge.append(list())
        # for k in range(count_qubit_edge):
        #     for kk in range(k + 1, count_qubit_edge):
        #         if (   (list_qubit_edge[k][0] == list_qubit_edge[kk][0]
        #                 or list_qubit_edge[k][0] == list_qubit_edge[kk][1])
        #             or (list_qubit_edge[k][1] == list_qubit_edge[kk][0]
        #                 or list_qubit_edge[k][1] == list_qubit_edge[kk][1]) ):
        #             list_overlap_edge[k].append(kk)
        #             list_overlap_edge[kk].append(k)
        #             if_overlap_edge[kk][k] = 1
        #             if_overlap_edge[k][kk] = 1
        # for k in range(count_qubit_edge):
        #     list_count_overlap_edge.append(len(list_overlap_edge[k]))

        # No swap for t<s
        for t in range(min(self.swap_duration - 1, bound_depth)):
            for k in range(count_qubit_edge):
                model.add(sigma[k][t] == False)
        # swap gates can not overlap with swap
        for t in range(self.swap_duration - 1, bound_depth):
            for k in range(count_qubit_edge):
                for tt in range(t - self.swap_duration + 1, t):
                    model.add(Or(Not(sigma[k][t]),
                        Not(sigma[k][tt])))
                for tt in range(t - self.swap_duration + 1, t + 1):
                    for i in range(2):
                        for kk in self.list_span_edge[list_qubit_edge[k][i]]:
                            if kk != k:
                                model.add(Or(Not(sigma[k][t]), 
                                            Not(sigma[kk][tt])))
                    # for kk in list_overlap_edge[k]:
                    #     model.add(Implies(sigma[k][t] == True,
                    #         sigma[kk][tt] == False))
                

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
                                model.add(Implies(And(time[l] == tt,
                                    Or(pi[list_gate_qubits[l][0]][tt] == list_qubit_edge[k][0],
                                        pi[list_gate_qubits[l][0]][tt] == list_qubit_edge[k][1])),
                                        sigma[k][t] == False             ))
                            elif l in self.list_gate_two:
                                model.add(Implies(And( \
                                    time[l] == tt, \
                                        Or(pi[list_gate_qubits[l][0]][tt] == list_qubit_edge[k][0], \
                                        pi[list_gate_qubits[l][0]][tt] == list_qubit_edge[k][1], \
                                        pi[list_gate_qubits[l][1]][tt] == list_qubit_edge[k][0], \
                                        pi[list_gate_qubits[l][1]][tt] == list_qubit_edge[k][1])), \
                                        sigma[k][t] == False           ))
        # else:
        #     for t in range(bound_depth-1):
        #             model.add(Or([And(sigma[k][t],
        #                              And([Not(sigma[kk][t]) for kk in itertools.chain(range(k), range(k+1, count_qubit_edge))]))
        #                          for k in range(count_qubit_edge)]))

     
    def _add_transformation_constraints(self, bound_depth, model, sigma, pi):
        list_qubit_edge = self.list_qubit_edge
        count_qubit_edge = len(list_qubit_edge)
        list_span_edge = self.list_span_edge
        # Mapping Not Transformations by SWAP Gates.
        for t in range(bound_depth - 1):
            for n in range(self.count_physical_qubit):
                for m in range(self.count_program_qubit):
                    model.add(
                        Implies(And(Not(Or([sigma[k][t] for k in list_span_edge[n]])),
                                pi[m][t] == n), pi[m][t + 1] == n))
        
        # Mapping Transformations by SWAP Gates.
        for t in range(bound_depth - 1):
            for k in range(count_qubit_edge):
                for m in range(self.count_program_qubit):
                    model.add(Implies(And(sigma[k][t] == True,
                        pi[m][t] == list_qubit_edge[k][0]),
                            pi[m][t + 1] == list_qubit_edge[k][1]))
                    model.add(Implies(And(sigma[k][t] == True,
                        pi[m][t] == list_qubit_edge[k][1]),
                            pi[m][t + 1] == list_qubit_edge[k][0]))
           
    def _add_atmostk_cnf(self, model, sigma, k, tight_bound_depth):
        from pysat.card import CardEnc
        num_sigma = 1+(tight_bound_depth+1)*len(self.list_qubit_edge)
        sigma_list = [i for i in range(1,num_sigma)]
        # print(len(sigma_list))
        print("Using encoding mehtod {}".format(self.card_encoding))
        cnf = CardEnc.atmost(lits = sigma_list, bound = k, encoding = self.card_encoding)
        # and_list = []
        ancillary = dict()
        # print(cnf.clauses)
        # print("sigma: ",len(sigma))
        # print("sigma[0]: ",len(sigma[0]))
        # print("tight_bound_depth: ",tight_bound_depth)
        # print("num_sigma: ", num_sigma)
        # max_idx = 0
        for c in cnf.clauses:
            # print("clause: ",c)
            or_list = []
            # print("c: ",c)
            for l in c:
                var = abs(l)
                # print(var)
                if var < num_sigma:
                    sigma_idx1 = (var - 1)//(1+tight_bound_depth)
                    sigma_idx2 = (var - 1) % (1+tight_bound_depth)
                    # print(sigma_idx1, sigma_idx2)
                    if l < 0:
                        or_list.append(Not(sigma[sigma_idx1][sigma_idx2]))
                    else:
                        or_list.append(sigma[sigma_idx1][sigma_idx2])
                else:
                    # print("anxi: ", anxillary)
                    # max_idx = max(max_idx, var)
                    # var = var + self.ancillary_var_counter
                    if var not in ancillary.keys():
                        ancillary[var] = Bool("anx_{}".format(var))
                    if l < 0:
                        or_list.append(Not(ancillary[var]))
                    else:
                        or_list.append(ancillary[var])
            # print("or_list: ", or_list)
            # input()
            # and_list.append(Or(or_list))
            model.add(Or(or_list))
        # self.ancillary_var_counter = max_idx + 1
        # print("clauses num: {}".format(len(cnf.clauses)))
        # print("anxillary num: {}".format(len(ancillary)))
        # print(and_list)
        # model.add(And(and_list))

    def _count_swap(self, model, sigma, result_depth):
        n_swap = 0
        for k in range(len(self.list_qubit_edge)):
            for t in range(result_depth):
                if model[sigma[k][t]]:
                    n_swap += 1
        return n_swap

    def _optimize_circuit_tran(self, tight_bound_depth, lsqc, pi, time, sigma, count_gate, bound_depth, use_sabre):
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
        bit_tight_bound_depth = None
        while not find_min_depth:
            bit_tight_bound_depth = BitVecVal(tight_bound_depth, length)
            print("Trying maximal depth = {}...".format(tight_bound_depth))
            start_time = datetime.datetime.now()
            # for depth optimization
            satisfiable = lsqc.check([UGE(bit_tight_bound_depth, time[l]) for l in range(count_gate)])
            print("Depth optimization time = {}, time including preprocessing = {}.".format(datetime.datetime.now() - start_time, timeit.default_timer()-self.start))
            if satisfiable == sat:
                find_min_depth = True
                model = lsqc.model()
                n_swap = self._count_swap(model, sigma, tight_bound_depth)
                upper_b_swap = min(n_swap-1, upper_b_swap)
                bound_swap_num = upper_b_swap
                lower_b_swap = tight_bound_depth
                # print("Find minimal depth {} with swap num {}".format(tight_bound_depth, model[count_swap].as_long()))
                print("Find minimal depth {} with swap num {}".format(tight_bound_depth, n_swap))
            else:
                # lsqc.pop()
                tight_bound_depth += 1
                if tight_bound_depth >= bound_depth:
                    print("FAIL to find solution with depth less than  {}.".format(bound_depth - 1))
                    break
        if not find_min_depth:
            return True, None
        
        # stats = lsqc.statistics()
        # print(stats)
        lsqc.add([UGE(bit_tight_bound_depth, time[l]) for l in range(count_gate)])

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
            # add atmost-k constraint
            self._add_atmostk_cnf(lsqc, sigma, bound_swap_num, tight_bound_depth)
            satisfiable = lsqc.check()
            print("status:{}, optimization time = {}, time including preprocessing = {}".format(satisfiable, datetime.datetime.now() - start_time, timeit.default_timer()-self.start))
            if satisfiable == sat:
                model = lsqc.model()
                # cur_swap = model[count_swap].as_long()
                cur_swap = self._count_swap(model, sigma, tight_bound_depth)
                # print(cur_swap, lower_b_swap)
                result = self._extract_results(model, time, pi, sigma, "qasm", "intermediate_result_swap_count_{}.qasm".format(cur_swap))
                if cur_swap > lower_b_swap:
                    upper_b_swap = cur_swap
                    if use_sabre:
                        bound_swap_num = upper_b_swap - 1
                    else:
                        bound_swap_num = (upper_b_swap + lower_b_swap) // 2
                else: 
                    find_min_swap = True
                    not_solved = False
                lsqc.pop()
            else:
                lsqc.pop()
                lower_b_swap = bound_swap_num + 1
                if upper_b_swap <= lower_b_swap:
                    while not find_min_swap:
                        start_time = datetime.datetime.now()
                        print("Trying min swap = {}...".format(lower_b_swap))
                        # TODO: add atmost-k constraint
                        lsqc.push()
                        self._add_atmostk_cnf(lsqc, sigma, lower_b_swap, tight_bound_depth)
                        satisfiable = lsqc.check()
                        # satisfiable = lsqc.check(PbLe([(sigma[k][t],1) for k in range(count_qubit_edge)
                        #      for t in range(bound_depth)], upper_b_swap) ):;q
                        print("status:{}, optimization time = {}, time including preprocessing = {}".format(satisfiable, datetime.datetime.now() - start_time, timeit.default_timer()-self.start))
                        if satisfiable == sat :
                            model = lsqc.model()
                            find_min_swap = True
                            not_solved = False
                        else:
                            lower_b_swap += 1
                        lsqc.pop()
                elif not use_sabre:
                    bound_swap_num = (upper_b_swap + lower_b_swap)//2
        return not_solved, model
    

    def _optimize_circuit_normal(self, tight_bound_depth, lsqc, pi, time, sigma, count_gate, bound_depth, use_sabre):
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
        bit_tight_bound_depth = None
        n_swap = -1
        step = 5
        while not find_min_depth:
            bit_tight_bound_depth = BitVecVal(tight_bound_depth-1, length)
            print("Trying maximal depth = {}...".format(tight_bound_depth))
            start_time = datetime.datetime.now()
            # for depth optimization
            satisfiable = lsqc.check([UGE(bit_tight_bound_depth, time[l]) for l in range(count_gate)])
            print("status:{}, Depth optimization time = {}, time including preprocessing = {}".format(satisfiable, datetime.datetime.now() - start_time, timeit.default_timer()-self.start))
            if satisfiable == sat:
                model = lsqc.model()
                n_swap = self._count_swap(model, sigma, tight_bound_depth)
                upper_b_swap = min(n_swap-1, upper_b_swap)
                bound_swap_num = upper_b_swap
                # print("Find minimal depth {} with swap num {}".format(tight_bound_depth, model[count_swap].as_long()))
                for i in range(1, step):
                    tight_bound_depth -= 1
                    bit_tight_bound_depth = BitVecVal(tight_bound_depth-1, length)
                    satisfiable = lsqc.check([UGE(bit_tight_bound_depth, time[l]) for l in range(count_gate)])
                    if satisfiable == sat:
                        model = lsqc.model()
                        n_swap = self._count_swap(model, sigma, tight_bound_depth-i)
                        upper_b_swap = min(n_swap-1, upper_b_swap)
                        bound_swap_num = upper_b_swap
                    else:
                        bit_tight_bound_depth += 1
                        break
                find_min_depth = True
                print("Find minimal depth {} with swap num {}".format(tight_bound_depth, n_swap))
            else:
                # lsqc.pop()
                if tight_bound_depth > 100: 
                    step = 10
                tight_bound_depth = step + tight_bound_depth
                if tight_bound_depth > bound_depth:
                    print("FAIL to find solution with depth less than  {}.".format(tight_bound_depth-step))
                    break
        if not find_min_depth:
            return True, None
        if not self.obj_is_swap:
            return False, model
        # stats = lsqc.statistics()
        # print(stats)
        lsqc.add([UGE(bit_tight_bound_depth, time[l]) for l in range(count_gate)])

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
            satisfiable = lsqc.check()
            # satisfiable = lsqc.check(PbLe([(sigma[k][t],1) for k in range(count_qubit_edge)
            #              for t in range(bound_depth)], bound_swap_num) )
            # lsqc.pop()
            print("status:{}, optimization time = {}, time including preprocessing = {}".format(satisfiable, datetime.datetime.now() - start_time, timeit.default_timer()-self.start))
            if satisfiable == sat:
                model = lsqc.model()
                # cur_swap = model[count_swap].as_long()
                cur_swap = self._count_swap(model, sigma, tight_bound_depth)
                result = self._extract_results(model, time, pi, sigma, "qasm", "intermediate_result_swap_count_{}.qasm".format(cur_swap))
                # print(cur_swap, lower_b_swap)
                if cur_swap > lower_b_swap:
                    upper_b_swap = cur_swap
                    bound_swap_num = upper_b_swap - 1
                    # if use_sabre:
                    #     bound_swap_num = upper_b_swap - 1
                    # else:
                    #     bound_swap_num = (upper_b_swap + lower_b_swap) // 2
                else: 
                    find_min_swap = True
                    not_solved = False
                lsqc.pop()
            else:
                lsqc.pop()
                lower_b_swap = bound_swap_num + 1
                if upper_b_swap <= lower_b_swap:
                    while not find_min_swap:
                        start_time = datetime.datetime.now()
                        print("Trying min swap = {}...".format(lower_b_swap))
                        # TODO: add atmost-k constraint
                        lsqc.push()
                        self._add_atmostk_cnf(lsqc, sigma, lower_b_swap, tight_bound_depth)
                        satisfiable = lsqc.check()
                        # satisfiable = lsqc.check(PbLe([(sigma[k][t],1) for k in range(count_qubit_edge)
                        #      for t in range(bound_depth)], upper_b_swap) ):;q
                        print("status:{}, optimization time = {}, time including preprocessing = {}".format(satisfiable, datetime.datetime.now() - start_time, timeit.default_timer()-self.start))
                        if satisfiable == sat :
                            model = lsqc.model()
                            find_min_swap = True
                            not_solved = False
                        else:
                            lower_b_swap += 1
                        lsqc.pop()
                elif not use_sabre:
                    bound_swap_num = (upper_b_swap + lower_b_swap)//2
        return not_solved, model

    
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
            result_time.append(model[time[l]].as_long())
            result_depth = max(result_depth, result_time[-1])
        result_depth += 1
        list_result_swap = []
        # print(m.evaluate(f(x)))
        for k in range(count_qubit_edge):
            for t in range(result_depth):
                if model[sigma[k][t]]:
                    list_result_swap.append((k, t))
                    print(f"SWAP on physical edge ({list_qubit_edge[k][0]},"\
                        + f"{list_qubit_edge[k][1]}) at time {t}")
        for l in range(count_gate):
            if len(list_gate_qubits[l]) == 1:
                qq = list_gate_qubits[l][0]
                tt = result_time[l]
                # print(f"Gate {l}: {list_gate_name[l]} {qq} on qubit "\
                #     + f"{model[pi[qq][tt]].as_long()} at time {tt}")
                print(f"Gate {l}: {list_gate_name[l]} {qq} on qubit "\
                    + f"{model.evaluate(pi[qq][tt])} at time {tt}")
            else:
                qq = list_gate_qubits[l][0]
                qqq = list_gate_qubits[l][1]
                tt = result_time[l]
                print(f"Gate {l}: {list_gate_name[l]} {qq}, {qqq} on qubits "\
                    + f"{model[pi[qq][tt]].as_long()} and "\
                    + f"{model[pi[qqq][tt]].as_long()} at time {tt}")
                
                # print(f"Gate {l}: {list_gate_name[l]} {qq}, {qqq} on qubits "\
                #     + f"{model.evaluate(f_map(qq,tt))} and "\
                #     + f"{model.evaluate(f_map(qqq,tt))} at time {tt}")
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
                            p0 = model[pi[qubits[0]][block]].as_long()
                            # p0 = model[f_map(qubits[0],block)].as_long()
                            real_time[tmp_gate] = \
                                list_depth_on_qubit[p0] + 1
                            list_depth_on_qubit[p0] = \
                                real_time[tmp_gate]
                        else:
                            p0 = model[pi[qubits[0]][block]].as_long()
                            p1 = model[pi[qubits[1]][block]].as_long()
                            # p0 = model.evaluate(f_map(qubits[0],block))
                            # p1 = model.evalutae(f_map(qubits[1],block))
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
            result_time = real_time
            # print(list_result_swap)
        else:
            real_time = result_depth

        print(f"result- additional SWAP count = {len(list_result_swap)}.")
        print(f"result- circuit depth = {result_depth}.")

        list_scheduled_gate_qubits = [[] for i in range(result_depth)]
        list_scheduled_gate_name = [[] for i in range(result_depth)]
        for l in range(count_gate):
            t = result_time[l]
            list_scheduled_gate_name[t].append(list_gate_name[l])
            if l in list_gate_single:
                q = q = model[pi[list_gate_qubits[l][0]][model[time[l]].as_long()]].as_long()
                list_scheduled_gate_qubits[t].append((q,))
            elif l in list_gate_two:
                [q0, q1] = list_gate_qubits[l]
                tmp_t = t
                if self.mode == Mode.transition:
                    tmp_t = model[time[l]].as_long()
                q0 = model[pi[q0][tmp_t]].as_long()
                q1 = model[pi[q1][tmp_t]].as_long()
                # q0 = model.evaluate(f_map(q0,tmp_t)).as_long()
                # q1 = model.evaluate(f_map(q1,tmp_t)).as_long()
                list_scheduled_gate_qubits[t].append((q0, q1))
            else:
                raise ValueError("Expect single-qubit or two-qubit gate.")

        tmp_depth = result_depth - 1
        if self.mode == Mode.transition:
            tmp_depth = tran_detph - 1
        final_mapping = [model[pi[m][tmp_depth]].as_long() for m in range(count_program_qubit)]
        initial_mapping = [model[pi[m][0]].as_long() for m in range(count_program_qubit)]

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

        if self.obj_is_swap:
            objective_value = len(list_result_swap)    
        objective_value = result_depth
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
                    initial_mapping,
                    objective_value)
        else:
            return (output_qasm(self.device, result_depth, list_scheduled_gate_name,
                                list_scheduled_gate_qubits, final_mapping,
                                True, output_file_name),
                    final_mapping,
                    initial_mapping)

    def get_swap_upper_bound(self, heuristic = "sabre"):
        if heuristic == "sabre":
            swap_num, depth = run_sabre(self.list_gate_qubits, self.list_qubit_edge, self.count_physical_qubit)
            print("Run heuristic compiler sabre to get upper bound for SWAP: {}, depth: {}".format(swap_num, depth))
        else:
            raise TypeError("Only support sabre.")
        return swap_num, depth
