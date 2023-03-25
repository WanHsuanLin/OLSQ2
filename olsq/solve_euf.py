import math
import datetime

from z3 import Bool, Implies, And, Or, sat, unsat, Solver, set_option, BitVec, ULT, ULE, UGE, BitVecVal, PbLe, Not, Function, BitVecSort, Distinct, BoolSort

from olsq.input import input_qasm
from olsq.output import output_qasm
from olsq.device import qcdevice
import pkgutil
from enum import Enum

TIMEOUT = 90000
# MEMORY_MAX_SIZE = 1000000000 * 58
MEMORY_MAX_SIZE = 0
MAX_TREAD_NUM = 8
VERBOSE = 3

class Mode(Enum):
    transition = 1
    normal = 2
    mix = 3

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
        elif mode == "mixed":
            self.mode = Mode.mix
        else:
            raise ValueError("Invalid Choice of Transition-Based or Not")

        # These values should be updated in setdevice(...)
        self.device = None
        self.count_physical_qubit = 0
        self.list_qubit_edge = []
        self.swap_duration = 0
        self.dict_gate_duration = dict()
        self.list_gate_duration = []

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
        if gate_duration != None:
            # self.list_gate_duration
            for gate_name in self.list_gate_name:
                self.list_gate_duration.append(self.dict_gate_duration[gate_name])
        else:
            self.list_gate_duration = [1]*len(self.list_gate_qubits)

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
        # bound_depth: generate constraints until t = bound_depth
        bound_depth = 21
        print("start adding constraints...")
        # variable setting 
        output_file_name = folder+"/"+str(self.count_physical_qubit)+"_"+str(self.count_program_qubit) + "_21" + ".txt"
        self.list_gate_dependency = collision_extracting(self.list_gate_qubits)
        f_map, f_map_inv, time, sigma = self._construct_variable(bound_depth)
        lsqc = Solver()
        self._add_injective_mapping_constraints(bound_depth, f_map, f_map_inv, lsqc)
        self._add_consistency_gate_constraints(f_map, time, lsqc)
        self._add_dependency_constraints(False, lsqc, time, bound_depth)
        self._add_swap_constraints(bound_depth, sigma, lsqc, True, time, f_map)
        self._add_transformation_constraints(bound_depth, lsqc, sigma, f_map)
        lsqc.add([UGE(bound_depth, time[l]) for l in range(len(self.list_gate_qubits))])
        constraints = lsqc.sexpr()
        f = open(output_file_name, "w")
        f.write(constraints)
        f.write("(check-sat)")
        f.close()
        return

    def solve(self, output_mode: str = None, output_file_name: str = None, memory_max_size=MEMORY_MAX_SIZE, verbose = VERBOSE):
        """Formulate an SMT, pass it to z3 solver, and output results.
        CORE OF OLSQ, EDIT WITH CARE.

        Args:
            preprossess_only: Only used to find the bound for SWAP
        
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
        if not self.input_dependency:
            self.list_gate_dependency = collision_extracting(self.list_gate_qubits)
        if self.mode == Mode.transition:
            print("Using transition based mode...")
            _, results = self._solve(False, None, output_mode, output_file_name, memory_max_size, verbose)
        elif self.mode == Mode.normal:
            print("Using normal mode...")
            _, results = self._solve(False, None, output_mode, output_file_name, memory_max_size, verbose)
        elif self.mode == Mode.mix:
            print("Using mixed mode...")
            print("Perform preprocessing to find swap bound...")
            swap_bound, _ = self._solve(True, None, None, memory_max_size, verbose)
            print(f"Finish proprocessing. SWAP bound is ({swap_bound[0]},{swap_bound[1]})")
            print("Start normal searching...")
            _, results = self._solve(False, swap_bound, output_mode, output_file_name, memory_max_size, verbose)
        else:
            raise ValueError( ("Wrong type") )
        return results 

    def _solve(self, preprossess_only, swap_bound, output_mode: str = None, output_file_name: str = None, memory_max_size=MEMORY_MAX_SIZE, verbose = VERBOSE):
        

        # pre-processing
        count_gate = len(self.list_gate_qubits)
        

        not_solved = True
        model = None
        iteration = 1

        # tight_bound_depth: use for depth constraint

        # bound_depth: generate constraints until t = bound_depth
        if preprossess_only or self.mode == Mode.transition:
            bound_depth = 8 * self.bound_depth
        else:
            bound_depth = 2 * self.bound_depth
        while not_solved:
            print("start adding constraints...")
            # variable setting 
            f_map, time, sigma = self._construct_variable(bound_depth)
            lsqc = Solver()
            set_option("memory_max_size", memory_max_size)
            set_option("verbose", verbose)

            # constraint setting
            self._add_injective_mapping_constraints(bound_depth, f_map, lsqc)

            # Consistency between Mapping and Space Coordinates.
            self._add_consistency_gate_constraints(f_map, time, lsqc)
            
            # Avoiding Collisions and Respecting Dependencies. 
            self._add_dependency_constraints(preprossess_only, lsqc, time, bound_depth)

            # # No swap for t<s
            # # swap gates can not overlap with swap
            if preprossess_only or self.mode == Mode.transition:
                self._add_swap_constraints(bound_depth, sigma, lsqc)
            else:
                self._add_swap_constraints(bound_depth, sigma, lsqc, True, time, f_map)
            # Mapping Not Transformations by SWAP Gates.
            # Mapping Transformations by SWAP Gates.
            self._add_transformation_constraints(bound_depth, lsqc, sigma, f_map)
            start_time = datetime.datetime.now()
            # TODO: iterate each swap num
            not_solved, model, n_swap = self._optimize_circuit(iteration, lsqc, preprossess_only, time, sigma, count_gate, bound_depth, swap_bound)
            if not_solved:
                bound_depth *= 2
                iteration += 1
            else:
                if preprossess_only:
                    swap_bound = (self.bound_depth-1 , n_swap)
                if swap_bound != None:
                    swap_bound = (swap_bound[0], n_swap)
                else:
                    swap_bound = (0, n_swap)
            result = self._extract_results(model, time, f_map, sigma, output_mode, output_file_name)
                # result = self._extract_results(model, time, f_map, sigma, space, output_mode, output_file_name)
                
        print(f"compilation time = {datetime.datetime.now() - start_time}.")
        return swap_bound, result

    def _construct_variable(self, bound_depth):
        # at cycle t, logical qubit q is mapped to pi[q][t]
        count_qubit_edge = len(self.list_qubit_edge)
        length = int(math.log2(bound_depth))+1
        
        num_phy_q = int(math.log2(self.count_physical_qubit))+1
        num_pro_q = int(math.log2(self.count_program_qubit))+1
        f_map = Function('mapping', BitVecSort(num_pro_q), BitVecSort(length), BitVecSort(num_phy_q))
        f_map_inv = Function('mapping_inv', BitVecSort(num_phy_q), BitVecSort(length), BitVecSort(num_pro_q))

        # time coordinate for gate l is time[l]
        time = [BitVec("time_{}".format(i), length) for i in range(len(self.list_gate_qubits))]

        # if at cycle t, a SWAP finishing on edge k, then sigma[k][t]=1
        sigma = [[Bool("ifswap_e{}_t{}".format(i, j))
            for j in range(bound_depth)] for i in range(count_qubit_edge)]

        return f_map, f_map_inv, time, sigma


    def _add_injective_mapping_constraints(self, bound_depth, f_map, f_map_inv, model):
        # Injective Mapping
        for t in range(bound_depth):
            # model.add(Distinct([f_map(m,t) for m in range(self.count_program_qubit)]))
            for m in range(self.count_program_qubit):
                model.add(UGE(f_map(m,t), 0))
                model.add(ULT(f_map(m,t), self.count_physical_qubit))
            for t in range(bound_depth):
                for m in range(self.count_program_qubit):
                    model.add(f_map_inv(f_map(m,t),t) == m)

    def _add_consistency_gate_constraints(self, f_map, time, model):
        # Consistency between Mapping and Space Coordinates.
        list_gate_qubits = self.list_gate_qubits
        count_gate = len(self.list_gate_qubits)
        list_gate_two = self.list_gate_two

        for l in range(count_gate):
            if l in list_gate_two:
                model.add(Or([Or(And(f_map(list_gate_qubits[l][0],time[l]) == edge[0], f_map(list_gate_qubits[l][1],time[l]) == edge[1]), \
                                And(f_map(list_gate_qubits[l][0],time[l]) == edge[1], f_map(list_gate_qubits[l][1],time[l]) == edge[0])) for edge in self.list_qubit_edge ] ))

    def _add_dependency_constraints(self, preprossess_only, model, time, bound_depth):
        list_gate_duration = self.list_gate_duration
        list_gate_dependency = self.list_gate_dependency
        count_gate = len(self.list_gate_qubits)
        if preprossess_only or self.mode == Mode.transition:
            for d in list_gate_dependency:
                # lsqc.add(time[d[0]] <= time[d[1]])
                model.add(ULE(time[d[0]],time[d[1]]))
        else:
            length = int(math.log2(bound_depth))+1
            
            for d in list_gate_dependency:
                model.add(ULT(time[d[0]], time[d[1]]))
            # add initial condition for gates
            # for l in range(count_gate):
            #     model.add(ULE(bit_duration_minus_one_list[l], time[l]))

    def _add_swap_constraints(self, bound_depth, sigma, model, normal = False, time = None, f_map = None):
        # if_overlap_edge takes in two edge indices _e_ and _e'_,
        # and returns whether or not they overlap
        list_qubit_edge = self.list_qubit_edge
        count_qubit_edge = len(list_qubit_edge)
        if_overlap_edge = [[0] * count_qubit_edge
            for k in range(count_qubit_edge)]
        # list_over_lap_edge takes in an edge index _e_,
        # and returnsthe list of edges that overlap with _e_
        list_overlap_edge = list()
        # list_count_overlap_edge is the list of lengths of
        # overlap edge lists of all the _e_
        list_count_overlap_edge = list()
        list_gate_qubits = self.list_gate_qubits
        for k in range(count_qubit_edge):
            list_overlap_edge.append(list())
        for k in range(count_qubit_edge):
            for kk in range(k + 1, count_qubit_edge):
                if (   (list_qubit_edge[k][0] == list_qubit_edge[kk][0]
                        or list_qubit_edge[k][0] == list_qubit_edge[kk][1])
                    or (list_qubit_edge[k][1] == list_qubit_edge[kk][0]
                        or list_qubit_edge[k][1] == list_qubit_edge[kk][1]) ):
                    list_overlap_edge[k].append(kk)
                    list_overlap_edge[kk].append(k)
                    if_overlap_edge[kk][k] = 1
                    if_overlap_edge[k][kk] = 1
        for k in range(count_qubit_edge):
            list_count_overlap_edge.append(len(list_overlap_edge[k]))

        # No swap for t<s
        for t in range(min(self.swap_duration - 1, bound_depth)):
            for k in range(count_qubit_edge):
                model.add(sigma[k][t] == False)
        # swap gates can not overlap with swap
        for t in range(self.swap_duration - 1, bound_depth):
            for k in range(count_qubit_edge):
                for tt in range(t - self.swap_duration + 1, t):
                    model.add(Implies(sigma[k][t] == True,
                        sigma[k][tt] == False))
                for tt in range(t - self.swap_duration + 1, t + 1):
                    for kk in list_overlap_edge[k]:
                        model.add(Implies(sigma[k][t] == True,
                            sigma[kk][tt] == False))
        
        if normal:
            # swap gates can not ovelap with other gates
            # the time interval should be modified
            count_gate = len(self.list_gate_qubits)
            for t in range(self.swap_duration - 1, bound_depth):
                for k in range(count_qubit_edge):
                    for l in range(count_gate):
                        for tt in range(t - self.swap_duration + 1, t + self.list_gate_duration[l]):
                            if l in self.list_gate_single:
                                model.add(Implies(And(time[l] == tt,
                                    Or(f_map(list_gate_qubits[l][0],time[l]) == list_qubit_edge[k][0],
                                        f_map(list_gate_qubits[l][0],time[l]) == list_qubit_edge[k][1])),
                                        sigma[k][t] == False             ))
                            elif l in self.list_gate_two:
                                model.add(Implies(And( \
                                    time[l] == tt, \
                                        Or(f_map(list_gate_qubits[l][0],time[l]) == list_qubit_edge[k][0], \
                                        f_map(list_gate_qubits[l][0],time[l]) == list_qubit_edge[k][1], \
                                        f_map(list_gate_qubits[l][1],time[l]) == list_qubit_edge[k][0], \
                                        f_map(list_gate_qubits[l][1],time[l]) == list_qubit_edge[k][1])), \
                                        sigma[k][t] == False           ))
                                # for kk in list_overlap_edge[k]:
                                #     model.add(Implies(And(
                                #         time[l] == tt, f_map(list_gate_qubits[l][1],time[l]) == kk),
                                #             sigma[k][t] == False       ))          

     
    def _add_transformation_constraints(self, bound_depth, model, sigma, f_map):
        # list_adjacency_qubit takes in a physical qubit index _p_, and
        # returns the list of indices of physical qubits adjacent to _p_
        list_adjacent_qubit = list()
        # list_span_edge takes in a physical qubit index _p_,
        # and returns the list of edges spanned from _p_
        list_span_edge = list()
        list_qubit_edge = self.list_qubit_edge
        count_qubit_edge = len(list_qubit_edge)
        for n in range(self.count_physical_qubit):
            list_adjacent_qubit.append(list())
            list_span_edge.append(list())
        for k in range(count_qubit_edge):
            list_adjacent_qubit[list_qubit_edge[k][0]].append(
                                                        list_qubit_edge[k][1])
            list_adjacent_qubit[list_qubit_edge[k][1]].append(
                                                        list_qubit_edge[k][0])
            list_span_edge[list_qubit_edge[k][0]].append(k)
            list_span_edge[list_qubit_edge[k][1]].append(k)

        # Mapping Not Transformations by SWAP Gates.
        for t in range(bound_depth - 1):
            for n in range(self.count_physical_qubit):
                for m in range(self.count_program_qubit):
                    model.add(
                        Implies(And(Not(Or([sigma[k][t] for k in list_span_edge[n]])),
                                f_map(m,t) == n), f_map(m,t + 1) == n))
        
        # Mapping Transformations by SWAP Gates.
        for t in range(bound_depth - 1):
            for k in range(count_qubit_edge):
                for m in range(self.count_program_qubit):
                    model.add(Implies(And(sigma[k][t] == True,
                        f_map(m,t) == list_qubit_edge[k][0]),
                            f_map(m,t + 1) == list_qubit_edge[k][1]))
                    model.add(Implies(And(sigma[k][t] == True,
                        f_map(m,t) == list_qubit_edge[k][1]),
                            f_map(m,t + 1) == list_qubit_edge[k][0]))
           

    def _count_swap(self, model, sigma, result_depth):
        n_swap = 0
        for k in range(len(self.list_qubit_edge)):
            for t in range(result_depth):
                if model[sigma[k][t]]:
                    n_swap += 1
        return n_swap

  
