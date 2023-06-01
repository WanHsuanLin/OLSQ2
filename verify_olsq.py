from olsq.input import input_qasm
from olsq.device import qcdevice
# from collections import defaultdict

# def verifier(input_circ, input_mode, output_circ, output_mode, device):
#def verifier(input_circ, output_circ, device = None, input_mode = None, output_mode = None):
class Verifier:
    def verify_result(self, input_circ, output_circ, device:qcdevice, all_commute: bool = False):
    # parse the input program
        input_mode = 'IR' if isinstance(input_circ, list) else 'QASM'
        output_mode = 'IR' if len(output_circ) == 6 else 'QASM'

        if input_mode == 'IR':
            count_logic_qubit = input_circ[0]
            list_logic_gate_qubits = input_circ[1]
            list_logic_gate_name = input_circ[2]
        else:
            input_program = input_qasm(input_circ)
            count_logic_qubit = input_program[0]
            list_logic_gate_qubits = input_program[1]
            list_logic_gate_name = input_program[2]
        # print(list_logic_gate_name)
        # print(list_logic_gate_qubits)

        # parse the output program
        if output_mode == 'IR':
            result_depth = output_circ[0]
            list_scheduled_gate_name = output_circ[1]
            list_scheduled_gate_qubits = output_circ[2]
            final_mapping = output_circ[3]
            initial_mapping = output_circ[4]
        elif output_mode == 'QASM':
            output_program = input_qasm(output_circ[0])
            final_mapping = output_circ[1]
            initial_mapping = output_circ[2]
        # [16, 
        # ([5, 4], [15, 14], [5, 9], [15, 11], [14, 13], [11, 7], [1, 5], [12, 13], [1, 0], [14, 13], [5, 9], [3, 7], [4, 0], [9, 13], [11, 7], [14, 15], [4, 8], [14, 10], [13, 12], [3, 7], [0, 1], [5, 9], [9, 8], [1, 2], [7, 6], [9, 10], [12, 8], [3, 2], [6, 10], [6, 2]), 
        # ('rzz(pi/4)', 'rzz(pi/4)', 'rzz(pi/4)', 'rzz(pi/4)', 'rzz(pi/4)', 'rzz(pi/4)', 'SWAP', 'SWAP', 'rzz(pi/4)', 'rzz(pi/4)', 'rzz(pi/4)', 'SWAP', 'rzz(pi/4)', 'rzz(pi/4)', 'rzz(pi/4)', 'SWAP', 'rzz(pi/4)', 'rzz(pi/4)', 'rzz(pi/4)', 'rzz(pi/4)', 'SWAP', 'SWAP', 'rzz(pi/4)', 'rzz(pi/4)', 'rzz(pi/4)', 'rzz(pi/4)', 'rzz(pi/4)', 'rzz(pi/4)', 'rzz(pi/4)', 'rzz(pi/4)')]
        # print(list_scheduled_gate_name)
        # print(list_scheduled_gate_qubits)

        qubits_to_gate = dict()
        # map the logic qubit pair to their gate number
        for i, n in enumerate(list_logic_gate_qubits):
            qubits_to_gate[tuple(n)] = i
        # print(qubits_to_gate)
        # condense 2d gate arrays into 1d for traversing
        result_gate_names = []
        if output_mode == 'IR':
            for i in range(len(list_scheduled_gate_name)):
                cur = list_scheduled_gate_name[i]
                for j in range(len(cur)):
                    result_gate_names.append(cur[j])
        # print(result_gate_names)
        elif output_mode == 'QASM':
            result_gate_names = list(output_program[2])
        
        result_gate_qubits = []
        if output_mode == 'IR':
            for i in range(len(list_scheduled_gate_qubits)):
                cur = list_scheduled_gate_qubits[i]
                for j in range(len(cur)):
                    result_gate_qubits.append(cur[j])
        elif output_mode == 'QASM':
            result_gate_qubits = list(output_program[1])
        # print(result_gate_qubits)

        # check some constraints
        if len(result_gate_names) != len(result_gate_qubits):
            return False
        if len(final_mapping) != len(initial_mapping) or count_logic_qubit != len(initial_mapping):
            return False
        # map the logical qubits to physical qubits using initial mapping
        mapping = dict()
        for i, n in enumerate(initial_mapping):
            mapping[n] = i
        # print(mapping)

        # define the connection edges
        connections = device.list_qubit_edge
        # print(edges)
        # connections = set()
        # for n in edges:
        #     connections.add(n)

        # construct the dependency graph based on the input program
        list_collision = list()
        for g in range(len(list_logic_gate_qubits)):
            for gg in range(g + 1, len(list_logic_gate_qubits)):
                if list_logic_gate_qubits[g][0] == list_logic_gate_qubits[gg][0]:
                        list_collision.append((g, gg))
                    
                if len(list_logic_gate_qubits[gg]) == 2:
                    if list_logic_gate_qubits[g][0] == list_logic_gate_qubits[gg][1]:
                        list_collision.append((g, gg))
                
                if len(list_logic_gate_qubits[g]) == 2:
                    if list_logic_gate_qubits[g][1] == list_logic_gate_qubits[gg][0]:
                        list_collision.append((g, gg))
                    if len(list_logic_gate_qubits[gg]) == 2:
                        if list_logic_gate_qubits[g][1] == list_logic_gate_qubits[gg][1]:
                            list_collision.append((g, gg))
        # print(list_collision)
        # dependency = defaultdict(set)
        # for para, dep in list_collision:
        #         dependency[dep].add(para)
        # print(dependency)
        # save all gates' parent and child gates
        # only look at the perspective gate's child gates to remove it from its parent gate
        # [()] list of set = gate_parents - parent gates of the current gate, [()] list of set = gate_children - children gates of the current gate
        gate_parents = [set () for _ in range(len(list_logic_gate_name))]
        gate_children = [set () for _ in range(len(list_logic_gate_name))]
        for pairs in list_collision:
            child, parent = pairs[0], pairs[1]
            gate_parents[parent].add(child)
            gate_children[child].add(parent)
        # print(gate_children)
        # print(gate_parents)

        # if all_commute is true, no need to check dependency
        # traverse the physical gates
        for i in range(len(result_gate_names)):
            if result_gate_names[i] == "SWAP":
                if len(result_gate_qubits[i]) != 2:
                    return False
                swap_qb1 = result_gate_qubits[i][0]
                swap_qb2 = result_gate_qubits[i][1]
                temp = mapping[swap_qb1]
                mapping[swap_qb1] = mapping[swap_qb2]
                mapping[swap_qb2] = temp
            else:
                if len(result_gate_qubits[i]) == 1:
                    logic_qubit = mapping[result_gate_qubits[i][0]]
                    # print((logic_qubit, ))
                    # print(qubits_to_gate)
                    # return False
                    gate_num = qubits_to_gate[(logic_qubit,)]
                    if gate_parents[gate_num]:
                        return False
                    # remove this element from rest of dict values and remove the dict items with 0 values from the dict
                    for n in gate_children[gate_num]:
                        gate_parents[n].remove(gate_num)
                    gate_children[gate_num].clear()

                elif len(result_gate_qubits[i]) == 2:
                    physical_qb1 = result_gate_qubits[i][0]
                    physical_qb2 = result_gate_qubits[i][1]
                    logic_qb1 = mapping[physical_qb1]
                    logic_qb2 = mapping[physical_qb2]
                    if (physical_qb1, physical_qb2) not in connections and (physical_qb2, physical_qb1) not in connections: # wrong, change it
                        return False
                    gate_num = qubits_to_gate[(logic_qb1, logic_qb2)]
                    if gate_parents[gate_num]:
                        return False
                    # remove this element from rest of dict values and remove the dict items with 0 values from the dict
                    for n in gate_children[gate_num]:
                        gate_parents[n].remove(gate_num)
                    gate_children[gate_num].clear()

        end_mapping = [-999 for i in range(count_logic_qubit)]
        for key, value in mapping.items():
            end_mapping[value] = key
        # print(end_mapping)
        # return True
        return True if end_mapping == final_mapping else False



# test information
# output_str = [10, [["rzz(pi/4)", "rzz(pi/4)"], ["rzz(pi/4)", "rzz(pi/4)", "rzz(pi/4)"], ["rzz(pi/4)", "SWAP", "SWAP"], 
#       ["rzz(pi/4)", "rzz(pi/4)", "rzz(pi/4)", "SWAP"], ["rzz(pi/4)", "rzz(pi/4)", "rzz(pi/4)"], 
#       ["rzz(pi/4)", "rzz(pi/4)", "rzz(pi/4)", "SWAP", "SWAP", "SWAP"], 
#       ["rzz(pi/4)", "rzz(pi/4)", "rzz(pi/4)", "rzz(pi/4)"], ["rzz(pi/4)", "rzz(pi/4)", "rzz(pi/4)"], ["rzz(pi/4)"], 
#       ["rzz(pi/4)"]],[[[9, 8], [3, 2]], [[9, 5], [3, 7], [2, 1]], [[7, 11], [0, 1], [9, 13]], 
#                       [[13, 12], [2, 1], [9, 5], [11, 15]], [[8, 12], [5, 1], [7, 11]], [[8, 4], [1, 0], [15, 11], 
#                       [3, 7], [5, 9], [12, 13]], [[7, 6], [5, 4], [13, 14], [11, 10]], [[5, 6], [0, 4], [15, 14]], 
#                       [[10, 6]], [[10, 14]]], [7, 1, 5, 10, 12, 2, 11, 4, 14, 0, 13, 9, 6, 3, 15, 8], 
#                       [3, 0, 13, 10, 9, 2, 15, 4, 14, 1, 12, 5, 6, 7, 11, 8]]

# [16, 
# ([5, 4], [15, 14], [5, 9], [15, 11], [14, 13], [11, 7], [1, 5], [12, 13], [1, 0], [14, 13], [5, 9], [3, 7], [4, 0], [9, 13], [11, 7], [14, 15], [4, 8], [14, 10], [13, 12], [3, 7], [0, 1], [5, 9], [9, 8], [1, 2], [7, 6], [9, 10], [12, 8], [3, 2], [6, 10], [6, 2]), 
# ('rzz(pi/4)', 'rzz(pi/4)', 'rzz(pi/4)', 'rzz(pi/4)', 'rzz(pi/4)', 'rzz(pi/4)', 'SWAP', 'SWAP', 'rzz(pi/4)', 'rzz(pi/4)', 'rzz(pi/4)', 'SWAP', 'rzz(pi/4)', 'rzz(pi/4)', 'rzz(pi/4)', 'SWAP', 'rzz(pi/4)', 'rzz(pi/4)', 'rzz(pi/4)', 'rzz(pi/4)', 'SWAP', 'SWAP', 'rzz(pi/4)', 'rzz(pi/4)', 'rzz(pi/4)', 'rzz(pi/4)', 'rzz(pi/4)', 'rzz(pi/4)', 'rzz(pi/4)', 'rzz(pi/4)')]

# def get_nnGrid(n: int, swap_duration):
#     my_coupling = []
#     for i in range(n):
#         for j in range(n):
#             if j < n-1:
#                 my_coupling.append((i*n +j,i*n +j+1))
#             if i < n-1:
#                 my_coupling.append((i*n +j,i*n +j+n))
   
#     device = qcdevice(name="grid", nqubits=n*n,
#         connection=my_coupling, swap_duration=swap_duration)
#     return device
#     # print(device.list_qubit_edge)

# input_circuit_file = open("/Users/xinlin/Documents/OLSQ_Project/OLSQ2/benchmark/qaoa/qaoa_16_0.qasm", "r").read()
# output_circuit_file = open("/Users/xinlin/Documents/OLSQ_Project/OLSQ2/example/grid_qaoa_16_0.qasm", "r").read()
# final = [14, 13, 9, 6, 0, 15, 7, 8, 2, 12, 1, 5, 10, 11, 3, 4]
# initial = [15, 12, 1, 6, 5, 14, 3, 8, 2, 13, 0, 9, 10, 11, 7, 4]
# output_circuit = [output_circuit_file, final, initial]

# verifier = Verifier()
# print(verifier.verify_result(input_circuit_file, output_circuit, get_nnGrid(4, 1)))

# print(verifier(circuit_file, output_str, get_nnGrid(4, 1)))
# # circuit_file = open("/Users/xinlin/Documents/OLSQ_Project/OLSQ2/benchmark/others/tof_4.qasm", "r").read()
# # verifier(circuit_str, 'default', output_str, 'IR')

# # input program qubit: [gate, gate, gate]
# # output program qubit: [gate, gate, gate] (count the gates in terms of program qubits)
# # 0 -> 9 x , h
# 0 -> 8: final mapping 
