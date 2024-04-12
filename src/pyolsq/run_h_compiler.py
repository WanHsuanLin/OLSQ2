
from qiskit.transpiler import CouplingMap, Layout
from qiskit import QuantumCircuit
from qiskit.transpiler import PassManager
from qiskit.transpiler.passes import SabreLayout, SabreSwap

def run_sabre(is_qasm, circuit_info, coupling, count_physical_qubit):
    # read qasm
    if not is_qasm:
        list_gate = circuit_info
        qc = QuantumCircuit(count_physical_qubit)
        for gate in list_gate:
            if len(gate) == 2:
                qc.cx(gate[0], gate[1])
            elif len(gate) == 1:
                qc.h(gate[0])
            else:
                raise TypeError("Currently only support one and two-qubit gate.")
    else:
        # print(circuit_info)
        # print("+____________+")
        qc = QuantumCircuit.from_qasm_str(circuit_info)
        # print(qc.qasm())
    qubit2gateId = dict()
    idx = 0
    for gate in qc.data:
        qubit_list = []
        for q in gate[1]:
            qubit_list.append(q.index)
        qubit_list.sort()
        qubit2gateId[tuple(qubit_list)] = idx
        idx += 1
    # qc.draw(scale=0.7, filename = "ori.png", output='mpl', style='color')
    device = CouplingMap(couplinglist = coupling, description="sabre_test")
    # print("qubit2gateId")
    # print(qubit2gateId)
    # initialize sabre
    # ["basic", "lookahead", "decay"]
    sbs = SabreSwap(coupling_map = device, heuristic = "lookahead", seed = 0)
    sbl = SabreLayout(coupling_map = device, seed = 0)
    pass_manager1 = PassManager(sbl)
    sabre_cir = pass_manager1.run(qc)
    # sabre_cir.draw(scale=0.7, filename="sabrecir1.png", output='mpl', style='color')
    pass_manager2 = PassManager(sbs)
    # print(sabre_cir._layout)
    # print(sabre_cir._layout.initial_layout.get_physical_bits())
    # print(sabre_cir._layout.initial_layout.get_virtual_bits())
    initial_mapping = [0] * count_physical_qubit
    # current_mapping_dict = dict() # physical qubit to logical qubit
    initial_mapping_sabre = sabre_cir._layout.initial_layout.get_virtual_bits()
    i = 0
    for q in initial_mapping_sabre:
        # print(q, ", ", initial_mapping_sabre[q])
        # current_mapping_dict[i] = initial_mapping_sabre[q]
        initial_mapping[initial_mapping_sabre[q]] = i
        i += 1
    sabre_cir = pass_manager2.run(sabre_cir)
    # print(sabre_cir.qasm())
    # sabre_cir.draw(scale=0.7, filename="sabrecir2.png", output='mpl', style='color')
    # print("initial_mapping:")
    # print(initial_mapping)
    # print(initial_mapping_sabre)
    # print(current_mapping_dict)
    count_swap = 0
    qubit_last_gate_time = [-1] * sabre_cir.num_qubits
    pyPlan = [[]]
    current_mapping = [i for i in range(count_physical_qubit)]
    # print("Print SABRE's results")
    for gate in sabre_cir._data:
        # print(gate[0].name)
        # print(gate[0].num_qubits)
        # print(gate[1])
        if gate[0].name == 'swap':
            count_swap += 1
            tmp = current_mapping[gate[1][0].index]
            current_mapping[gate[1][0].index] = current_mapping[gate[1][1].index]
            current_mapping[gate[1][1].index] = tmp
            # print("current mapping:")
            # print(current_mapping)
        max_time = -1
        for q in gate[1]:
            max_time = max(max_time, qubit_last_gate_time[q.index])
        max_time += 1
        qubit_list = []
        for q in gate[1]:
            qubit_last_gate_time[q.index] = max_time
            qubit_list.append(current_mapping[q.index])
        if max_time >= len(pyPlan):
            pyPlan.append([])
        if gate[0].name == 'swap':
            pyPlan[max_time].append(-1)
        else:
            qubit_list.sort()
            # print("qubit_list:")
            # print(qubit_list)
            pyPlan[max_time].append(qubit2gateId[tuple(qubit_list)])
            # pyPlan[max_time].append(qubit2gateId[tuple(gate[1])])
        # print("plan:")
        # print(pyPlan)
    
    return count_swap, sabre_cir.depth(), initial_mapping, pyPlan