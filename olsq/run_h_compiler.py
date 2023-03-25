from qiskit.transpiler import CouplingMap, Layout
from qiskit import QuantumCircuit
from qiskit.transpiler import PassManager
from qiskit.transpiler.passes import SabreLayout, SabreSwap
from qiskit.converters import *
from qiskit.transpiler.passes import Unroller

def run_sabre(benchmark, circuit_info, coupling, count_physical_qubit):
    # read qasm
    if benchmark == "olsq":
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
        qc = QuantumCircuit.from_qasm_str(circuit_info)
    # qc.draw(scale=0.7, filename = "cir_for_tket.png", output='mpl', style='color')
    device = CouplingMap(couplinglist = coupling, description="sabre_test")
    
    # initialize sabre
    # ["basic", "lookahead", "decay"]
    sbs = SabreSwap(coupling_map = device, heuristic = "lookahead", seed = 0)
    sbl = SabreLayout(coupling_map = device, seed = 0)
    pass_manager1 = PassManager(sbl)
    sabre_cir = pass_manager1.run(qc)
    # sabre_cir.draw(scale=0.7, filename="sabrecir1.png", output='mpl', style='color')
    pass_manager2 = PassManager(sbs)
    # print(sabre_cir._layout)
    # print(sabre_cir._layout.get_physical_bits())
    # print(sabre_cir._layout.get_virtual_bits())
    initial_mapping = []
    initial_mapping_sabre = sabre_cir._layout.get_virtual_bits()
    for q in initial_mapping_sabre:
        initial_mapping.append(initial_mapping_sabre[q])
    sabre_cir = pass_manager2.run(sabre_cir)
    # sabre_cir.draw(scale=0.7, filename="sabrecir2.png", output='mpl', style='color')
    
    count_swap = 0
    for gate in sabre_cir.data:
        # print(gate[0].name)
        # print(gate[0].num_qubits)
        if gate[0].name == 'swap':
            count_swap += 1
    
    if benchmark != "qaoa" and benchmark != "olsq":
        unroller = Unroller(basis=['x', 'u', 'cx', 'rz', 't', 'h', 'tdg'])
        sabre_cir = unroller.run(circuit_to_dag(sabre_cir))


    return count_swap, sabre_cir.depth(), initial_mapping
