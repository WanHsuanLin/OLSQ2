from qiskit.transpiler import CouplingMap, Layout
from qiskit import QuantumCircuit
from qiskit.transpiler import PassManager
from qiskit.transpiler.passes import SabreLayout, SabreSwap
from qiskit.converters import *
from qiskit.transpiler.passes import Unroller

def run_sabre(circuit_info, coupling, count_physical_qubit):
    # read qasm
    list_gate = circuit_info
    qc = QuantumCircuit(count_physical_qubit)
    for gate in list_gate:
        if len(gate) == 2:
            qc.cx(gate[0], gate[1])
        elif len(gate) == 1:
            qc.h(gate[0])
        else:
            raise TypeError("Currently only support one and two-qubit gate.")
    
    # qc.draw(scale=0.7, filename = "cir_for_tket.png", output='mpl', style='color')
    device = CouplingMap(couplinglist = coupling, description="sabre_test")
    
    # initialize sabre
    # ["basic", "lookahead", "decay"]
    sbs = SabreSwap(coupling_map = device, heuristic = "lookahead", seed = 0)
    sbl = SabreLayout(coupling_map = device, seed = 0)
    pass_manager1 = PassManager(sbl)
    sabre_cir = pass_manager1.run(qc)
    pass_manager2 = PassManager(sbs)
    sabre_cir = pass_manager2.run(sabre_cir)
    # sabre_cir.draw(scale=0.7, filename="sabrecir2.png", output='mpl', style='color')
    
    count_swap = 0
    for gate in sabre_cir.data:
        # print(gate[0].name)
        # print(gate[0].num_qubits)
        if gate[0].name == 'swap':
            count_swap += 1

    return count_swap, sabre_cir.depth()
