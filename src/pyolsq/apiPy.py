import os
import sys
sys.path.append(os.getcwd())
# sys.path.insert(0, '/Users/wanhsuan/Desktop/Github/OLSQ-dev/build')
sys.path.insert(0, '/home/wanhsuan/OLSQ-dev/build')
from src.pyolsq.input import input_qasm
from src.pyolsq.run_h_compiler import run_sabre
from olsqPy import Device, Circuit, OLSQ, addGate, setEdge, setDependency, setInitialMapping

def createCircuit(name, program, is_qasm= True, gate_duration: dict = None):
    """Translate input program to Circuit
    Args:
        name: circuit name
        program: a qasm string, or a list of the three things in IR.
        input_mode: (optional) can be "IR" if the input has ben
            translated to OLSQ IR; can be "benchmark" to use one of
            the benchmarks.  Default mode assumes qasm input.
    """
    
    if not is_qasm:
        count_program_qubit = program[0]
        list_gate_qubits = program[1]
        list_gate_name = program[2]
    else:
        program = input_qasm(program)
        count_program_qubit = program[0]
        list_gate_qubits = program[1]
        list_gate_name = program[2]

    circuit = Circuit(name, program[0], len(program[1]))
    if gate_duration != None:
        for g, n in zip(list_gate_qubits, list_gate_name):
            addGate(circuit, n, g, gate_duration[n])
    else:
        for g, n in zip(list_gate_qubits, list_gate_name):
            addGate(circuit, n, g)
    return circuit

def createDevice(name: str, nqubits: int = None, connection: list = None):
    """Pass in parameters from the given device.  If in TB mode,
        swap_duration is set to 1 without modifying the device.

    Args:
        device: a qcdevice object for OLSQ
    """
    device = Device(name, nqubits, len(connection))
    setEdge(device, connection)
    return device

def useSabre(for_swap: bool, for_mapping: bool, olsq: OLSQ, circuit:Circuit, nqubits, connection: list, program, is_qasm):
    swap_num, depth, initial_mapping = run_sabre(is_qasm, program, connection, nqubits)
    print("[Info] Run heuristic compiler SABRE to get upper bound for SWAP: {}, depth: {}".format(swap_num, depth))
    if for_swap:
        olsq.setSabreForSwap(True, swap_num)
    if for_mapping:
        olsq.useCircuitInitalMapping()
        setInitialMapping(circuit, initial_mapping)

