[![dac](https://img.shields.io/badge/Published-DAC'23-brightgreen.svg?style=for-the-badge)]()

# OLSQ2: Scalable Optimal Layout Synthesis for NISQ Quantum Processors

Many quantum computers have constraints on the connections between qubits.
However, a quantum program may not conform to these constraints.
Thus, it is necessary to perform 'layout synthesis for quantum computing', LSQC, which transforms quantum programs prior to execution so that the connectivity issues are resolved.

OLSQ2 can solve LSQC optimally with respect to depth and number of SWAP gates.
There is also a transition-based mode (TB) to speed it up with little loss of optimality.

For more details on the theory and the experiments, please refer to [the paper]().
Below is a brief tutorial on how to use the package.

## Installation

Clone this repo:
```
git clone git@github.com:WanHsuanLin/OLSQ2.git
```
Please make sure that you have `pySAT` version `>=0.1.7` and `z3-solver` version `>=4.8.15.0` in your Python environment.
To reproduce the results reported in [the paper](), please install `z3-solver` with version `4.8.15.0`

## Initialization

```
from olsq import OLSQ

# initiate olsq with depth as objective, in normal mode
lsqc_solver = OLSQ("depth", "normal")
```

There are four argument in the constructor of OLSQ: `obj_is_swap`, `mode`, `encoding`, and `swap_up_bound`.
- `obj_is_swap`: `True` to set SWAP count as objective or `False` to set depth as objective. When optimizing SWAP count, OLSQ2 will save each intermediate compilation result to a qasm file with the name "intermediate_result_swap_count_{SWAP_COUNT}.qasm"
- `mode`:  `"normal"` or `"transition"`. The latter stands for TB-OLSQ in the paper, which is usually much faster with little loss of optimality.
- `encoding`: Different strategies for [pySAT](https://pysathq.github.io/docs/html/api/card.html#pysat.card.CardEnc) to encode cardinality constraint by CNF. Options: `1`, `2`, `3`, `6`, `7`, `8` and `9`.
- `swap_up_bound`:  Users can specify the starting point for SWAP optimization.

## Setting the device

To perform LSQC, we need to know the connections between the qubits, which is information about the physical device.
We are going to use the `setdevice` method by directly construct a device with some properties.

```
from olsq.device import qcdevice

# directly construct a device from properties needed by olsq
lsqc_solver.setdevice( qcdevice(name="dev", nqubits=5, 
     connection=[(0, 1), (1, 2), (1, 3), (3, 4)], swap_duration=3) )
```

## Setting the Input Program

Apart from the device, we need the quantum program/circuit to execute, which can be set with the `setprogram` method.

OLSQ has an intermediate representation (IR) of quantum programs. (For details, refer to [a later part](#olsq-ir) of this tutorial.)
In general, there are four ways to set the program: 
1. Use OLSQ IR
2. Use a string in QASM format

```
circuit_str = "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[3];\nh q[2];\n" \
              "cx q[1], q[2];\ntdg q[2];\ncx q[0], q[2];\nt q[2];\n" \
              "cx q[1], q[2];\ntdg q[2];\ncx q[0], q[2];\nt q[1];\nt q[2];\n" \
              "cx q[0], q[1];\nh q[2];\nt q[0];\ntdg q[1];\ncx q[0], q[1];\n"

# input the quantum program as a QASM string
lsqc_solver.setprogram(circuit_str)
```

The example above is a Toffoli gate.
We can also load an QASM file of it.
```
# load one of the QASM files from olsq/benchmarks
lsqc_solver.setprogram("toffoli", input_mode="benchmark")

# load your own QASM file
# circuit_file = open("my-qasm-file", "r").read()

lsqc_solver.setprogram(circuit_file)

# Toffoli Gate:
#                                                        ┌───┐      
# q_0: ───────────────────■─────────────────────■────■───┤ T ├───■──
#                         │             ┌───┐   │  ┌─┴─┐┌┴───┴┐┌─┴─┐
# q_1: ───────■───────────┼─────────■───┤ T ├───┼──┤ X ├┤ TDG ├┤ X ├
#      ┌───┐┌─┴─┐┌─────┐┌─┴─┐┌───┐┌─┴─┐┌┴───┴┐┌─┴─┐├───┤└┬───┬┘└───┘
# q_2: ┤ H ├┤ X ├┤ TDG ├┤ X ├┤ T ├┤ X ├┤ TDG ├┤ X ├┤ T ├─┤ H ├──────
#      └───┘└───┘└─────┘└───┘└───┘└───┘└─────┘└───┘└───┘ └───┘      
"""
```

## Solving and Output

It can be seen that in the Toffoli gate above, there are two-qubit gates on pair `(q_0,q_1)`, `(q_1,q_2)`, and `(q_2,q_0)`.
However, there are no such triangles on device `ourense`.
This means that no matter how the qubits in the program are mapped to physical qubits, we need to insert SWAP gates.

```
# solve LSQC
result = lsqc_solver.solve()
```

The `solve` method can take two optional arguemnts
- `use_sabre`: `True` to use SABRE to get the upper bound of the SWAP count for SWAP optimization.
- `output_mode`: can be `"IR"`. Refer [here](#olsq-ir) on what would be returned in this case.
- `output_file_name`
If `output_mode` is default, the return is a tuple of three things:
- A string representing the output quantum program in QASM format.
If `output_file_name` is provided, then the QASM string would be written to that file.
- final_mapping: from each program qubit to the corresponding physical qubit at the end of execution.
- objective_value

The result of the Toffoli example is shown below.
Note that a SWAP gate, decomposed into three CX gates, has been inserted.
```
# a LSQC solution to the Toffoli gate on device 'ourense'
#                                                  ┌───┐     ┌───┐┌───┐ ┌───┐      ┌─┐      
# q_0: ───────────────────■─────────────────────■──┤ X ├──■──┤ X ├┤ T ├─┤ H ├──────┤M├──────
#      ┌───┐┌───┐┌─────┐┌─┴─┐┌───┐┌───┐┌─────┐┌─┴─┐└─┬─┘┌─┴─┐└─┬─┘└───┘ ├───┤      └╥┘┌─┐   
# q_1: ┤ H ├┤ X ├┤ TDG ├┤ X ├┤ T ├┤ X ├┤ TDG ├┤ X ├──■──┤ X ├──■────■───┤ T ├───■───╫─┤M├───
#      └───┘└─┬─┘└─────┘└───┘└───┘└─┬─┘└┬───┬┘└───┘     └───┘     ┌─┴─┐┌┴───┴┐┌─┴─┐ ║ └╥┘┌─┐
# q_2: ───────■─────────────────────■───┤ T ├─────────────────────┤ X ├┤ TDG ├┤ X ├─╫──╫─┤M├
#                                       └───┘                     └───┘└─────┘└───┘ ║  ║ └╥┘
# q_3: ─────────────────────────────────────────────────────────────────────────────╫──╫──╫─
#                                                                                   ║  ║  ║
# q_4: ─────────────────────────────────────────────────────────────────────────────╫──╫──╫─
#                                                                                   ║  ║  ║
# c: 5/═════════════════════════════════════════════════════════════════════════════╩══╩══╩═
#                                                                                   2  0  1
```


## TB-OLSQ2

The transition-based mode is enabled if chosen at the initiation of `OLSQ`.
Roughly speaking, we only use a kind of coarse-grain time in this mode, so the runtime is much shorter.
The returned QASM string and `final_mapping` should be similar to what they were before.
Only if the objective is `"depth"`, the objective value would be very different from the normal mode.
There is only one SWAP inserted, so there are only two coarse-grain time steps, separated by the SWAP, whereas there are 14 time steps if using exact time.

## OLSQ IR

OLSQ IR contains three things:
1. `count_program_qubit`: the number of qubits in the program.
2. `gates`: a list of tuples representing qubit(s) acted on by a gate, each tuple has one index if it is a single-qubit gate, two indices if it is a two-qubit gate.
3. `gate_spec`: list of type/name of each gate, which is not important to OLSQ, and only needed when generating output.

```
# For the following circuit
# q_0: ───────────────────■───
#                         │  
# q_1: ───────■───────────┼───
#      ┌───┐┌─┴─┐┌─────┐┌─┴─┐
# q_2: ┤ H ├┤ X ├┤ TDG ├┤ X ├─
#      └───┘└───┘└─────┘└───┘ 

# count_program_qubit = 3
# gates = ((2,), (1,2), (2,), (0,1))
# gate_spec = ("h", "cx", "tdg", "cx")
```

If in the `solve` method, `output_mode` is set to `"IR"`, the return is a tuple of five things
1. `result_depth`: depth of the resulting quantum program
2. `list_scheduled_gate_name`: similar to `gate_spec` in the IR
3. `list_scheduled_gate_qubits`: similar to `gates` in the IR
4. `final_mapping`
5. `objective_value`

## Example: run_olsq.py

run_olsq.py is an example program to use OLSQ2/TB-OLSQ2 to perform layout synthesis.
```
# compile an qaoa circuit on a 5-by-5 grid quantum device by TB-OLSQ2 using swap as objective and SABRE's result for the starting point of optimization. The output file is IR and will be store in example/.
python3 run_olsq.py --dt grid --d 4 --f example/ --qf benchmark/qaoa/qaoa_16_0.qasm --swap --sabre --tran
# The output files (Final IR output file and the intermediate qasm file) of running the command are in example/.

# compile an qaoa circuit on sycamore quantum device by TB-OLSQ2 using swap as objective and store the output IR file in the current directory
python3 run_olsq.py --dt sycamore --f . --qf benchmark/qaoa/qaoa_16_0.qasm --tran
```
- `--tran`: Use TB-OLSQ2.
- `--swap`: Set SWAP count as objective.
- `--dt $(str)`: Type of the quantum device: ourense, sycamore, rochester, tokyo, aspen-4, eagle, or grid. When using a grid architecure, add `--d $(int)` to specify the grid length.
- `--d $(int)`: Grid length of the grid architecture
- `--qasm $(str)`: Input QASM file name
- `--f $(str)`: The location to stroe the output IR file. Default: current directory
- `--swap_duration $(int)`: SWAP duration. Default: 1 
- `--sabre`: Use sabre to get SWAP upper bound
- `--encoding $(int)`: Different encoding strategies to convert cardinality constraint to CNF. seqcounter = 1, sortnetwrk  = 2, cardnetwrk = 3, totalizer = 6, mtotalizer = 7. kmtotalizer = 8, native = 9
- `--swap_bound $(int)`: Specify user-defined SWAP count upper bound for SWAP optimization


## BibTeX Citation
```
```
