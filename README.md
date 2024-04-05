[![iccad](https://img.shields.io/badge/Published-ICCAD'20-brightgreen.svg?style=for-the-badge)](https://ieeexplore.ieee.org/document/9256696)
[![arXiv](https://img.shields.io/badge/arXiv-2007.15671-brightgreen.svg?style=for-the-badge)](https://arxiv.org/abs/2007.15671)
[![Unitary Fund](https://img.shields.io/badge/Supported%20By-UNITARY%20FUND-brightgreen.svg?style=for-the-badge)](http://unitary.fund)

# OLSQ2: Scalable Optimal Layout Synthesis for NISQ Quantum Processors

Many quantum computers have constraints on the connections between qubits.
However, a quantum program may not conform to these constraints.
Thus, it is necessary to perform 'quantum layout synthesis', QLS, which transforms quantum programs prior to execution so that the connectivity issues are resolved.

OLSQ2 can solve QLS optimally with respect to depth and number of SWAP gates.
There is also a transition-based mode (TB) to speed it up with little loss of optimality.

For more details on the theory and the experiments, please refer to [the paper]().
Below is a brief tutorial on how to use the package.

## Setup

```
git clone --recursive https://github.com/WanHsuanLin/OLSQ2.git
git checkout cpp
mkdir build && cd build
cmake ..
```

- If you have an existing pybind11 installation, you do not need `--recursive`.
- You must have an existing z3 installation, which you can specify using `CMAKE_PREFIX_PATH` or through other cmake mechanisms.
- Please make sure that you have `networkx` version `>=2.5` and `z3-solver` version `>=4.8.9.0` in your Python environment.

## Example: run_olsq.py

run_olsq.py is an example program to use OLSQ2/TB-OLSQ2 to perform layout synthesis.
```
# compile an qaoa circuit on a 5-by-5 grid quantum device by TB-OLSQ2 using swap as objective and SABRE's result for the starting point of optimization. The output file is IR and will be stored in example/.
python3 run_olsq.py --dt grid --d 5 --f example/ --qf benchmark/qaoa/qaoa_16_0.qasm --swap --sabre --tran
# The output files (Final IR output file and the intermediate qasm file) of running the command are in example/.

# compile an qaoa circuit on sycamore quantum device by TB-OLSQ2 using swap as objective and store the output IR file in the current directory
python3 run_olsq.py --dt sycamore --f . --qf benchmark/qaoa/qaoa_16_0.qasm --tran
```
- `--tran`: Use TB-OLSQ2.
- `--swap`: Set SWAP count as objective.
- `--dt $(str)`: Type of the quantum device: ourense, sycamore, rochester, tokyo, aspen-4, eagle, or grid. When using a grid architecture, add `--d $(int)` to specify the grid length.
- `--d $(int)`: Grid length of the grid architecture
- `--qasm $(str)`: Input QASM file name
- `--f $(str)`: The location to store the output IR file. Default: current directory
- `--swap_duration $(int)`: SWAP duration. Default: 1 
- `--sabre`: Use sabre to get SWAP upper bound
- `--encoding $(int)`: Different encoding strategies to convert cardinality constraint to CNF. seqcounter = 1, sortnetwrk  = 2, cardnetwrk = 3, totalizer = 6, mtotalizer = 7. kmtotalizer = 8, native = 9
- `--swap_bound $(int)`: Specify user-defined SWAP count upper bound for SWAP optimization


## BibTeX Citation
```
@INPROCEEDINGS{10247760,
  author={Lin, Wan-Hsuan and Kimko, Jason and Tan, Bochen and Bjørner, Nikolaj and Cong, Jason},
  booktitle={2023 60th ACM/IEEE Design Automation Conference (DAC)}, 
  title={Scalable Optimal Layout Synthesis for NISQ Quantum Processors}, 
  year={2023},
  volume={},
  number={},
  pages={1-6},
  doi={10.1109/DAC56929.2023.10247760}}
```
