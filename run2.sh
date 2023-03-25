timeout 1d python3 -u run_olsq.py --dt grid --d 11 --f logs/scalability --swap --tran --b qaoa --qf  ../quantum_cir_benchmark/qaoa/qaoa_110_0.qasm >> logs/scalability/121_110.log
timeout 1d python3 -u run_olsq.py --dt grid --d 11 --f logs/scalability --swap --tran --b qaoa --qf  ../quantum_cir_benchmark/qaoa/qaoa_120_0.qasm >> logs/scalability/121_120.log
