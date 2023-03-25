
timeout 1d python3 -u run_olsq.py --dt grid --d 12 --f logs/scalability --b qaoa --qf  ../quantum_cir_benchmark/qaoa/qaoa_140_0.qasm --swap --tran --sabre >> logs/scalability/144_130.log
timeout 1d python3 -u run_olsq.py --dt grid --d 12 --f logs/scalability --b qaoa --qf  ../quantum_cir_benchmark/qaoa/qaoa_140_0.qasm --swap --tran --sabre >> logs/scalability/144_140.log
timeout 1d python3 -u run_olsq.py --dt grid --d 13 --f logs/scalability --b qaoa --qf  ../quantum_cir_benchmark/qaoa/qaoa_160_0.qasm --swap --tran --sabre >> logs/scalability/169_150.log