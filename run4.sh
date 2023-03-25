timeout 1d python3 -u run_olsq.py --dt grid --d 6 --f logs/scalability --swap --tran --sabre --b qaoa --use_sabre_mapping --qf  ../quantum_cir_benchmark/qaoa/qaoa_26_0.qasm >> logs/scalability/sabre_mapping_36_26.log
timeout 1d python3 -u run_olsq.py --dt grid --d 6 --f logs/scalability --swap --tran --sabre --b qaoa --use_sabre_mapping --qf  ../quantum_cir_benchmark/qaoa/qaoa_36_0.qasm >> logs/scalability/sabre_mapping_36_36.log
timeout 1d python3 -u run_olsq.py --dt grid --d 7 --f logs/scalability --swap --tran --sabre --b qaoa --use_sabre_mapping --qf  ../quantum_cir_benchmark/qaoa/qaoa_40_0.qasm >> logs/scalability/sabre_mapping_49_40.log




