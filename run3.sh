timeout 1d python3 -u run_olsq.py --dt grid --d 7 --f logs/scalability --swap --sabre --b qaoa --use_sabre_mapping --qf  ../quantum_cir_benchmark/qaoa/qaoa_40_0.qasm >> logs/scalability/sabre_mapping_exact_49_40.log
timeout 1d python3 -u run_olsq.py --dt grid --d 8 --f logs/scalability --swap --sabre --b qaoa --use_sabre_mapping --qf  ../quantum_cir_benchmark/qaoa/qaoa_50_0.qasm >> logs/scalability/sabre_mapping_exact_64_50.log
timeout 1d python3 -u run_olsq.py --dt grid --d 8 --f logs/scalability --swap --sabre --b qaoa --use_sabre_mapping --qf  ../quantum_cir_benchmark/qaoa/qaoa_60_0.qasm >> logs/scalability/sabre_mapping_exact_64_60.log
timeout 1d python3 -u run_olsq.py --dt grid --d 9 --f logs/scalability --swap --sabre --b qaoa --use_sabre_mapping --qf  ../quantum_cir_benchmark/qaoa/qaoa_70_0.qasm >> logs/scalability/sabre_mapping_exact_81_70.log
timeout 1d python3 -u run_olsq.py --dt grid --d 9 --f logs/scalability --swap --sabre --b qaoa --use_sabre_mapping --qf  ../quantum_cir_benchmark/qaoa/qaoa_80_0.qasm >> logs/scalability/sabre_mapping_exact_81_80.log




