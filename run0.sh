# timeout 1d python3 -u run_olsq.py --dt sycamore --f dac23/swap_opt --b arith --qf  ../quantum_cir_benchmark/arith_circuits/tof_5.qasm --sabre --tran --swap >> dac23/log/swap_opt/syc_arith_tof5_no_var_tran.log
# timeout 1d python3 -u run_olsq.py --dt aspen-4 --f dac23/swap_opt --b queko --qf  ../QUEKO-benchmark/BNTF/16QBT_45CYC_TFL_0.qasm --swap --sabre --tran >> dac23/log/swap_opt/aspen_queko_45_no_var_tran.log
# timeout 1d python3 -u run_olsq.py --dt sycamore --f dac23/swap_opt --b arith --qf  ../quantum_cir_benchmark/qft/qft_8.qasm --swap --sabre --tran >> dac23/log/swap_opt/syc_qft_no_var_tran.log

timeout 1d python3 -u run_olsq.py --dt grid --d 7 --f logs/scalability --swap --sabre --b qaoa --qf  ../quantum_cir_benchmark/qaoa/qaoa_40_0.qasm >> logs/scalability/exact_49_40.log
timeout 1d python3 -u run_olsq.py --dt grid --d 8 --f logs/scalability --swap --sabre --b qaoa --qf  ../quantum_cir_benchmark/qaoa/qaoa_50_0.qasm >> logs/scalability/exact_64_50.log
timeout 1d python3 -u run_olsq.py --dt grid --d 8 --f logs/scalability --swap --sabre --b qaoa --qf  ../quantum_cir_benchmark/qaoa/qaoa_60_0.qasm >> logs/scalability/exact_64_60.log
timeout 1d python3 -u run_olsq.py --dt grid --d 9 --f logs/scalability --swap --sabre --b qaoa --qf  ../quantum_cir_benchmark/qaoa/qaoa_70_0.qasm >> logs/scalability/exact_81_70.log
timeout 1d python3 -u run_olsq.py --dt grid --d 9 --f logs/scalability --swap --sabre --b qaoa --qf  ../quantum_cir_benchmark/qaoa/qaoa_80_0.qasm >> logs/scalability/exact_81_80.log




