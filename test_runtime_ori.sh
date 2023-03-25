
# python3 run_olsq_ori.py --dt grid --d 7 --f dac23/depth/tb_olsq_new --b qaoa --qf  ../quantum_cir_benchmark/qaoa/qaoa_16_0.qasm --dump 
# python3 run_olsq_ori.py --dt grid --d 7 --f dac23/depth/tb_olsq_new --b qaoa --qf  ../quantum_cir_benchmark/qaoa/qaoa_18_0.qasm --dump 
# python3 run_olsq_ori.py --dt grid --d 7 --f dac23/depth/tb_olsq_new --b qaoa --qf  ../quantum_cir_benchmark/qaoa/qaoa_20_0.qasm --dump 
# python3 run_olsq_ori.py --dt grid --d 7 --f dac23/depth/tb_olsq_new --b qaoa --qf  ../quantum_cir_benchmark/qaoa/qaoa_22_0.qasm --dump 
# python3 run_olsq_ori.py --dt grid --d 7 --f dac23/depth/tb_olsq_new --b qaoa --qf  ../quantum_cir_benchmark/qaoa/qaoa_24_0.qasm --dump 

python3 run_olsq_ori.py --dt grid --d 5 --f dac23/depth/tb_olsq_new --b qaoa --qf  ../quantum_cir_benchmark/qaoa/qaoa_16_0.qasm --dump
python3 run_olsq_ori.py --dt grid --d 5 --f dac23/depth/tb_olsq_new --b qaoa --qf  ../quantum_cir_benchmark/qaoa/qaoa_18_0.qasm --dump
python3 run_olsq_ori.py --dt grid --d 5 --f dac23/depth/tb_olsq_new --b qaoa --qf  ../quantum_cir_benchmark/qaoa/qaoa_20_0.qasm --dump
python3 run_olsq_ori.py --dt grid --d 5 --f dac23/depth/tb_olsq_new --b qaoa --qf  ../quantum_cir_benchmark/qaoa/qaoa_22_0.qasm --dump
python3 run_olsq_ori.py --dt grid --d 5 --f dac23/depth/tb_olsq_new --b qaoa --qf  ../quantum_cir_benchmark/qaoa/qaoa_24_0.qasm --dump

file_path="dac23/depth/tb_olsq_new"
log_path="dac23/log/depth/tb_olsq_new"

# for file in $1/*; do
for file in $file_path/*; do
    echo $file
    substring=".txt"
    log_file_tmp="${file%"$substring"}"
    log_file="${log_file_tmp#"$file_path/"}"
    
    # timeout 1d z3 -st -memory:62000 -smt2 $file >> $2/$log_file.log
    timeout 1d z3 -st -memory:62000 -smt2 $file >> $log_path/$log_file.log
    #timeout 1d z3 tactic.default_tactic=sat sat.euf=true smt.ematching=false -st -memory:62000 -smt2 $file >> log/$log_file.log
done


# timeout 1d z3 sat.cardinality.solver=false sat.cardinality.encoding=grouped -st -memory:62000 -smt2 cardinality_constraint_z3/25_20_21_30.txt >> cardinality_constraint_z3_log/25_20_21_30_grouped.log;
# timeout 1d z3 sat.cardinality.solver=false sat.cardinality.encoding=circuit -st -memory:62000 -smt2 cardinality_constraint_z3/25_20_21_30.txt >> cardinality_constraint_z3_log/25_20_21_30_circuit.log;
# timeout 1d z3 sat.cardinality.solver=false sat.cardinality.encoding=bimander -st -memory:62000 -smt2 cardinality_constraint_z3/25_20_21_30.txt >> cardinality_constraint_z3_log/25_20_21_30_bimander.log;
# timeout 1d z3 sat.cardinality.solver=false sat.cardinality.encoding=ordered -st -memory:62000 -smt2 cardinality_constraint_z3/25_20_21_30.txt >> cardinality_constraint_z3_log/25_20_21_30_ordered.log;
# timeout 1d z3 sat.cardinality.solver=false sat.cardinality.encoding=unate -st -memory:62000 -smt2 cardinality_constraint_z3/25_20_21_30.txt >> cardinality_constraint_z3_log/25_20_21_30_unate.log

# timeout 1d z3 sat.cardinality.solver=false sat.pb.solver=circuit -st -memory:62000 -smt2 cardinality_constraint_z3/25_20_21_30.txt >> cardinality_constraint_z3_log/25_20_21_30_pb_circuit.log;
# timeout 1d z3 sat.cardinality.solver=false sat.pb.solver=sorting -st -memory:62000 -smt2 cardinality_constraint_z3/25_20_21_30.txt >> cardinality_constraint_z3_log/25_20_21_30__pb_sorting.log;
# timeout 1d z3 sat.cardinality.solver=false sat.pb.solver=totalizer -st -memory:62000 -smt2 cardinality_constraint_z3/25_20_21_30.txt >> cardinality_constraint_z3_log/25_20_21_30_pb_totalizer.log;
# timeout 1d z3 sat.cardinality.solver=false sat.pb.solver=binary_merge -st -memory:62000 -smt2 cardinality_constraint_z3/25_20_21_30.txt >> cardinality_constraint_z3_log/25_20_21_30_pb_binary_merge.log;
# timeout 1d z3 sat.cardinality.solver=false sat.pb.solver=segmented -st -memory:62000 -smt2 cardinality_constraint_z3/25_20_21_30.txt >> cardinality_constraint_z3_log/25_20_21_30_pb_segmented.log

