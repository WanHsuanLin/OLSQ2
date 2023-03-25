

file_path="tmp"
log_path="int_euf_bv_log"


# for file in $1/*; do
for file in $file_path/*; do
    echo $file
    substring=".txt"
    log_file_tmp="${file%"$substring"}"
    log_file="${log_file_tmp#"$file_path/"}"
    
    # timeout 1d z3 -st -memory:62000 -smt2 $file >> $2/$log_file.log
    timeout 1d z3 -st -memory:62000 -v:10 -log -smt2 $file &> $log_path/$log_file.log
    #timeout 1d z3 tactic.default_tactic=sat sat.euf=true smt.ematching=false -st -memory:62000 -smt2 $file >> log/$log_file.log
done


