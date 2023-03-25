
for file in $1/*; do
    echo $file
    substring=".txt"
    log_file_tmp="${file%"$substring"}"
    log_file="${log_file_tmp#"$file_path/"}"
    
    timeout 1d z3 -st -memory:62000 -smt2 $file >> $2/$log_file.log
    #timeout 1d z3 tactic.default_tactic=sat sat.euf=true smt.ematching=false -st -memory:62000 -smt2 $file >> log/$log_file.log
done
