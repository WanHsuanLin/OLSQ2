timeout 1d python3 -u ../satmap/src/satmap.py ../quantum_cir_benchmark/arith_circuits/barenco_tof_4.qasm --arch ../OLSQ-dev/architecture/syc.txt --output_path ../OLSQ-dev/dac23/satmap_results/ >> ../OLSQ-dev/dac23/log/swap_opt/satmap_syc_btof_4.log
timeout 1d python3 -u ../satmap/src/satmap.py ../quantum_cir_benchmark/arith_circuits/barenco_tof_5.qasm --arch ../OLSQ-dev/architecture/syc.txt --output_path ../OLSQ-dev/dac23/satmap_results/ >> ../OLSQ-dev/dac23/log/swap_opt/satmap_syc_btof_5.log

timeout 1d python3 -u ../satmap/src/satmap.py ../QUEKO-benchmark/BNTF/54QBT_05CYC_QSE_0.qasm --arch ../OLSQ-dev/architecture/aspen.txt --output_path ../OLSQ-dev/dac23/satmap_results/ --timeout 86400 >> ../OLSQ-dev/dac23/log/swap_opt/satmap_k0_syc_queko_05.log
timeout 1d python3 -u ../satmap/src/satmap.py ../QUEKO-benchmark/BNTF/54QBT_15CYC_QSE_0.qasm --arch ../OLSQ-dev/architecture/aspen.txt --output_path ../OLSQ-dev/dac23/satmap_results/ --timeout 86400 >> ../OLSQ-dev/dac23/log/swap_opt/satmap_k0_syc_queko_15.log
timeout 1d python3 -u ../satmap/src/satmap.py ../QUEKO-benchmark/BNTF/54QBT_25CYC_QSE_0.qasm --arch ../OLSQ-dev/architecture/aspen.txt --output_path ../OLSQ-dev/dac23/satmap_results/ --timeout 86400 >> ../OLSQ-dev/dac23/log/swap_opt/satmap_k0_syc_queko_25.log
timeout 1d python3 -u ../satmap/src/satmap.py ../QUEKO-benchmark/BNTF/54QBT_35CYC_QSE_0.qasm --arch ../OLSQ-dev/architecture/aspen.txt --output_path ../OLSQ-dev/dac23/satmap_results/ --timeout 86400 >> ../OLSQ-dev/dac23/log/swap_opt/satmap_k0_syc_queko_35.log
timeout 1d python3 -u ../satmap/src/satmap.py ../QUEKO-benchmark/BNTF/54QBT_45CYC_QSE_0.qasm --arch ../OLSQ-dev/architecture/aspen.txt --output_path ../OLSQ-dev/dac23/satmap_results/ --timeout 86400 >> ../OLSQ-dev/dac23/log/swap_opt/satmap_k0_syc_queko_45.log

timeout 1d python3 -u ../satmap/src/satmap.py ../quantum_cir_benchmark/qaoa/qaoa_24_0.qasm --arch ../OLSQ-dev/architecture/syc.txt --output_path ../OLSQ-dev/dac23/satmap_results/ >> ../OLSQ-dev/dac23/log/swap_opt/satmap_syc_qaoa_24.log
timeout 1d python3 -u ../satmap/src/satmap.py ../quantum_cir_benchmark/qaoa/qaoa_28_0.qasm --arch ../OLSQ-dev/architecture/syc.txt --output_path ../OLSQ-dev/dac23/satmap_results/ >> ../OLSQ-dev/dac23/log/swap_opt/satmap_syc_qaoa_28.log


