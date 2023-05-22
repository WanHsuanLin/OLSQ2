/***********************************************************************
  File        [ olsq.cpp ]
  System      [ OLSQ: optimal quantum layout synthesis tool]
  Package     [ olsq ]
  Synopsis    [ OLSQ class implementation ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 22, Nov., 2022 ]
***********************************************************************/
#include "olsq/olsq.hpp"

OLSQ_NAMESPACE_CPP_START

void OLSQ::setDependency(vector<pair<unsigned_t, unsigned_t> >& vDependencies){
    for (pair<unsigned_t, unsigned_t>& p : vDependencies){
        assert(isValidGateIdx(p.first));
        assert(isValidGateIdx(p.second));
        addDependency(p.first, p.second);
    }
}

void OLSQ::run(string const & fileName){
    _timer.start(TimeUsage::FULL);
    fprintf(stdout, "[Info] OLSQ Layout Synthesis                        \n");
    if(!_olsqParam.is_given_dependency){
        // _timer.start(TimeUsage::PARTIAL);
        constructDependency();
        fprintf(stdout, "[Info] Constructing dependency                        \n");
        // _timer.showUsage("Constructing dependency", TimeUsage::PARTIAL);
        if(_verbose == 2){
            printDependency();
        }
    }
    runSMT();
    if(_verbose > 0){
        _pCircuit->printCircuitLayout();
    }
    if (_olsqParam.is_transition){
        asapScheduling();
    }
    // cout << "finish scheduling" << endl;
    _pCircuit->printCircuitLayout();
    outputQASM(fileName);
}

void OLSQ::runSMT(){
    fprintf(stdout, "[Info] OLSQ Layout Synthesis                        \n");
    if(!_olsqParam.is_transition){
        if (!_olsqParam.is_given_depth){
            _olsqParam.min_depth = extract_longest_chain();
            fprintf(stdout, "[Info] Longest chain = %d\n", _olsqParam.min_depth);
        }
        if (_olsqParam.max_depth == 0){
            _olsqParam.max_depth = 2 * _olsqParam.min_depth;
        }
    }
    bool solve = false;
    unsigned_t iter = 0;
    while (!solve){
        fprintf(stdout, "[Info] Iter %d: Solving with depth range (%d, %d)            \n", iter, _olsqParam.min_depth, _olsqParam.max_depth);
        _timer.start(TimeUsage::PARTIAL);
        fprintf(stdout, "[Info] Iter %d: Generating formulation                        \n", iter);
        generateFormulationZ3();
        _timer.showUsage("Generating formulation", TimeUsage::PARTIAL);
        // _timer.showUsage("Generating formulation", TimeUsage::FULL);
        _timer.start(TimeUsage::PARTIAL);
        fprintf(stdout, "[Info] Iter %d: Optimizing model                             \n", iter);
        solve = optimize();
        if(!solve)
            increase_depth_bound();
        ++iter;
    }
}

void OLSQ::generateFormulationZ3(){
    fprintf(stdout, "[Info]          constructing variables                       \n");
    constructVariableZ3();
    fprintf(stdout, "[Info]          constructing injective mapping constraint    \n");
    addInjectiveMappingConstraintsZ3();
    fprintf(stdout, "[Info]          constructing valid two-qubit gate constraint \n");
    addValidTwoQubitGateConstraintsZ3();
    fprintf(stdout, "[Info]          constructing dependency constraint           \n");
    addDependencyConstraintsZ3();
    // string s = Z3_solver_to_string(_smt.c, _smt.smtSolver);
    // cout << s;
    // getchar();
    fprintf(stdout, "[Info]          constructing swap overlapping constraint     \n");
    addSwapConstraintsZ3();
    // s = Z3_solver_to_string(_smt.c, _smt.smtSolver);
    // cout << s;
    // getchar();
    fprintf(stdout, "[Info]          constructing mapping transformation constraint\n");
    addTransformationConstraintsZ3();
    // s = Z3_solver_to_string(_smt.c, _smt.smtSolver);
    // cout << s;
    // getchar();
}

void OLSQ::constructVariableZ3(){
    unsigned_t bit_length_pi, bit_length_time;
    bit_length_pi = ceil(log2(_device.nQubit() + 1));
    bit_length_time = ceil(log2(_olsqParam.max_depth));
    // fprintf(stdout, "[Info]          bit_length_pi: %d, bit_length_time: %d\n", bit_length_pi, bit_length_time);
    // _smt.vvPi.resize(_olsqParam.max_depth, vector<expr> (_pCircuit->nProgramQubit()));
    _smt.vvPi.reserve(_olsqParam.max_depth);
    _smt.vTg.reserve(_pCircuit->nGate());
    // _smt.vvSigma.resize(_olsqParam.max_depth, vector<expr> (_device.nEdge()));
    _smt.vvSigma.reserve(_olsqParam.max_depth);

    string s;
    unsigned_t i, j;

    // Create a bit-vector sort of size 1.
    BitwuzlaSort *sortbvpi = bitwuzla_mk_bv_sort(_smt.pSolver, bit_length_pi);
    BitwuzlaSort *sortbvtime = bitwuzla_mk_bv_sort(_smt.pSolver, bit_length_time);
    BitwuzlaSort *sortbool = bitwuzla_mk_bool_sort(_smt.pSolver);

    for (i = 0; i < _olsqParam.max_depth; ++i){
        _smt.vvPi.emplace_back(vector<expr> ());
        for (j = 0; j < _pCircuit->nProgramQubit(); ++j){
            s = _pCircuit->name() + "_map_t" + to_string(i) + "_q" + to_string(j);
            // _smt.vvPi[i].emplace_back(_smt.c.bv_const(s.c_str(), bit_length_pi));
            _smt.vvPi[i].emplace_back(bitwuzla_mk_const(_smt.pSolver, sortbvpi, s));
        }
    }

    for (i = 0; i < _pCircuit->nGate(); ++i){
        s = _pCircuit->name() + "_time_" + to_string(i);
        // _smt.vTg.emplace_back(_smt.c.bv_const(s.c_str(), bit_length_time));
        _smt.vTg.emplace_back(bitwuzla_mk_const(_smt.pSolver, sortbvtime, s));
    }

    for (i = 0; i < _olsqParam.max_depth; ++i){
        _smt.vvSigma.emplace_back(vector<expr> ());
        for (j = 0; j < _device.nEdge(); ++j){
            s = _pCircuit->name() + "_ifswap_t" + to_string(i) + "_e" + to_string(j);
            // _smt.vvSigma[i].emplace_back(_smt.c.bool_const(s.c_str()));
            _smt.vvSigma[i].emplace_back(bitwuzla_mk_const(_smt.pSolver, sortbvtime, s));
        }
    }

}

void OLSQ::addInjectiveMappingConstraintsZ3(){
    unsigned_t i = 0, j, k, qId;
    int_t pId;
    BitwuzlaSort *sortbvpi = bitwuzla_mk_bv_sort(_smt.pSolver, ceil(log2(_device.nQubit() + 1)));
    if(_olsqParam.is_given_initial_mapping){
        for(qId = 0;  qId < _pCircuit->nProgramQubit(); ++qId){
            pId = _pCircuit->initialMapping(qId);
            assert(pId > -1);
            // _smt.smtSolver.add(_smt.vvPi[0][qId] == pId);  
            BitwuzlaTerm * pIdBv = bitwuzla_mk_bv_value_uint64(_smt.pSolver, sortbvpi, pId);
            bitwuzla_assert(_smt.pSolver, bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[0][qId], pId));
        }
        ++i;
    }
    BitwuzlaTerm * zero = bitwuzla_mk_bv_zero(_smt.pSolver, sortbvpi);
    BitwuzlaTerm * nqubit = bitwuzla_mk_bv_value_uint64(_smt.pSolver, sortbvpi, _device.nQubit());
    for (; i < _olsqParam.max_depth; ++i){
        for (j = 0; j < _pCircuit->nProgramQubit(); ++j){
            // _smt.smtSolver.add(ule(0, _smt.vvPi[i][j]));
            // _smt.smtSolver.add(ult(_smt.vvPi[i][j], _device.nQubit()));
            bitwuzla_assert(_smt.pSolver, bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_BV_ULE, zero, _smt.vvPi[i][j]));
            bitwuzla_assert(_smt.pSolver, bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_BV_ULT, _smt.vvPi[i][j], nqubit));
            for (k = 0; k < j; ++k){
                // _smt.smtSolver.add(_smt.vvPi[i][j] != _smt.vvPi[i][k]);
                bitwuzla_assert(_smt.pSolver, bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_DISTINCT, _smt.vvPi[i][j], _smt.vvPi[i][k]));
            }
        }
        // expr_vector vDistinctPi(_smt.c);
        // unsigned_t t, e;
        // for (j = 0; j < _pCircuit->nProgramQubit(); ++j){
        //     vDistinctPi.push_back(_smt.vvPi[i][j]);
        // }
        // _smt.smtSolver.add(distinct(vDistinctPi));
    }
}

void OLSQ::addValidTwoQubitGateConstraintsZ3(){
    Gate gate;
    Edge edge;
    unsigned_t i, t, j;
    BitwuzlaSort *sortbvtime = bitwuzla_mk_bv_sort(_smt.pSolver, ceil(log2(_olsqParam.max_depth)));
    BitwuzlaSort *sortbvpi = bitwuzla_mk_bv_sort(_smt.pSolver, ceil(log2(_device.nQubit() + 1)));
    for ( i = 0; i < _pCircuit->nGate();  ++i ){
        Gate & gate = _pCircuit->gate(i);
        if ((gate).nTargetQubit() == 2){
            for (t = 0; t < _olsqParam.max_depth; ++t){
                // expr clause = !(_smt.vTg[i] == (int_t)t);  
                BitwuzlaTerm * bvt = bitwuzla_mk_bv_value_uint64(_smt.pSolver, sortbvtime, t);
                BitwuzlaTerm * clause = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vTg[i], bvt);
                clause = bitwuzla_mk_term1(_smt.pSolver, BITWUZLA_KIND_NOT, clause);
                BitwuzlaTerm *pro1EqPhy1, *pro2EqPhy2, *pro1EqPhy2, *pro2EqPhy1, *cond1, *cond2;
                for ( j = 0; j < _device.nEdge();  ++j ){
                    Edge & edge = _device.edge(j);
                    BitwuzlaTerm * q1Bv = bitwuzla_mk_bv_value_uint64(_smt.pSolver, sortbvpi, edge.qubitId1());
                    BitwuzlaTerm * q2Bv = bitwuzla_mk_bv_value_uint64(_smt.pSolver, sortbvpi, edge.qubitId2());
                    pro1EqPhy1 = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[t][gate.targetProgramQubit(0)], q1Bv);
                    pro2EqPhy2 = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[t][gate.targetProgramQubit(1)], q2Bv);
                    pro1EqPhy2 = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[t][gate.targetProgramQubit(0)], q2Bv);
                    pro2EqPhy1 = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[t][gate.targetProgramQubit(1)], q1Bv);
                    cond1 = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_AND, pro1EqPhy1, pro2EqPhy2);
                    cond2 = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_AND, pro2EqPhy1, pro1EqPhy2);
                    clause = bitwuzla_mk_term3(_smt.pSolver, BITWUZLA_KIND_OR, clause, cond1, cond2);
                    // clause = clause || ((_smt.vvPi[t][gate.targetProgramQubit(0)] == (int_t)edge.qubitId1()) && 
                    //                 (_smt.vvPi[t][gate.targetProgramQubit(1)] == (int_t)edge.qubitId2()) )
                    //                 || ((_smt.vvPi[t][gate.targetProgramQubit(1)] == (int_t)edge.qubitId1()) && 
                    //                 (_smt.vvPi[t][gate.targetProgramQubit(0)] == (int_t)edge.qubitId2()) );
                }                       
                // _smt.smtSolver.add(clause);
                bitwuzla_assert(_smt.pSolver, clause);
            }
        }
    }    
}

void OLSQ::addDependencyConstraintsZ3(){
    unsigned_t i;
    Gate g;
    BitwuzlaSort * sortbvtime = bitwuzla_mk_bv_sort(_smt.pSolver, ceil(log2(_olsqParam.max_depth)));
    BitwuzlaTerm * zero = bitwuzla_mk_bv_zero(_smt.pSolver, sortbvtime);
    for (i = 0; i < _pCircuit->nGate(); ++i){
        // _smt.smtSolver.add(ule(0, _smt.vTg[i]));
        bitwuzla_assert(_smt.pSolver, bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_BV_ULE, zero, _smt.vTg[i]));
    }
    // printDependency();
    // cerr << "line 212" << endl;
    // cerr << _vpGateDependency.size() << endl;
    if (_olsqParam.is_transition){
        for ( i = 0; i < _vpGateDependency.size();  ++i ){
            // cerr << i << ": " << _vpGateDependency[i].first << ", " <<_vpGateDependency[i].second << endl;
            // cerr << i << ": " << _smt.vTg[_vpGateDependency[i].first] << ", " << _smt.vTg[_vpGateDependency[i].second] << endl;
            // _smt.smtSolver.add(ule(_smt.vTg[_vpGateDependency[i].first], _smt.vTg[_vpGateDependency[i].second]));
            bitwuzla_assert(_smt.pSolver, bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_BV_ULE, _smt.vTg[_vpGateDependency[i].first], _smt.vTg[_vpGateDependency[i].second]));
        }
    }
    else{
        for ( i = 0; i < _vpGateDependency.size();  ++i ){
            // _smt.smtSolver.add(ult(_smt.vTg[_vpGateDependency[i].first], _smt.vTg[_vpGateDependency[i].second]));
            bitwuzla_assert(_smt.pSolver, bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_BV_ULT, _smt.vTg[_vpGateDependency[i].first], _smt.vTg[_vpGateDependency[i].second]));
        }
    }
}

void OLSQ::addSwapConstraintsZ3(){
    // No swap for t<s
    unsigned_t i, j, t, e, tt, q1, q2;
    Qubit qubit;
    unsigned_t bound = (_olsqParam.swap_duration < _olsqParam.max_depth) ? _olsqParam.swap_duration - 1 : _olsqParam.max_depth;
    for (i = 0; i < bound; ++i){
        for (j = 0; j < _device.nEdge(); ++j){
            // _smt.smtSolver.add(!_smt.vvSigma[i][j]);
            bitwuzla_assert(_smt.pSolver, bitwuzla_mk_term1(_smt.pSolver, BITWUZLA_KIND_NOT, _smt.vvSigma[i][j]));
        }
    }
    // swap gates can not overlap with swap in space
    for (t = _olsqParam.swap_duration -1; t < _olsqParam.max_depth; ++t){
        for (e = 0; e < _device.nEdge(); ++e){
            q1 = _device.edge(e).qubitId1();
            q2 = _device.edge(e).qubitId2();
            for (unsigned_t ee : _device.qubit(q1).vSpanEdge){
                for (tt = t - _olsqParam.swap_duration + 1; tt < t + 1; ++tt){
                    if (ee < e){
                        // cout << "add nonoverlap constraint for " << ee << " and " << e << endl;
                        // _smt.smtSolver.add(!_smt.vvSigma[t][e] || !_smt.vvSigma[tt][ee]);  
                        bitwuzla_assert(_smt.pSolver, 
                            bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_OR, 
                                bitwuzla_mk_term1(_smt.pSolver, BITWUZLA_KIND_NOT, _smt.vvSigma[t][e]),
                                bitwuzla_mk_term1(_smt.pSolver, BITWUZLA_KIND_NOT, _smt.vvSigma[tt][ee])));
                        // sOverlapEdge.insert((make_pair(ee,e)))  
                    }
                }
            }
            for (unsigned_t ee : _device.qubit(q2).vSpanEdge){
                for (tt = t - _olsqParam.swap_duration + 1; tt < t + 1; ++tt){
                    // if (ee < e && sOverlapEdge.count((make_pair(ee,e))) == 0)
                    if (ee < e){
                        // cout << "add nonoverlap constraint for " << ee << " and " << e << endl;
                        // _smt.smtSolver.add(!_smt.vvSigma[t][e] || !_smt.vvSigma[tt][ee]);    
                        bitwuzla_assert(_smt.pSolver, 
                            bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_OR, 
                                bitwuzla_mk_term1(_smt.pSolver, BITWUZLA_KIND_NOT, _smt.vvSigma[t][e]),
                                bitwuzla_mk_term1(_smt.pSolver, BITWUZLA_KIND_NOT, _smt.vvSigma[tt][ee])));
                        // sOverlapEdge.insert((make_pair(ee,e)))  
                    }
                }
            }
        }
    }
    if (!_olsqParam.is_transition){
        BitwuzlaSort *sortbvtime = bitwuzla_mk_bv_sort(_smt.pSolver, ceil(log2(_olsqParam.max_depth)));
        BitwuzlaSort *sortbvpi = bitwuzla_mk_bv_sort(_smt.pSolver, ceil(log2(_device.nQubit() + 1)));
        // swap gates can not overlap with swap in time
        for (t = _olsqParam.swap_duration -1; t < _olsqParam.max_depth; ++t){
            for (e = 0; e < _device.nEdge(); ++e){
                for (tt = t - _olsqParam.swap_duration + 1; tt < t; ++tt){
                    // _smt.smtSolver.add(!_smt.vvSigma[t][e] || !_smt.vvSigma[tt][e]);
                    bitwuzla_assert(_smt.pSolver, 
                            bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_OR, 
                                bitwuzla_mk_term1(_smt.pSolver, BITWUZLA_KIND_NOT, _smt.vvSigma[t][e]),
                                bitwuzla_mk_term1(_smt.pSolver, BITWUZLA_KIND_NOT, _smt.vvSigma[tt][e])));
                }
            }
        }
        // swap gates can not ovelap with other gates
        // the time interval should be modified
        for (t = _olsqParam.swap_duration -1; t < _olsqParam.max_depth; ++t){
            for (e = 0; e < _device.nEdge(); ++e){
               for (i = 0; i < _pCircuit->nGate(); ++i){
                    Gate& gate = _pCircuit->gate(i);
                    // cout << _pCircuit->nGate() << endl;
                    for (tt = t - _olsqParam.swap_duration + 1; tt < t + 1; ++tt){
                        // cout << gate.name() << " " << gate.nTargetQubit() << endl;
                        BitwuzlaTerm * bvt = bitwuzla_mk_bv_value_uint64(_smt.pSolver, sortbvtime, tt);
                        BitwuzlaTerm * q1Bv = bitwuzla_mk_bv_value_uint64(_smt.pSolver, sortbvpi, _device.edge(e).qubitId1());
                        BitwuzlaTerm * q2Bv = bitwuzla_mk_bv_value_uint64(_smt.pSolver, sortbvpi, _device.edge(e).qubitId2());
                        BitwuzlaTerm * clause = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vTg[i], bvt);
                        clause = bitwuzla_mk_term1(_smt.pSolver, BITWUZLA_KIND_NOT, clause);
                        BitwuzlaTerm *pro1EqPhy1, *pro2EqPhy2, *pro1EqPhy2, *pro2EqPhy1, *cond;
                    // for ( j = 0; j < _device.nEdge();  ++j ){
                    //     Edge & edge = _device.edge(j);
                    //     BitwuzlaTerm * q1Bv = bitwuzla_mk_bv_value_uint64(_smt.pSolver, sortbvpi, edge.qubitId1());
                    //     BitwuzlaTerm * q2Bv = bitwuzla_mk_bv_value_uint64(_smt.pSolver, sortbvpi, edge.qubitId2());
                    //     pro1EqPhy1 = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[t][gate.targetProgramQubit(0)], q1Bv);
                    //     pro2EqPhy2 = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[t][gate.targetProgramQubit(1)], q2Bv);
                    //     pro1EqPhy2 = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[t][gate.targetProgramQubit(0)], q2Bv);
                    //     pro2EqPhy1 = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[t][gate.targetProgramQubit(1)], q1Bv);
                    //     cond1 = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_AND, pro1EqPhy1, pro2EqPhy2);
                    //     cond2 = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_AND, pro2EqPhy1, pro1EqPhy2);
                        pro1EqPhy1 = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[tt][gate.targetProgramQubit(0)], q1Bv);
                        pro1EqPhy2 = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[tt][gate.targetProgramQubit(0)], q2Bv);
                        if (gate.nTargetQubit() == 1){
                            cond = bitwuzla_mk_term3(_smt.pSolver, BITWUZLA_KIND_OR, pro1EqPhy1, pro1EqPhy2, 
                                                bitwuzla_mk_term1(_smt.pSolver, BITWUZLA_KIND_NOT, _smt.vvSigma[t][e]));
                            bitwuzla_assert(_smt.pSolver, 
                                bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_AND, clause, cond));
                            // _smt.smtSolver.add( !((_smt.vTg[i] == (int_t)tt) && 
                            //                     ((_smt.vvPi[tt][gate.targetProgramQubit(0)] == (int_t)_device.edge(e).qubitId1())  ||
                            //                     (_smt.vvPi[tt][gate.targetProgramQubit(0)] == (int_t)_device.edge(e).qubitId2()) )) || !_smt.vvSigma[t][e]);
                        }
                        else if(gate.nTargetQubit() == 2){
                            pro2EqPhy2 = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[tt][gate.targetProgramQubit(1)], q2Bv);
                            pro2EqPhy1 = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[tt][gate.targetProgramQubit(1)], q1Bv);
                            cond = bitwuzla_mk_term3(_smt.pSolver, BITWUZLA_KIND_OR, pro1EqPhy1, pro1EqPhy2, 
                                                bitwuzla_mk_term1(_smt.pSolver, BITWUZLA_KIND_NOT, _smt.vvSigma[t][e]));
                            cond = bitwuzla_mk_term3(_smt.pSolver, BITWUZLA_KIND_OR, pro2EqPhy1, pro2EqPhy2, cond);
                            bitwuzla_assert(_smt.pSolver, 
                                bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_AND, clause, cond));
                            // _smt.smtSolver.add( !((_smt.vTg[i] == (int_t)tt) && 
                            //                     ((_smt.vvPi[tt][gate.targetProgramQubit(0)] == (int_t)_device.edge(e).qubitId1())  ||
                            //                     (_smt.vvPi[tt][gate.targetProgramQubit(1)] == (int_t)_device.edge(e).qubitId2())  ||
                            //                     (_smt.vvPi[tt][gate.targetProgramQubit(1)] == (int_t)_device.edge(e).qubitId1())  ||
                            //                     (_smt.vvPi[tt][gate.targetProgramQubit(0)] == (int_t)_device.edge(e).qubitId2()) )) || !_smt.vvSigma[t][e]);
                        }
                        else{
                            assert(false);
                        }
                    } 
                } 
            }
        }
    }
}

void OLSQ::addTransformationConstraintsZ3(){
    unsigned_t i, j, t, e, nSpanEdge;
    // Mapping Not Transformations by SWAP Gates.
    // fprintf(stdout, "[Info]          Mapping Not Transformations by SWAP Gates.                       \n");
    BitwuzlaSort *sortbvpi = bitwuzla_mk_bv_sort(_smt.pSolver, ceil(log2(_device.nQubit() + 1)));
    for (t = _olsqParam.swap_duration -1; t < _olsqParam.max_depth - 1; ++t){
        for (i = 0; i < _pCircuit->nProgramQubit(); ++i){
            for (j = 0; j < _device.nQubit(); ++j){
                Qubit& qubit = _device.qubit(j);
                nSpanEdge = qubit.vSpanEdge.size();
                // fprintf(stdout, "[Info]          t = %d, i = %d, j = %d, %d, %d\n", t, i, j, nSpanEdge, qubit.vSpanEdge[0]);
                // expr clause = _smt.vvSigma[t][qubit.vSpanEdge[0]];
                BitwuzlaTerm * clause = _smt.vvSigma[t][qubit.vSpanEdge[0]];
                for (e = 1; e < nSpanEdge; ++e){
                    // clause = clause || _smt.vvSigma[t][qubit.vSpanEdge[e]];
                    clause = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_OR, clause, _smt.vvSigma[t][qubit.vSpanEdge[e]]);
                }
                // clause = !clause && (_smt.vvPi[t][i] == (int_t)j);
                BitwuzlaTerm * qBv = bitwuzla_mk_bv_value_uint64(_smt.pSolver, sortbvpi, j);
                clause = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_AND,
                            bitwuzla_mk_term1(_smt.pSolver, BITWUZLA_KIND_NOT, clause), 
                            bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[t][i], qBv));
                // _smt.smtSolver.add( !clause || _smt.vvPi[t+1][i] == (int_t)j );
                bitwuzla_assert(_smt.pSolver, 
                                bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_OR,
                                    bitwuzla_mk_term1(_smt.pSolver, BITWUZLA_KIND_NOT, clause), 
                                    bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[t+1][i], qBv)));
            }
        }
    }
    // Mapping Transformations by SWAP Gates.
    // fprintf(stdout, "[Info]          Mapping Transformations by SWAP Gates.                       \n");
    for (t = _olsqParam.swap_duration -1; t < _olsqParam.max_depth - 1; ++t){
        for (i = 0; i < _pCircuit->nProgramQubit(); ++i){
            for (e = 0; e < _device.nEdge(); ++e){
                BitwuzlaTerm * q1Bv = bitwuzla_mk_bv_value_uint64(_smt.pSolver, sortbvpi, _device.edge(e).qubitId1());
                BitwuzlaTerm * q2Bv = bitwuzla_mk_bv_value_uint64(_smt.pSolver, sortbvpi, _device.edge(e).qubitId2());
                // _smt.smtSolver.add( !(_smt.vvSigma[t][e] && (_smt.vvPi[t][i] == (int_t)_device.edge(e).qubitId1()))
                //                     || (_smt.vvPi[t+1][i] == (int_t)_device.edge(e).qubitId2()));
                BitwuzlaTerm * cond = bitwuzla_mk_term1(_smt.pSolver, BITWUZLA_KIND_NOT,
                                        bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_AND,
                                            _smt.vvSigma[t][e]), 
                                            bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[t][i], q1Bv));
                bitwuzla_assert(_smt.pSolver, bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_OR,
                                                cond, 
                                                bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[t+1][i], q2Bv)));
                // _smt.smtSolver.add( !(_smt.vvSigma[t][e] && (_smt.vvPi[t][i] == (int_t)_device.edge(e).qubitId2()))
                //                     || (_smt.vvPi[t+1][i] == (int_t)_device.edge(e).qubitId1()));
                cond = bitwuzla_mk_term1(_smt.pSolver, BITWUZLA_KIND_NOT,
                        bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_AND,
                        _smt.vvSigma[t][e]), 
                        bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[t][i], q2Bv));
                bitwuzla_assert(_smt.pSolver, bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_OR,
                                                cond, 
                                                bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[t+1][i], q1Bv)));
            }
        }
    }
}

void OLSQ::addDepthConstraintsZ3(){
    // cout << " add depth constraints" << endl;
    BitwuzlaSort * sortbvtime = bitwuzla_mk_bv_sort(_smt.pSolver, ceil(log2(_olsqParam.max_depth)));
    BitwuzlaTerm * depthBv = bitwuzla_mk_bv_value_uint64(_smt.pSolver, sortbvtime, _olsqParam.min_depth);
    for (BitwuzlaSort* t: _smt.vTg){
        // _smt.smtSolver.add(ult(t, _olsqParam.min_depth));
        bitwuzla_assert(_smt.pSolver, bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_BV_ULT, t, depthBv));
    }
    // cout << "finish add depth constraints" << endl;
}

void OLSQ::addSwapCountConstraintsZ3(unsigned_t bound){
    PB2CNF pb2cnf;
    vector<int_t> vSigmaLit;
    unsigned_t firstFreshVariable = 1 + (_olsqParam.min_depth+1)*_device.nEdge(), newFreashVariable;
    for (unsigned_t i = 1; i < firstFreshVariable; ++i){
        vSigmaLit.emplace_back(i);
    }
    vector< vector<int_t> > formula;
    newFreashVariable = pb2cnf.pb2cnf.encodeAtMostK(vSigmaLit, bound, formula, firstFreshVariable) + 1;
    map<unsigned_t, BitwuzlaTerm*> mAncillary;
    vector<BitwuzlaTerm *> vOrs;
    unsigned_t sigmaT, sigmaE, var;
    BitwuzlaSort *sortbool = bitwuzla_mk_bool_sort(_smt.pSolver);
    string s;
    for(vector<int_t>& clause : formula){
        vOrs.clear();
        for(int& lit : clause){
            var = abs(lit);
            if(var < firstFreshVariable){
                sigmaT = (var - 1);//(1+tight_bound_depth)
                sigmaE = (var - 1) % (1+_olsqParam.min_depth);
                if (lit < 0){
                    vOrs.emplace_back(bitwuzla_mk_term(_smt.pSolver, BITWUZLA_KIND_BV_NOT, _smt.vvSigma[sigmaT][sigmaE]));
                }
                else{
                    vOrs.emplace_back(_smt.vvSigma[sigmaT][sigmaE]);
                }
            }
            else{
                if(mAncillary.find(var) == mAncillary.end()){
                    s = "anc_" + to_string(var);
                    mAncillary[var] = bitwuzla_mk_const(_smt.pSolver, sortbool, s)
                }
                if (lit < 0){
                    vOrs.emplace_back(bitwuzla_mk_term(_smt.pSolver, BITWUZLA_KIND_BV_NOT, mAncillary[var]));
                }
                else{
                    vOrs.emplace_back(mAncillary[var]);
                }
            }
        }    
        BitwuzlaTerm * clause = vOrs[0];
        for (unsigned_t i = 1; i < vOrs.size(); ++i){
            clause = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_OR, vOrs, clause);
        }
        bitwuzla_assert(_smt.pSolver, clause);
    }
}

bool OLSQ::checkModel(){
    BitwuzlaResult status = bitwuzla_check_sat(_smt.pSolver);
    // cout << status << endl;
    // cout << _smt.smtSolver.reason_unknown() << endl;
    // getchar();
    // return true;
    if (status == BITWUZLA_SAT){
        return true;
    }
    else{
        return false;
    }
}

bool OLSQ::optimize(){
    bool success_optimize;
    // cout << "_olsqParam.is_optimize_swap: " << _olsqParam.is_optimize_swap << endl;
    if(!_olsqParam.is_optimize_swap){
        success_optimize = optimizeDepth();
    }
    else{
        success_optimize = optimizeSwap();
    }
    return success_optimize;
}

bool OLSQ::optimizeDepth(){
    bool success, find_min_depth = false, has_jump = false;
    unsigned_t step;
    unsigned_t i;
    if (_olsqParam.is_transition){
        step = 1;
    }
    else{
        step = (_olsqParam.min_depth > 100) ? 10 : 5;
    }
    while(!find_min_depth && _olsqParam.min_depth < _olsqParam.max_depth ){
        fprintf(stdout, "[Info]          trying to optimize for depth bound %d            \n", _olsqParam.min_depth);
        _timer.start(TimeUsage::PARTIAL);
        bitwuzla_push(_smt.pSolver, 1);
        addDepthConstraintsZ3();
        success = checkModel();
        // if (_verbose == 2 ){
        //     fprintf(stdout, "[Info]          SMT solving statistics                         \n");
        //     cout << _smt.smtSolver.statistics() << endl;
        // }
        // string s = Z3_solver_to_string(_smt.c, _smt.smtSolver);
        // cout << s;
        // getchar();
        fprintf(stdout, "[Info]          optimization results: %s                         \n", success ? "success" : "fail");
        _timer.showUsage("optimizing depth", TimeUsage::PARTIAL);
        _timer.showUsage("optimizing depth", TimeUsage::FULL);
        if (success){
            extractModel();
            // getchar();
            bitwuzla_pop(_smt.pSolver, 1);
            for (i = 1; i < step && !find_min_depth && has_jump; ++i){
                fprintf(stdout, "[Info]          trying to optimize for depth bound %d            \n", _olsqParam.min_depth);
                _timer.start(TimeUsage::PARTIAL);
                --_olsqParam.min_depth;
                bitwuzla_push(_smt.pSolver, 1);
                addDepthConstraintsZ3();
                success = checkModel();
                fprintf(stdout, "[Info]          optimization results: %s                         \n", success ? "success" : "fail");
                _timer.showUsage("optimizing depth", TimeUsage::PARTIAL);
                _timer.showUsage("optimizing depth", TimeUsage::FULL);
                if (success) 
                    extractModel();
                else
                    find_min_depth = true;
                bitwuzla_pop(_smt.pSolver, 1);
            }
            find_min_depth = true;
        }
        else{
            bitwuzla_pop(_smt.pSolver, 1);
            _olsqParam.min_depth += step;
            has_jump = true;
        }
    }
    if (find_min_depth){
        return true;
    }
    else{
        return false;
    }
}

bool OLSQ::optimizeSwap(){
    bool success_optimize;
    success_optimize = optimizeDepth();
    // getchar();
    if (!success_optimize){
        return false;
    }

    unsigned_t lower_swap_bound = (_olsqParam.is_transition) ? _olsqParam.min_depth : 0;
    unsigned_t upper_swap_bound = (_olsqParam.is_use_SABRE_for_swap) ? _olsqParam.sabre_swap_bound : _pCircuit->nGate();
    upper_swap_bound = (_pCircuit->nSwapGate() < upper_swap_bound) ? _pCircuit->nSwapGate() : upper_swap_bound;
    // cout << "swap lower bound = " << lower_swap_bound << " , swap upper bound = " << upper_swap_bound << endl;
    // getchar();
    bool reduce_swap = true;
    // getchar();
    // cout << "lower_swap_bound: " << lower_swap_bound << endl;
    // cout << "upper_swap_bound: " << upper_swap_bound << endl;
    // cout << reduce_swap << endl;
    // cout << _timer.fullRealTime() << endl;
    // cout << _olsqParam.timeout << endl;
    bool firstRun = true;
    while (reduce_swap && _timer.fullRealTime() < _olsqParam.timeout){
        // cout << "enter loop" << endl;
        addDepthConstraintsZ3();
        reduce_swap = optimizeSwapForDepth(lower_swap_bound, upper_swap_bound, firstRun);
        if (_olsqParam.is_transition){
            _olsqParam.min_depth = _olsqParam.min_depth + 5;
            _olsqParam.max_depth = _olsqParam.min_depth + 6;
        }
        else{
            _olsqParam.min_depth = _olsqParam.min_depth + 10;
            _olsqParam.max_depth = _olsqParam.min_depth + 11;
        }
        upper_swap_bound = _pCircuit->nSwapGate() - 1;
        firstRun = false;
        // getchar();
        if(reduce_swap){
            fprintf(stdout, "[Info] Successfully reduce SWAP count. Go to next run.            \n");
            fprintf(stdout, "[Info] Solving with depth %d            \n", _olsqParam.min_depth);
            _timer.start(TimeUsage::PARTIAL);
            fprintf(stdout, "[Info] Generating formulation                        \n");
            generateFormulationZ3();
            _timer.showUsage("Generating formulation", TimeUsage::PARTIAL);
            // _timer.showUsage("Generating formulation", TimeUsage::FULL);
            _timer.start(TimeUsage::PARTIAL);
        }
    }
    return true;
}

bool OLSQ::optimizeSwapForDepth(unsigned_t lower_swap_bound, unsigned_t upper_swap_bound, bool firstRun){
    unsigned_t swap_bound = upper_swap_bound;
    bool find_min_swap = false;
    bool success;
    // if (upper_swap_bound <= lower_swap_bound)
    //     find_min_swap = true;
    while (!find_min_swap && lower_swap_bound <= swap_bound && swap_bound <= upper_swap_bound && _timer.fullRealTime() < _olsqParam.timeout){
        fprintf(stdout, "[Info]          trying to optimize for swap bound %d            \n", swap_bound);
        _timer.start(TimeUsage::PARTIAL);
        bitwuzla_push(_smt.pSolver, 1);
        addSwapCountConstraintsZ3(swap_bound);
        // string s = Z3_solver_to_string(_smt.c, _smt.smtSolver);
        // cout << s;
        // getchar();
        success = checkModel();
        if (_verbose == 2 ){
            fprintf(stdout, "[Info]          SMT solving statistics                         \n");
            cout << _smt.smtSolver.statistics() << endl;
        }
        fprintf(stdout, "[Info]          optimization results: %s                         \n", success ? "success" : "fail");
        _timer.showUsage("optimizing swap", TimeUsage::PARTIAL);
        _timer.showUsage("optimizing swap", TimeUsage::FULL);
        if (success){
            extractModel();
            --swap_bound;
        }
        else{
            // cout << "First run " << firstRun << endl;
            // cout << "upper_swap_bound " << upper_swap_bound << endl;
            // cout << "swap_bound " << swap_bound << endl;
            if (swap_bound < upper_swap_bound){
                find_min_swap = true;
            }
            else if (!firstRun && swap_bound == upper_swap_bound){
                break;
            }
            else if (firstRun && swap_bound == upper_swap_bound){
                // cout << "line 518: increase swap bound" << endl;
                ++upper_swap_bound;
                lower_swap_bound = upper_swap_bound;
                swap_bound = upper_swap_bound;
            }
        }
        bitwuzla_pop(_smt.pSolver, 1);
    }
    return find_min_swap;
}

void OLSQ::extractModel(){
    fprintf(stdout, "[Info] Extract Model Info                              \n");
    // model m = _smt.smtSolver.get_model();
    unsigned_t circuitDepth = 0, i, gateTime, q, j, swapId, qId, t, e;
    // collect gate execution time
    string s;
    for ( i = 0; (i < _pCircuit->nGate());  ++i){
        Gate &gate = _pCircuit->gate(i);
        const char *rstr = bitwuzla_get_bv_value(_smt.pSolver, _smt.vTg[i]);
        s = rstr;
        gateTime = stoi(s, nullptr, 2); 
        circuitDepth = (circuitDepth < gateTime) ? gateTime : circuitDepth;
        gate.setExecutionTime(gateTime);
        for ( j = 0; j < gate.nTargetQubit();  ++j ){
            const char *rstr = bitwuzla_get_bv_value(_smt.pSolver, _smt.vvPi[gateTime][gate.targetProgramQubit(j)]);
            s = rstr;
            gate.addTargetPhysicalQubit(j, stoi(s, nullptr, 2));
            // gate.addTargetPhysicalQubit(j, m.eval(_smt.vvPi[gateTime][gate.targetProgramQubit(j)]).get_numeral_int64());
            // cout << m.eval(_smt.vvPi[gateTime][q]) << " " << m.eval(_smt.vvPi[gateTime][q]).get_numeral_int64() << endl;
        }
        if (_verbose > 0){
            fprintf(stdout, "        - Gate %d, name: %s, duration: %d, time: %d, target program qubit: ", i, gate.name().c_str(), gate.duration(), gate.executionTime());
            for ( j = 0; j < gate.nTargetQubit();  ++j ){
                fprintf(stdout, "%d ", gate.targetProgramQubit(j));
            }
            fprintf(stdout, ", target physical qubit: ");
            for ( j = 0; j < gate.nTargetQubit();  ++j ){
                fprintf(stdout, "%d ", gate.targetPhysicalQubit(j));
            }
            fprintf(stdout, "\n");
        }
    }
    // get SWAP gate
    _pCircuit->clearSwap();
    vector<unsigned_t> swapTargetQubit(2,0);
    for (t = _olsqParam.swap_duration -1; t < circuitDepth; ++t){
        for (e = 0; e < _device.nEdge(); ++e){
            // cout << "t = " << t << ", e = " << e << " " << m.eval(_smt.vvSigma[t][e]).bool_value() << endl;
            const char *rstr = bitwuzla_get_bv_value(_smt.pSolver, _smt.vTg[i]);
            s = rstr;
            // if (m.eval(_smt.vvSigma[t][e]).bool_value() == 1){
            if (stoi(s, nullptr, 2)){
                swapTargetQubit[0] = _device.edge(e).qubitId1();
                swapTargetQubit[1] = _device.edge(e).qubitId2();
                _pCircuit->addSwapGate(swapTargetQubit, _olsqParam.swap_duration);
                swapId = _pCircuit->nSwapGate() - 1;
                Gate & gate = _pCircuit->swapGate(swapId);
                gate.setExecutionTime(t);
                if (_verbose > 0){
                    fprintf(stdout, "        - SWAP Gate %d, duration: %d, time: %d, target qubit: %d %d\n", swapId + _pCircuit->nGate(), gate.duration(), gate.executionTime(), gate.targetPhysicalQubit(0), gate.targetPhysicalQubit(1));
                }
            }
        }
    }
    // set initial and final mapping
    if (_verbose > 0){
        fprintf(stdout, "        - Qubit mapping: \n");
        for (t = 0; t <= circuitDepth; ++t){
        fprintf(stdout, "        - Time %d: ",t);
            for (i = 0; i < _pCircuit->nProgramQubit(); ++i){
                fprintf(stdout, "%d->%d ", i, (int_t)m.eval(_smt.vvPi[t][i]).get_numeral_int64());
            }
        fprintf(stdout, "\n");
        }
    }
    
    for (i = 0; i < _pCircuit->nProgramQubit(); ++i){
        const char *rstr1 = bitwuzla_get_bv_value(_smt.pSolver, _smt.vvPi[0][i]);
        const char *rstr2 = bitwuzla_get_bv_value(_smt.pSolver, _smt.vvPi[circuitDepth][i]);
        s = rstr1;
        _pCircuit->setInitialMapping(i,  stoi(s, nullptr, 2));
        s = rstr2;
        _pCircuit->setInitialMapping(i,  stoi(s, nullptr, 2));
        // _pCircuit->setInitialMapping(i, m.eval(_smt.vvPi[0][i]).get_numeral_int64());
        // _pCircuit->setFinalMapping(i, m.eval(_smt.vvPi[circuitDepth][i]).get_numeral_int64());
    }
    _pCircuit->setCircuitDepth(circuitDepth+1);
}

void OLSQ::asapScheduling(){
    vector<int_t> vPushForwardDepth(_device.nQubit(), -1);
    int_t gateExecutionTime;
    unsigned_t block, i, j, q0, q1, qId, maxTime = 0;
    Gate gate;
    set<unsigned_t> sGateId;
    for (block = 0; block < _pCircuit->circuitDepth(); ++block){
        // cout << "block " << block << endl;
        for (i = 0; i < _pCircuit->nGate(); ++i){
            Gate & gate = _pCircuit->gate(i);
            if (gate.executionTime() == block && sGateId.count(gate.idx()) == 0){
                // cout << "gate " << gate.idx() << endl;
                // cout << "vPushForwardDepth[0] " << vPushForwardDepth[gate.targetPhysicalQubit(0)] << endl;
                // cout << "vPushForwardDepth[1] " << vPushForwardDepth[gate.targetPhysicalQubit(1)] << endl;
                gateExecutionTime = vPushForwardDepth[gate.targetPhysicalQubit(0)];
                for (j = 1; j < gate.nTargetQubit(); ++j){
                    qId = gate.targetPhysicalQubit(j);
                    gateExecutionTime = (gateExecutionTime < vPushForwardDepth[qId]) ? vPushForwardDepth[qId] : gateExecutionTime;
                    // cout << "gate' execution time " << gateExecutionTime << endl;
                }
                ++gateExecutionTime;
                for (j = 0; j < gate.nTargetQubit(); ++j){
                    qId = gate.targetPhysicalQubit(j);
                    vPushForwardDepth[qId] = gateExecutionTime;
                }
                // cout << "gate execution time " << gateExecutionTime << endl;
                // cout << "vPushForwardDepth[0] " << vPushForwardDepth[gate.targetPhysicalQubit(0)] << endl;
                // cout << "vPushForwardDepth[1] " << vPushForwardDepth[gate.targetPhysicalQubit(1)] << endl;
                maxTime = (maxTime < gateExecutionTime) ? gateExecutionTime : maxTime;
                gate.setExecutionTime(gateExecutionTime);
                sGateId.insert(gate.idx());
            }
        }
        if (block < _pCircuit->circuitDepth() - 1){
            for (j = 0; j < _pCircuit->nSwapGate(); ++j){
                Gate & gate = _pCircuit->swapGate(j);
                if (gate.executionTime() == block && sGateId.count(gate.idx()) == 0){
                    // cout << "swap gate " << gate.idx() << endl;
                    // cout << "vPushForwardDepth[0] " << vPushForwardDepth[gate.targetPhysicalQubit(0)] << endl;
                    // cout << "vPushForwardDepth[1] " << vPushForwardDepth[gate.targetPhysicalQubit(1)] << endl;
                    q0 = gate.targetPhysicalQubit(0);
                    q1 = gate.targetPhysicalQubit(1);
                    gateExecutionTime = (vPushForwardDepth[q0] < vPushForwardDepth[q1]) ? vPushForwardDepth[q1] : vPushForwardDepth[q0];
                    ++gateExecutionTime;
                    vPushForwardDepth[q0] = gateExecutionTime;
                    vPushForwardDepth[q1] = gateExecutionTime;
                    maxTime = (maxTime < gateExecutionTime) ? gateExecutionTime : maxTime;
                    gate.setExecutionTime(gateExecutionTime);
                    sGateId.insert(gate.idx());
                }
            }
        }
    }
    ++maxTime;
    _pCircuit->setCircuitDepth(maxTime);
}

void OLSQ::outputQASM(string const & fileName){
    vector< vector<unsigned_t> > vvTimeGate(_pCircuit->circuitDepth(), vector<unsigned_t>());
    Gate gate;
    unsigned_t i, t, qId, j;
    string line;
    for (i = 0; i < _pCircuit->nGate(); ++i){
        Gate & gate = _pCircuit->gate(i);
        vvTimeGate[gate.executionTime()].emplace_back(i);
    }
        for (i = 0; i < _pCircuit->nSwapGate(); ++i){
        Gate & gate = _pCircuit->swapGate(i);
        vvTimeGate[gate.executionTime()].emplace_back(i);
    }
    line = "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[" + to_string(_device.nQubit()) + "];\ncreg c[" + to_string(_device.nQubit()) + "];\n";
    for (t = 0; t < _pCircuit->circuitDepth(); ++t){
        for (i = 0; i < vvTimeGate[t].size(); ++i){
            if (vvTimeGate[t][i] < _pCircuit->nGate()){
                gate = _pCircuit->gate(vvTimeGate[t][i]);
            }
            else{
                gate = _pCircuit->swapGate(vvTimeGate[t][i] - _pCircuit->nGate());
            }
            
            if (gate.nTargetQubit() == 1){
                line = line + gate.name() + " q[" + to_string(gate.targetPhysicalQubit(0)) + "];\n";
            }
            else{
                line = line + gate.name() + " q[" + to_string(gate.targetPhysicalQubit(0)) + "], q[" + to_string(gate.targetPhysicalQubit(1)) + "];\n";
            }
        }
    }
    FILE* fout = fopen(fileName.c_str(), "w");
    fprintf(fout, "%s", line.c_str());
    fclose(fout);
}

void OLSQ::increase_depth_bound(){
    _olsqParam.max_depth *= _olsqParam.max_depth_expand_factor; 
}

unsigned_t OLSQ::extract_longest_chain(){
    vector<unsigned_t> vPushForwardDepth(_pCircuit->nProgramQubit(),0);
    Gate gate;
    unsigned_t max, qId, i, j;
    for (i = 0; i < _pCircuit->nGate(); ++i){
        Gate & gate = _pCircuit->gate(i);
        if (gate.nTargetQubit() == 1){
            ++vPushForwardDepth[gate.targetProgramQubit(0)];
        }
        else{
            max = 0;
            for ( j = 0; j < gate.nTargetQubit();  ++j ){
                qId = gate.targetProgramQubit(j);
                max = (max > vPushForwardDepth[qId]) ? max : vPushForwardDepth[qId];
            }
            ++max;
            for ( j = 0; j < gate.nTargetQubit();  ++j ){
                vPushForwardDepth[gate.targetProgramQubit(j)] = max;
            }
        }
    }
    return *std::max_element(vPushForwardDepth.begin(), vPushForwardDepth.end());
}

void OLSQ::constructDependency(){
    _vpGateDependency.clear();
    vector<int_t> vQubitLastGate(_pCircuit->nProgramQubit(), -1);
    Gate gate;
    unsigned_t qId, i, j;
    for (i = 0; i < _pCircuit->nGate(); ++i){
        Gate & gate = _pCircuit->gate(i);
        for ( j = 0; j < gate.nTargetQubit();  ++j ){
                qId = gate.targetProgramQubit(j);
            if (vQubitLastGate[qId] > -1){
                addDependency(vQubitLastGate[qId], i);
            }
            vQubitLastGate[qId] = i;
        }
    }
}

void OLSQ::printDependency(){
    pair<unsigned_t, unsigned_t> p;
    unsigned_t i;
    fprintf(stdout, "[Info] Circuit Dependency Info                              \n");
    fprintf(stdout, "       ------------------------------------------\n");
    fprintf(stdout, "       - Dependency list\n");
    for ( i = 0; i < _vpGateDependency.size();  ++i ){
        fprintf(stdout, "        - Gate %d depends on gate %d\n", _pCircuit->gate(_vpGateDependency[i].second).idx(), _pCircuit->gate(_vpGateDependency[i].first).idx());
    }
}

OLSQ_NAMESPACE_CPP_END
