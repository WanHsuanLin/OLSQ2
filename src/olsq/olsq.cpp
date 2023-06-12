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

// TODO: add restricted qubit region: (1) if keeping unsat, relax region. (2) transforming constraint. (3) swap constraint,

OLSQ_NAMESPACE_CPP_START

void OLSQ::setDependency(vector<pair<unsigned_t, unsigned_t> >& vDependencies){
    for (pair<unsigned_t, unsigned_t>& p : vDependencies){
        assert(isValidGateIdx(p.first));
        assert(isValidGateIdx(p.second));
        addDependency(p.first, p.second);
    }
}

void OLSQ::run(string const & fileName){
    _outputQASMFile = fileName;
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
            increaseDepthBound();
            fprintf(stdout, "[Info] Longest chain = %d\n", _olsqParam.min_depth);
        }
    }
    bool solve = false;
    unsigned_t iter = 0;
    if (!_olsqParam.is_transition && _olsqParam.use_window_range_for_gate){
        constructGateTimeWindow();
    }
    if (_olsqParam.is_given_mapping_region){
        collectQubitRegion();
    }
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
        if(!solve){
            increaseDepthBound();
            _smt.reset();
        }
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
    fprintf(stdout, "[Info]          constructing swap overlapping constraint     \n");
    addSwapConstraintsZ3();
    fprintf(stdout, "[Info]          constructing mapping transformation constraint\n");
    addTransformationConstraintsZ3();
}

void OLSQ::constructVariableZ3(){
    unsigned_t bit_length_pi, bit_length_time;
    bit_length_pi = ceil(log2(_device.nQubit() + 1));
    bit_length_time = _olsqParam.max_depth_bit;
    _smt.vvPi.reserve(_olsqParam.max_depth);
    _smt.vTg.reserve(_pCircuit->nGate());
    _smt.vvSigma.reserve(_olsqParam.max_depth);

    string s;
    unsigned_t i, j;

    // Create a bit-vector sort of size 1.
    const BitwuzlaSort *sortbvpi = bitwuzla_mk_bv_sort(_smt.pSolver, bit_length_pi);
    const BitwuzlaSort *sortbvtime = bitwuzla_mk_bv_sort(_smt.pSolver, bit_length_time);
    const BitwuzlaSort *sortbool = bitwuzla_mk_bool_sort(_smt.pSolver);

    for (i = 0; i < _olsqParam.max_depth; ++i){
        _smt.vvPi.emplace_back(vector<const BitwuzlaTerm*> ());
        for (j = 0; j < _pCircuit->nProgramQubit(); ++j){
            s = "map_t" + to_string(i) + "_q" + to_string(j);
            _smt.vvPi[i].emplace_back(bitwuzla_mk_const(_smt.pSolver, sortbvpi, s.c_str()));
        }
    }

    for (i = 0; i < _pCircuit->nGate(); ++i){
        s = "time_" + to_string(i);
        _smt.vTg.emplace_back(bitwuzla_mk_const(_smt.pSolver, sortbvtime, s.c_str()));
    }

    for (i = 0; i < _olsqParam.max_depth; ++i){
        _smt.vvSigma.emplace_back(vector<const BitwuzlaTerm*> ());
        for (j = 0; j < _device.nEdge(); ++j){
            s = "ifswap_t" + to_string(i) + "_e" + to_string(j);
            _smt.vvSigma[i].emplace_back(bitwuzla_mk_const(_smt.pSolver, sortbool, s.c_str()));
        }
    }

}

void OLSQ::addInjectiveMappingConstraintsZ3(unsigned_t boundOffset){
    unsigned_t i, j, k, qId;
    unsigned_t end = _olsqParam.min_depth, begin = 0;
    if(boundOffset > 0){
        begin = end;
        end += boundOffset;
    }
    int_t pId;
    const BitwuzlaSort *sortbvpi = bitwuzla_mk_bv_sort(_smt.pSolver, ceil(log2(_device.nQubit() + 1)));
    const BitwuzlaTerm * zero = bitwuzla_mk_bv_zero(_smt.pSolver, sortbvpi);
    const BitwuzlaTerm * nqubit = bitwuzla_mk_bv_value_uint64(_smt.pSolver, sortbvpi, _device.nQubit());
    vector<set<int_t>> * pvsQubitRegion = _pCircuit->pvQubitRegion();
    set<int>::iterator it;
    if(_olsqParam.is_given_mapping_region){
        for (j = 0; j < _pCircuit->nProgramQubit(); ++j){
            for (i = begin; i < end; ++i){
                it = (*pvsQubitRegion)[j].begin();
                const BitwuzlaTerm * bvp = bitwuzla_mk_bv_value_uint64(_smt.pSolver, sortbvpi, (*it));
                const BitwuzlaTerm * clause = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[i][j], bvp);
                ++it;
                for (; it != (*pvsQubitRegion)[j].end(); ++it) {
                    bvp = bitwuzla_mk_bv_value_uint64(_smt.pSolver, sortbvpi, (*it));
                    clause = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_OR, 
                                clause, 
                                bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[i][j], bvp));
                }
                bitwuzla_assert(_smt.pSolver, clause);
            }
            for (k = 0; k < j; ++k){
                set<int> intersectSet;
                set_intersection((*pvsQubitRegion)[j].begin(), (*pvsQubitRegion)[j].end(), (*pvsQubitRegion)[k].begin(), (*pvsQubitRegion)[k].end(),
                            std::inserter(intersectSet, intersectSet.begin()));
                if(intersectSet.size() > 0){
                    for (i = begin; i < end; ++i){
                        bitwuzla_assert(_smt.pSolver, bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_DISTINCT, _smt.vvPi[i][j], _smt.vvPi[i][k]));
                    }
                }
            }
        }
    }
    else{
        for (i = begin; i < end; ++i){
            for (j = 0; j < _pCircuit->nProgramQubit(); ++j){
                bitwuzla_assert(_smt.pSolver, bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_BV_ULE, zero, _smt.vvPi[i][j]));
                bitwuzla_assert(_smt.pSolver, bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_BV_ULT, _smt.vvPi[i][j], nqubit));
                for (k = 0; k < j; ++k){
                    bitwuzla_assert(_smt.pSolver, bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_DISTINCT, _smt.vvPi[i][j], _smt.vvPi[i][k]));
                }
            }
        }
    }
}

void OLSQ::addValidTwoQubitGateConstraintsZ3(unsigned_t boundOffset){
    Gate gate;
    Edge edge;
    unsigned_t i, t, j;
    vector<set<int_t>> * pvsQubitRegion = _pCircuit->pvQubitRegion();
    unsigned_t end = _olsqParam.min_depth, begin = 0;
    const BitwuzlaSort *sortbvtime = bitwuzla_mk_bv_sort(_smt.pSolver, _olsqParam.max_depth_bit);
    const BitwuzlaSort *sortbvpi = bitwuzla_mk_bv_sort(_smt.pSolver, ceil(log2(_device.nQubit() + 1)));
    if(boundOffset > 0){
        begin = end;
        end += boundOffset;
    }
    for ( i = 0; i < _pCircuit->nGate();  ++i ){
        Gate & gate = _pCircuit->gate(i);
        if ((gate).nTargetQubit() == 2){
            if(_olsqParam.use_window_range_for_gate){
                begin = _vpGateTimeWindow[i].first;
                end = _vpGateTimeWindow[i].second + 1;
                if(boundOffset > 0){
                    begin = end;
                    end += boundOffset;
                }
                cerr << "construct valid two qubit constraint for gate "<< i <<" from " << begin << " to " << end << endl;
            }
            for (t = begin; t < end; ++t){
                const BitwuzlaTerm * bvt = bitwuzla_mk_bv_value_uint64(_smt.pSolver, sortbvtime, t);
                if(_olsqParam.is_given_mapping_region){
                    const BitwuzlaTerm * clause = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vTg[i], bvt);
                    const BitwuzlaTerm *pro1EqPhy1, *pro2EqPhy2, *cond;
                    for (int_t j : (*pvsQubitRegion)[gate.targetProgramQubit(0)]){
                        for (int_t k : (*pvsQubitRegion)[gate.targetProgramQubit(1)]){
                            if(_device.isAdjacent(j, k)){
                                const BitwuzlaTerm * q1Bv = bitwuzla_mk_bv_value_uint64(_smt.pSolver, sortbvpi, j);
                                const BitwuzlaTerm * q2Bv = bitwuzla_mk_bv_value_uint64(_smt.pSolver, sortbvpi, k);
                                pro1EqPhy1 = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[t][gate.targetProgramQubit(0)], q1Bv);
                                pro2EqPhy2 = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[t][gate.targetProgramQubit(1)], q2Bv);
                                cond = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_AND, pro1EqPhy1, pro2EqPhy2);
                                clause = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_OR, clause, cond);
                            }
                        }
                    }
                    bitwuzla_assert(_smt.pSolver, clause);
                }
                else{
                    const BitwuzlaTerm * clause = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vTg[i], bvt);
                    clause = bitwuzla_mk_term1(_smt.pSolver, BITWUZLA_KIND_NOT, clause);
                    const BitwuzlaTerm *pro1EqPhy1, *pro2EqPhy2, *pro1EqPhy2, *pro2EqPhy1, *cond1, *cond2;
                    for ( j = 0; j < _device.nEdge();  ++j ){
                        Edge & edge = _device.edge(j);
                        const BitwuzlaTerm * q1Bv = bitwuzla_mk_bv_value_uint64(_smt.pSolver, sortbvpi, edge.qubitId1());
                        const BitwuzlaTerm * q2Bv = bitwuzla_mk_bv_value_uint64(_smt.pSolver, sortbvpi, edge.qubitId2());
                        pro1EqPhy1 = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[t][gate.targetProgramQubit(0)], q1Bv);
                        pro2EqPhy2 = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[t][gate.targetProgramQubit(1)], q2Bv);
                        pro1EqPhy2 = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[t][gate.targetProgramQubit(0)], q2Bv);
                        pro2EqPhy1 = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[t][gate.targetProgramQubit(1)], q1Bv);
                        cond1 = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_AND, pro1EqPhy1, pro2EqPhy2);
                        cond2 = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_AND, pro2EqPhy1, pro1EqPhy2);
                        clause = bitwuzla_mk_term3(_smt.pSolver, BITWUZLA_KIND_OR, clause, cond1, cond2);
                    }                       
                    bitwuzla_assert(_smt.pSolver, clause);
                }
            }
        }
    }    
}

void OLSQ::addDependencyConstraintsZ3(){
    unsigned_t i;
    Gate g;
    const BitwuzlaSort * sortbvtime = bitwuzla_mk_bv_sort(_smt.pSolver, _olsqParam.max_depth_bit);
    const BitwuzlaTerm * zero = bitwuzla_mk_bv_zero(_smt.pSolver, sortbvtime);
    for (i = 0; i < _pCircuit->nGate(); ++i){
        bitwuzla_assert(_smt.pSolver, bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_BV_ULE, zero, _smt.vTg[i]));
    }
    if (_olsqParam.is_transition){
        for ( i = 0; i < _vpGateDependency.size();  ++i ){
            bitwuzla_assert(_smt.pSolver, bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_BV_ULE, _smt.vTg[_vpGateDependency[i].first], _smt.vTg[_vpGateDependency[i].second]));
        }
    }
    else{
        for ( i = 0; i < _vpGateDependency.size();  ++i ){
            bitwuzla_assert(_smt.pSolver, bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_BV_ULT, _smt.vTg[_vpGateDependency[i].first], _smt.vTg[_vpGateDependency[i].second]));
        }
    }
}

void OLSQ::addSwapConstraintsZ3(unsigned_t boundOffset){
    // No swap for t<s
    unsigned_t i, j, t, e, tt, q1, q2;
    Qubit qubit;
    unsigned_t end = _olsqParam.min_depth, begin = 0;
    unsigned_t bound = (_olsqParam.swap_duration < _olsqParam.max_depth) ? _olsqParam.swap_duration: _olsqParam.max_depth;
    if(_olsqParam.is_transition){
        --bound;
    }
    // unsigned_t begin = _olsqParam.swap_duration -1, end = _olsqParam.min_depth;
    if(boundOffset > 0){
        begin = end;
        end += boundOffset;
    }
    for (i = 0; i < bound; ++i){
        for (j = 0; j < _device.nEdge(); ++j){
            bitwuzla_assert(_smt.pSolver, bitwuzla_mk_term1(_smt.pSolver, BITWUZLA_KIND_NOT, _smt.vvSigma[i][j]));
        }
    }
    // cout << "begin: " << begin << ", end: " << end << endl;
    // swap gates can not overlap with swap in space
    for (t = begin; t < end; ++t){
        for (e = 0; e < _device.nEdge(); ++e){
            q1 = _device.edge(e).qubitId1();
            q2 = _device.edge(e).qubitId2();
            for (unsigned_t ee : _device.qubit(q1).vSpanEdge){
                for (tt = t - _olsqParam.swap_duration + 1; tt < t + 1; ++tt){
                    if (ee < e){
                        bitwuzla_assert(_smt.pSolver, 
                            bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_OR, 
                                bitwuzla_mk_term1(_smt.pSolver, BITWUZLA_KIND_NOT, _smt.vvSigma[t][e]),
                                bitwuzla_mk_term1(_smt.pSolver, BITWUZLA_KIND_NOT, _smt.vvSigma[tt][ee])));
                    }
                }
            }
            for (unsigned_t ee : _device.qubit(q2).vSpanEdge){
                for (tt = t - _olsqParam.swap_duration + 1; tt < t + 1; ++tt){
                    if (ee < e){
                        bitwuzla_assert(_smt.pSolver, 
                            bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_OR, 
                                bitwuzla_mk_term1(_smt.pSolver, BITWUZLA_KIND_NOT, _smt.vvSigma[t][e]),
                                bitwuzla_mk_term1(_smt.pSolver, BITWUZLA_KIND_NOT, _smt.vvSigma[tt][ee])));
                    }
                }
            }
        }
    }
    if (!_olsqParam.is_transition){
        begin = 0;
        end = _olsqParam.min_depth;
        const BitwuzlaSort *sortbvtime = bitwuzla_mk_bv_sort(_smt.pSolver, _olsqParam.max_depth_bit);
        const BitwuzlaSort *sortbvpi = bitwuzla_mk_bv_sort(_smt.pSolver, ceil(log2(_device.nQubit() + 1)));
        // swap gates can not overlap with swap in time
        for (t = begin; t < end; ++t){
            for (e = 0; e < _device.nEdge(); ++e){
                for (tt = t - _olsqParam.swap_duration + 1; tt < t; ++tt){
                    bitwuzla_assert(_smt.pSolver, 
                            bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_OR, 
                                bitwuzla_mk_term1(_smt.pSolver, BITWUZLA_KIND_NOT, _smt.vvSigma[t][e]),
                                bitwuzla_mk_term1(_smt.pSolver, BITWUZLA_KIND_NOT, _smt.vvSigma[tt][e])));
                }
            }
        }
        // swap gates can not ovelap with other gates
        // the time interval should be modified
        for (i = 0; i < _pCircuit->nGate(); ++i){
            if(_olsqParam.use_window_range_for_gate){
                begin = _vpGateTimeWindow[i].first;
                end = _vpGateTimeWindow[i].second + 1;
                if(boundOffset > 0){
                    begin = end;
                    end += boundOffset;
                }
            }
            for (t = begin; t < end; ++t){
                for (e = 0; e < _device.nEdge(); ++e){
                    Gate& gate = _pCircuit->gate(i);
                    for (tt = t; tt < t + _olsqParam.swap_duration; ++tt){
                        const BitwuzlaTerm * bvt = bitwuzla_mk_bv_value_uint64(_smt.pSolver, sortbvtime, tt);
                        const BitwuzlaTerm * q1Bv = bitwuzla_mk_bv_value_uint64(_smt.pSolver, sortbvpi, _device.edge(e).qubitId1());
                        const BitwuzlaTerm * q2Bv = bitwuzla_mk_bv_value_uint64(_smt.pSolver, sortbvpi, _device.edge(e).qubitId2());
                        const BitwuzlaTerm * clause = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vTg[i], bvt);
                        const BitwuzlaTerm * clauseOr = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_OR, 
                                                    bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[tt][gate.targetProgramQubit(0)], q1Bv), //pro1EqPhy1,
                                                    bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[tt][gate.targetProgramQubit(0)], q2Bv));
                        const BitwuzlaTerm * clauseAnd;
                        const BitwuzlaTerm *cond;
                        cond = bitwuzla_mk_term1(_smt.pSolver, BITWUZLA_KIND_NOT, _smt.vvSigma[t][e]);
                        if (gate.nTargetQubit() == 1){
                            clauseAnd =bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_AND, clause, clauseOr);
                            clauseAnd = bitwuzla_mk_term1(_smt.pSolver, BITWUZLA_KIND_NOT, clauseAnd);
                            bitwuzla_assert(_smt.pSolver, 
                                bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_OR, clauseAnd, cond));
                        }
                        else if(gate.nTargetQubit() == 2){
                            clauseOr = bitwuzla_mk_term3(_smt.pSolver, BITWUZLA_KIND_OR, 
                                                    bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[tt][gate.targetProgramQubit(1)], q1Bv), // pro2EqPhy1,
                                                    bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[tt][gate.targetProgramQubit(1)], q2Bv), // pro2EqPhy2,
                                                    clauseOr);
                            clauseAnd =bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_AND, clause, clauseOr);
                            clauseAnd = bitwuzla_mk_term1(_smt.pSolver, BITWUZLA_KIND_NOT, clauseAnd);
                            bitwuzla_assert(_smt.pSolver, 
                                bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_OR, clauseAnd, cond));
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

void OLSQ::addTransformationConstraintsZ3(unsigned_t boundOffset){
    unsigned_t i, j, t, e, nSpanEdge;
    unsigned_t end = _olsqParam.min_depth, begin = 0;
    if(boundOffset > 0){
        begin = end;
        end += boundOffset;
    }
    // Mapping Not Transformations by SWAP Gates.
    // fprintf(stdout, "[Info]          Mapping Not Transformations by SWAP Gates.                       \n");
    const BitwuzlaSort *sortbvpi = bitwuzla_mk_bv_sort(_smt.pSolver, ceil(log2(_device.nQubit() + 1)));
    for (t = begin; t < end; ++t){
        for (i = 0; i < _pCircuit->nProgramQubit(); ++i){
            for (j = 0; j < _device.nQubit(); ++j){
                Qubit& qubit = _device.qubit(j);
                nSpanEdge = qubit.vSpanEdge.size();
                const BitwuzlaTerm * clause = _smt.vvSigma[t][qubit.vSpanEdge[0]];
                for (e = 1; e < nSpanEdge; ++e){
                    clause = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_OR, clause, _smt.vvSigma[t][qubit.vSpanEdge[e]]);
                }
                const BitwuzlaTerm * qBv = bitwuzla_mk_bv_value_uint64(_smt.pSolver, sortbvpi, j);
                clause = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_AND,
                            bitwuzla_mk_term1(_smt.pSolver, BITWUZLA_KIND_NOT, clause), 
                            bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[t][i], qBv));
                bitwuzla_assert(_smt.pSolver, 
                                bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_OR,
                                    bitwuzla_mk_term1(_smt.pSolver, BITWUZLA_KIND_NOT, clause), 
                                    bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[t+1][i], qBv)));
            }
        }
    }
    // Mapping Transformations by SWAP Gates.
    // fprintf(stdout, "[Info]          Mapping Transformations by SWAP Gates.                       \n");
    for (t = begin; t < end; ++t){
        for (i = 0; i < _pCircuit->nProgramQubit(); ++i){
            for (e = 0; e < _device.nEdge(); ++e){
                const BitwuzlaTerm * q1Bv = bitwuzla_mk_bv_value_uint64(_smt.pSolver, sortbvpi, _device.edge(e).qubitId1());
                const BitwuzlaTerm * q2Bv = bitwuzla_mk_bv_value_uint64(_smt.pSolver, sortbvpi, _device.edge(e).qubitId2());
                const BitwuzlaTerm * cond = bitwuzla_mk_term1(_smt.pSolver, BITWUZLA_KIND_NOT,
                                        bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_AND,
                                            _smt.vvSigma[t][e], 
                                            bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[t][i], q1Bv)));
                bitwuzla_assert(_smt.pSolver, bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_OR,
                                                cond, 
                                                bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[t+1][i], q2Bv)));
                cond = bitwuzla_mk_term1(_smt.pSolver, BITWUZLA_KIND_NOT,
                        bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_AND,
                        _smt.vvSigma[t][e], 
                        bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[t][i], q2Bv)));
                bitwuzla_assert(_smt.pSolver, bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_OR,
                                                cond, 
                                                bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_EQUAL, _smt.vvPi[t+1][i], q1Bv)));
            }
        }
    }
}

void OLSQ::addDepthConstraintsZ3(){
    const BitwuzlaSort * sortbvtime = bitwuzla_mk_bv_sort(_smt.pSolver, _olsqParam.max_depth_bit);
    if(!_olsqParam.use_window_range_for_gate){
        const BitwuzlaTerm * depthBv = bitwuzla_mk_bv_value_uint64(_smt.pSolver, sortbvtime, _olsqParam.min_depth);
        for (const BitwuzlaTerm* t: _smt.vTg){
            bitwuzla_assert(_smt.pSolver, bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_BV_ULT, t, depthBv));
        }
    }
    else{
        for(unsigned_t i = 0; i < _vpGateTimeWindow.size(); ++i){
        // for (const BitwuzlaTerm* t: _smt.vTg){
            const BitwuzlaTerm * depthBv = bitwuzla_mk_bv_value_uint64(_smt.pSolver, sortbvtime, _vpGateTimeWindow[i].second);
            bitwuzla_assert(_smt.pSolver, bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_BV_ULE,  _smt.vTg[i], depthBv));
        }
    }
}

void OLSQ::addSwapCountConstraintsZ3(unsigned_t bound){
    PB2CNF pb2cnf;
    vector<int_t> vSigmaLit;
    unsigned_t firstFreshVariable = 1 + (_olsqParam.min_depth+1)*_device.nEdge(), newFreashVariable;
    for (unsigned_t i = 1; i < firstFreshVariable; ++i){
        vSigmaLit.emplace_back(i);
    }
    vector< vector<int_t> > formula;
    newFreashVariable = pb2cnf.encodeAtMostK(vSigmaLit, bound, formula, firstFreshVariable) + 1;
    map<unsigned_t, const BitwuzlaTerm*> mAncillary;
    vector<const BitwuzlaTerm *> vOrs;
    unsigned_t sigmaT, sigmaE, var;
    const BitwuzlaSort *sortbool = bitwuzla_mk_bool_sort(_smt.pSolver);
    string s;
    for(vector<int_t>& clause : formula){
        vOrs.clear();
        for(int& lit : clause){
            var = abs(lit);
            if(var < firstFreshVariable){
                sigmaT = (var - 1) / _device.nEdge();
                sigmaE = (var - 1) % _device.nEdge();
                if (lit < 0){ 
                    vOrs.emplace_back(bitwuzla_mk_term1(_smt.pSolver, BITWUZLA_KIND_BV_NOT, _smt.vvSigma[sigmaT][sigmaE]));
                }
                else{
                    vOrs.emplace_back(_smt.vvSigma[sigmaT][sigmaE]);
                }
            }
            else{
                if(mAncillary.find(var) == mAncillary.end()){
                    s = "anc_" + to_string(var);
                    mAncillary[var] = bitwuzla_mk_const(_smt.pSolver, sortbool, s.c_str());
                }
                if (lit < 0){
                    vOrs.emplace_back(bitwuzla_mk_term1(_smt.pSolver, BITWUZLA_KIND_BV_NOT, mAncillary[var]));
                }
                else{
                    vOrs.emplace_back(mAncillary[var]);
                }
            }
        }    
        const BitwuzlaTerm * cnf = vOrs[0];
        for (unsigned_t i = 1; i < vOrs.size(); ++i){
            cnf = bitwuzla_mk_term2(_smt.pSolver, BITWUZLA_KIND_OR, vOrs[i], cnf);
        }
        bitwuzla_assert(_smt.pSolver, cnf);
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
        step = (_olsqParam.min_depth > 100) ? 10 : 1;
    }
    // FILE *ptr = fopen("constraint.txt","w");
    // bitwuzla_dump_formula(_smt.pSolver, "smt2", ptr);
    
    while(!find_min_depth && _olsqParam.min_depth < _olsqParam.max_depth ){
        fprintf(stdout, "[Info]          trying to optimize for depth bound %d            \n", _olsqParam.min_depth);
        _timer.start(TimeUsage::PARTIAL);
        bitwuzla_push(_smt.pSolver, 1);
        addDepthConstraintsZ3();
        success = checkModel();
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
            if(_olsqParam.min_depth + step < _olsqParam.max_depth){
                updateSMT(step);
                if(step > 1)
                    has_jump = true;
            }
            _olsqParam.min_depth += step;
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

    unsigned_t lower_swap_bound = 0;
    unsigned_t upper_swap_bound = (_olsqParam.is_use_SABRE_for_swap) ? _olsqParam.sabre_swap_bound : _pCircuit->nGate();
    upper_swap_bound = (_pCircuit->nSwapGate() < upper_swap_bound) ? _pCircuit->nSwapGate() : upper_swap_bound;
    bool reduce_swap = true;
    bool firstRun = true;
    unsigned_t step = (_olsqParam.is_transition) ? 3 : 9; 
    while (reduce_swap && _timer.fullRealTime() < _olsqParam.timeout){
        // cout << "enter loop" << endl;
        addDepthConstraintsZ3();
        reduce_swap = optimizeSwapForDepth(lower_swap_bound, upper_swap_bound, firstRun);
        upper_swap_bound = _pCircuit->nSwapGate() - 1;
        firstRun = false;
        // getchar();
        if(reduce_swap){
            fprintf(stdout, "[Info] Successfully reduce SWAP count. Go to next run.            \n");
            fprintf(stdout, "[Info] Solving with depth %d            \n", _olsqParam.min_depth + step);
            _timer.start(TimeUsage::PARTIAL);
            fprintf(stdout, "[Info] Generating formulation                        \n");
            if(_olsqParam.min_depth + step < _olsqParam.max_depth){
                updateSMT(step);
            }
            else{
                increaseDepthBound();
                generateFormulationZ3();
            }
            _olsqParam.min_depth += step;
            _timer.showUsage("Generating formulation", TimeUsage::PARTIAL);
            _timer.start(TimeUsage::PARTIAL);
        }
    }
    return true;
}

bool OLSQ::optimizeSwapForDepth(unsigned_t lower_swap_bound, unsigned_t upper_swap_bound, bool firstRun){
    unsigned_t swap_bound = upper_swap_bound;
    bool find_min_swap = false;
    bool success;
    while (!find_min_swap && lower_swap_bound <= swap_bound && swap_bound <= upper_swap_bound && _timer.fullRealTime() < _olsqParam.timeout){
        fprintf(stdout, "[Info]          trying to optimize for swap bound %d            \n", swap_bound);
        _timer.start(TimeUsage::PARTIAL);
        bitwuzla_push(_smt.pSolver, 1);
        addSwapCountConstraintsZ3(swap_bound);
        success = checkModel();
        fprintf(stdout, "[Info]          optimization results: %s                         \n", success ? "success" : "fail");
        _timer.showUsage("optimizing swap", TimeUsage::PARTIAL);
        _timer.showUsage("optimizing swap", TimeUsage::FULL);
        if (success){
            extractModel();
            if(swap_bound > _pCircuit->nSwapGate()){
                swap_bound = _pCircuit->nSwapGate();
            }
            --swap_bound;
        }
        else{
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
    unsigned_t circuitDepth = 0, i, gateTime, q, j, swapId, qId, t, e, tt;
    string s;
    // collect initial mapping
    vector<int_t> vInitialMapping(_pCircuit->nProgramQubit(), -1);

    // collect gate execution time
    vector<int_t> vQubitFirstGateTime(_pCircuit->nProgramQubit(), -1);
    vector<int_t> vQubitLastGateTime(_pCircuit->nProgramQubit(), -1);
    for ( i = 0; (i < _pCircuit->nGate());  ++i){
        Gate &gate = _pCircuit->gate(i);
        const char *rstr = bitwuzla_get_bv_value(_smt.pSolver, _smt.vTg[i]);
        s = rstr;
        gateTime = stoi(s, nullptr, 2); 
        circuitDepth = (circuitDepth < gateTime) ? gateTime : circuitDepth;
        gate.setExecutionTime(gateTime);
        for ( j = 0; j < gate.nTargetQubit();  ++j ){
            if(vQubitFirstGateTime[gate.targetProgramQubit(j)] == -1){
                vQubitFirstGateTime[gate.targetProgramQubit(j)] = gateTime;
            }
            vQubitLastGateTime[gate.targetProgramQubit(j)] = (vQubitLastGateTime[gate.targetProgramQubit(j)] < (int_t)gateTime) ? (int_t)gateTime: vQubitLastGateTime[gate.targetProgramQubit(j)];
        }
        for ( j = 0; j < gate.nTargetQubit();  ++j ){
            const char *rstr = bitwuzla_get_bv_value(_smt.pSolver, _smt.vvPi[gateTime][gate.targetProgramQubit(j)]);
            s = rstr;
            gate.addTargetPhysicalQubit(j, stoi(s, nullptr, 2));
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

    // collect qubit mapping
    vector<vector<int_t> > vvPhy2Pro(circuitDepth,  vector<int_t>(_device.nQubit(), -1));
    for (t = 0; t < circuitDepth; ++t){
        for (i = 0; i < _pCircuit->nProgramQubit(); ++i){
            const char *rstr = bitwuzla_get_bv_value(_smt.pSolver, _smt.vvPi[t][i]);
            s = rstr;
            vvPhy2Pro[t][stoi(s, nullptr, 2)] = i;
        }
    }

    // get SWAP gate
    _pCircuit->clearSwap();
    vector<unsigned_t> swapTargetQubit(2,0);
    bool cond1, cond2;
    int_t tmp, bound = circuitDepth - _olsqParam.swap_duration;
    vector<vector<bool> > vvTimeSwap(circuitDepth, vector<bool>(_device.nEdge(), 0));
    for (t = _olsqParam.swap_duration -1; t < circuitDepth; ++t){
        for (e = 0; e < _device.nEdge(); ++e){
            const char *rstr = bitwuzla_get_bv_value(_smt.pSolver, _smt.vvSigma[t][e]);
            s = rstr;
            vvTimeSwap[t][e] = stoi(s, nullptr, 2);
        }
    }
    for (e = 0; e < _device.nEdge(); ++e){
        for (t = _olsqParam.swap_duration -1; t < bound; ++t){
            if(vvTimeSwap[t][e] && vvTimeSwap[t + _olsqParam.swap_duration][e]){
                // only cancel two consecutive swap
                vvTimeSwap[t][e] = 0;
                vvTimeSwap[t + _olsqParam.swap_duration][e] = 0;
            }
        }
    }
    for (t = _olsqParam.swap_duration -1; t < circuitDepth; ++t){
        for (e = 0; e < _device.nEdge(); ++e){
            const char *rstr = bitwuzla_get_bv_value(_smt.pSolver, _smt.vvSigma[t][e]);
            s = rstr;
            if (vvTimeSwap[t][e]){
                cond1 = !((vvPhy2Pro[t][_device.edge(e).qubitId1()] == -1) && (vvPhy2Pro[t][_device.edge(e).qubitId2()] == -1));
                cond2 = cond1;
                // check if it is before the first gates of these qubits
                if(vvPhy2Pro[t][_device.edge(e).qubitId1()] != -1){
                    cond1 = (vQubitFirstGateTime[vvPhy2Pro[t][_device.edge(e).qubitId1()]] <= t);
                }
                if(vvPhy2Pro[t][_device.edge(e).qubitId2()] != -1){
                    cond1 = (vQubitFirstGateTime[vvPhy2Pro[t][_device.edge(e).qubitId2()]] <= t) || cond1;
                }
                // check if it is after the last gates of these qubits
                if(vvPhy2Pro[t][_device.edge(e).qubitId1()] != -1){
                    cond2 = (vQubitLastGateTime[vvPhy2Pro[t][_device.edge(e).qubitId1()]] > t);
                }
                if(vvPhy2Pro[t][_device.edge(e).qubitId2()] != -1){
                    cond2 = (vQubitLastGateTime[vvPhy2Pro[t][_device.edge(e).qubitId2()]] > t) || cond2;
                }
                // cerr << "cond1: " << cond1 << ", cond2: " << cond2 << endl;
                if(cond1 && cond2){
                    swapTargetQubit[0] = _device.edge(e).qubitId1();
                    swapTargetQubit[1] = _device.edge(e).qubitId2();
                    swapId = _pCircuit->nSwapGate();
                    _pCircuit->addSwapGate(swapId, swapTargetQubit, _olsqParam.swap_duration);
                    Gate & gate = _pCircuit->swapGate(swapId);
                    gate.setExecutionTime(t);
                    if (_verbose > 0){
                        fprintf(stdout, "        - SWAP Gate %d, duration: %d, time: %d, target qubit: %d %d\n", swapId + _pCircuit->nGate(), gate.duration(), gate.executionTime(), gate.targetPhysicalQubit(0), gate.targetPhysicalQubit(1));
                    }
                    if(vvPhy2Pro[t][swapTargetQubit[0]] != -1){
                        vQubitFirstGateTime[vvPhy2Pro[t][swapTargetQubit[0]]] = (vQubitFirstGateTime[vvPhy2Pro[t][swapTargetQubit[0]]] > t) ? t: vQubitFirstGateTime[vvPhy2Pro[t][swapTargetQubit[0]]];
                        vQubitLastGateTime[vvPhy2Pro[t][swapTargetQubit[0]]] = (vQubitLastGateTime[vvPhy2Pro[t][swapTargetQubit[0]]] < t) ? t + 1: vQubitLastGateTime[vvPhy2Pro[t][swapTargetQubit[0]]];
                    }
                    if(vvPhy2Pro[t][swapTargetQubit[0]] != -1){
                        vQubitFirstGateTime[vvPhy2Pro[t][swapTargetQubit[1]]] = (vQubitFirstGateTime[vvPhy2Pro[t][swapTargetQubit[1]]] > t) ? t: vQubitFirstGateTime[vvPhy2Pro[t][swapTargetQubit[1]]];
                        vQubitLastGateTime[vvPhy2Pro[t][swapTargetQubit[1]]] = (vQubitLastGateTime[vvPhy2Pro[t][swapTargetQubit[1]]] < t) ? t + 1: vQubitLastGateTime[vvPhy2Pro[t][swapTargetQubit[1]]];
                    }
                }
                // else{
                //     fprintf(stdout, "        - SWAP Gate, time: %d, target qubit: %d %d\n", t, _device.edge(e).qubitId1(), _device.edge(e).qubitId2());
                // }
            }
        }
    }
    // collect qubit mapping
    if (_verbose > 0){
        fprintf(stdout, "        - Qubit mapping: \n");
        for (t = 0; t <= circuitDepth; ++t){
        fprintf(stdout, "        - Time %d: ",t);
            for (i = 0; i < _pCircuit->nProgramQubit(); ++i){
                const char *rstr = bitwuzla_get_bv_value(_smt.pSolver, _smt.vvPi[t][i]);
                s = rstr;
                fprintf(stdout, "%d->%d ", i, stoi(s, nullptr, 2));
            }
        fprintf(stdout, "\n");
        }
    }
    
    for (i = 0; i < _pCircuit->nProgramQubit(); ++i){
        _pCircuit->setInitialMapping(i, vInitialMapping[i]);
        // cerr << "vQubitLastGateTime[" << i << "]:" << vQubitLastGateTime[i] << endl;
        // cerr << "vQubitFirstGateTime[" << i << "]:" << vQubitFirstGateTime[i] << endl;
        const char *rstr1 = bitwuzla_get_bv_value(_smt.pSolver, _smt.vvPi[vQubitFirstGateTime[i]][i]);
        s = rstr1;
        _pCircuit->setInitialMapping(i,  stoi(s, nullptr, 2));
        const char *rstr2 = bitwuzla_get_bv_value(_smt.pSolver, _smt.vvPi[vQubitLastGateTime[i]][i]);
        s = rstr2;
        _pCircuit->setFinalMapping(i,  stoi(s, nullptr, 2));
    }
    _pCircuit->setCircuitDepth(circuitDepth+1);
    outputQASM(_outputQASMFile);
}

void OLSQ::asapScheduling(){
    vector<int_t> vPushForwardDepth(_device.nQubit(), -1);
    int_t gateExecutionTime;
    unsigned_t block, i, j, q0, q1, qId, maxTime = 0;
    Gate gate;
    set<unsigned_t> sGateId;
    for (block = 0; block < _pCircuit->circuitDepth(); ++block){
        for (i = 0; i < _pCircuit->nGate(); ++i){
            Gate & gate = _pCircuit->gate(i);
            if (gate.executionTime() == block && sGateId.count(gate.idx()) == 0){
                gateExecutionTime = vPushForwardDepth[gate.targetPhysicalQubit(0)];
                for (j = 1; j < gate.nTargetQubit(); ++j){
                    qId = gate.targetPhysicalQubit(j);
                    gateExecutionTime = (gateExecutionTime < vPushForwardDepth[qId]) ? vPushForwardDepth[qId] : gateExecutionTime;
                }
                ++gateExecutionTime;
                for (j = 0; j < gate.nTargetQubit(); ++j){
                    qId = gate.targetPhysicalQubit(j);
                    vPushForwardDepth[qId] = gateExecutionTime;
                }
                maxTime = (maxTime < gateExecutionTime) ? gateExecutionTime : maxTime;
                gate.setExecutionTime(gateExecutionTime);
                sGateId.insert(gate.idx());
            }
        }
        if (block < _pCircuit->circuitDepth() - 1){
            for (j = 0; j < _pCircuit->nSwapGate(); ++j){
                Gate & gate = _pCircuit->swapGate(j);
                if (gate.executionTime() == block && sGateId.count(gate.idx()) == 0){
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



void OLSQ::increaseDepthBound(){
    while(_olsqParam.min_depth >= _olsqParam.max_depth){
        _olsqParam.max_depth_bit += _olsqParam.max_depth_expand_factor; 
        _olsqParam.max_depth = _olsqParam.max_depth << _olsqParam.max_depth_expand_factor;
    }
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

void OLSQ::updateSMT(unsigned_t d){
    fprintf(stdout, "[Info]          constructing injective mapping constraint    \n");
    addInjectiveMappingConstraintsZ3(d);
    fprintf(stdout, "[Info]          constructing valid two-qubit gate constraint \n");
    addValidTwoQubitGateConstraintsZ3(d); 
    fprintf(stdout, "[Info]          constructing swap overlapping constraint     \n");
    addSwapConstraintsZ3(d);
    fprintf(stdout, "[Info]          constructing mapping transformation constraint\n");
    addTransformationConstraintsZ3(d);
    if(_olsqParam.use_window_range_for_gate){
        updateGateTimeWindow(d);
    }
}

void OLSQ::constructGateTimeWindow(){
    vector<int_t> vQubitArriveTime(_pCircuit->nProgramQubit(), 0);
    vector<int_t> vQubitRequiredTime(_pCircuit->nProgramQubit(), _olsqParam.min_depth-1);
    unsigned_t qId, j, t;
    int_t i;
    for ( i = 0; i < _pCircuit->nGate(); ++i){
        Gate & gate = _pCircuit->gate(i);
        t = 0;
        for ( j = 0; j < gate.nTargetQubit();  ++j ){
            qId = gate.targetProgramQubit(j);
            t = (vQubitArriveTime[qId] < t) ? t : vQubitArriveTime[qId];
        }
        _vpGateTimeWindow.emplace_back(make_pair(t, 0));
        ++t;
        for ( j = 0; j < gate.nTargetQubit();  ++j ){
            vQubitArriveTime[gate.targetProgramQubit(j)] = t;
        }
        // cerr << "qubit arrive time: " << endl;
        // for (auto qat : vQubitArriveTime){
        //     cerr << qat << " ";
        // }
        // cerr << endl;
    }
    // cerr << "_olsqParam.min_depth: " << _olsqParam.min_depth << endl;
    for ( i = _pCircuit->nGate() - 1; i >= 0; --i){
        Gate & gate = _pCircuit->gate(i);
        t = _olsqParam.min_depth - 1;
        // cerr << "before t: " << t << endl;
        for ( j = 0; j < gate.nTargetQubit();  ++j ){
            qId = gate.targetProgramQubit(j);
            t = (vQubitRequiredTime[qId] > t) ? t : vQubitRequiredTime[qId];
            // cerr << "update t: " << t << endl;
            
        }
        // cerr << "t: " << t << endl;
        _vpGateTimeWindow[i].second = t;
        --t;
        for ( j = 0; j < gate.nTargetQubit();  ++j ){
            vQubitRequiredTime[gate.targetProgramQubit(j)] = t;
        }
        
    }
    printGateTimeWindow();
}

void OLSQ::updateGateTimeWindow(unsigned_t d){
    for(int_t i = 0; i < _vpGateTimeWindow.size(); ++i){
        _vpGateTimeWindow[i].second += d;
    }
    // printGateTimeWindow();
}

void OLSQ::printGateTimeWindow(){
    fprintf(stdout, "[Info] Gate Time Window Info                              \n");
    fprintf(stdout, "       ------------------------------------------\n");
    fprintf(stdout, "       - Min depth: %d                           \n", _olsqParam.min_depth);
    fprintf(stdout, "       - Gate time window:\n");
    for (unsigned_t i = 0; i < _vpGateTimeWindow.size();  ++i ){
        fprintf(stdout, "        - Gate %d: [%d, %d]\n", i, _vpGateTimeWindow[i].first, _vpGateTimeWindow[i].second);
    }
}


void OLSQ::collectQubitRegion(){
    for(unsigned_t i = 0; i < _pCircuit->nProgramQubit(); ++i){
        // collect circuit qubit in the current device
        bfsSearch(i);
    }
}

void OLSQ::bfsSearch(unsigned_t q){
    // collect qubit region from mapping
    // cerr << "start bfs search for pro q " << qId << endl;
    // for(unsigned_t k = 0; k < _pDevice->nQubit(); ++k){
    //     sQubitRegion.insert(k);
    // }
    // for (unsigned_t e = 0; e < _pDevice->nEdge(); ++e){
    //     sSwapRegion.insert(e);
    // }
    // return;
    map<int_t, int_t> mQubitPathLength;
    set<int_t>& sQubitRegion = _pCircuit->qubitRegion(q);
    set<int_t>::iterator it = sQubitRegion.begin();
    ++it;
    for (; it != sQubitRegion.end(); ++it) {
        mQubitPathLength[*it] = -1;

    }
    // cerr << "init qubit region: ";
    // for (int_t s : sQubitRegion){
    //     cerr << s << " ";
    // }
    // cerr << endl;
    set<int_t> sVisitedQubitRegion;
    unsigned_t i, j, p, nIdx;
    vector<Node> vNode;
    priority_queue<int_t, std::vector<int>, std::greater<int> > priorityQ;
    vector<int_t> vBacktraceNode;
    unsigned_t cost, nSpanEdge, maxDis = 0;
    vNode.emplace_back(-1, 0, (*sQubitRegion.begin()), -1, 0);
    priorityQ.push(0);
    while(priorityQ.size() > 0 && (sVisitedQubitRegion.size() != sQubitRegion.size() || vNode[priorityQ.top()].dis <= maxDis)){
        nIdx = priorityQ.top();
        priorityQ.pop();
        // cerr << "expand node " << vNode[nIdx].idx << " phy q: " << vNode[nIdx].qIdx << " parent idx: " << vNode[nIdx].parentIdx << ", dis: " << vNode[nIdx].dis << endl;
        Qubit& qubit = _device.qubit(vNode[nIdx].qIdx);
        if(sQubitRegion.find(vNode[nIdx].qIdx) != sQubitRegion.end() && (mQubitPathLength[vNode[nIdx].qIdx] == -1 || mQubitPathLength[vNode[nIdx].qIdx] >= vNode[nIdx].dis)){
            // cerr << "find node " << vNode[nIdx].idx << " qubit " << vNode[nIdx].qIdx << " with distance " << vNode[nIdx].dis << " parent idx " << vNode[nIdx].parentIdx << endl;
            maxDis = (vNode[nIdx].dis > maxDis) ? vNode[nIdx].dis : maxDis;
            sVisitedQubitRegion.insert(vNode[nIdx].qIdx);
            vBacktraceNode.emplace_back(vNode[nIdx].idx);
        }
        // sTraversedQubit.insert(node.qIdx);
        nSpanEdge = qubit.vSpanEdge.size();
        for (j = 0; j < nSpanEdge; ++j){
            p = _device.edge(qubit.vSpanEdge[j]).qubitId1();
            if (qubit.idx == p){
                p = _device.edge(qubit.vSpanEdge[j]).qubitId2();
            }
                // if (sTraversedQubit.find(q) == sTraversedQubit.end()){
                // sTraversedQubit.insert(q);
            // cerr << "node idx " << vNode[nIdx].idx << endl;
            vNode.emplace_back(vNode[nIdx].idx, vNode.size(), p, qubit.vSpanEdge[j], vNode[nIdx].dis+1);
            // cerr << "add node " << vNode.back().idx << " qubit " << vNode.back().qIdx << " with distance " << vNode.back().dis << " parent idx " << vNode.back().parentIdx << endl;
            // queue.emplace_back(vNode.back().idx);
            priorityQ.push(vNode.back().idx);
        }
        // cerr << "top node dis: " << vNode[priorityQ.top()].dis << endl;
        // getchar();
    }
    // backtrack 
    // cerr << "start backtrace" << endl;
    int_t pIdx;
    for(int_t nIdx : vBacktraceNode){
        pIdx = nIdx;
        // cerr << "find node " << vNode[nIdx].idx << " qubit " << vNode[nIdx].qIdx << " with distance " << vNode[nIdx].dis << " parent idx " << vNode[nIdx].parentIdx << endl;
        while(pIdx > 0){
            // cerr << "pIdx: " << pIdx << endl;
            sQubitRegion.insert(vNode[pIdx].qIdx);
            // cerr << "add edge " << vNode[pIdx].parentEIdx << endl;
            pIdx = vNode[pIdx].parentIdx;
        }
    }

    expandRegion(q);
    
    // if(sQubitRegion.size() == 1){
    //     q = *(sQubitRegion.begin());
    //     Qubit& qubit = _pDevice->qubit(q);
    //     nSpanEdge = qubit.vSpanEdge.size();
    //     for(i = 0; i < nSpanEdge; ++i){
    //         sQubitRegion.insert(i);
    //         sSwapRegion.insert(qubit.vSpanEdge[i]);
    //     }
    // }
    // cerr << "final qubit region: ";
    // for (int_t s : sQubitRegion){
    //     cerr << s << " ";
    // }
    // cerr << endl;
}


void OLSQ::expandRegion(unsigned_t q){
    unsigned_t nSpanEdge;
    int_t p;
    set<int_t> sExpandQubit;
    set<int_t>& sQubitRegion = _pCircuit->qubitRegion(q);
    for(unsigned_t i = 0; i < sQubitRegion.size(); ++i){
        set<int_t>::iterator it;
        sExpandQubit.clear();
        for (it = sQubitRegion.begin(); it != sQubitRegion.end(); ++it) {
            p = *(it);
            Qubit& qubit = _device.qubit(p);
            nSpanEdge = qubit.vSpanEdge.size();
            for(i = 0; i < nSpanEdge; ++i){
                sExpandQubit.insert(i);
            }
        }
        _pCircuit->addQubitRegion(q, sExpandQubit);
    }
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
        vvTimeGate[gate.executionTime()].emplace_back(i+_pCircuit->nGate());
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
    line = line + "\n// measurement\n";
    for (i = 0; i < _pCircuit->nProgramQubit(); ++i){
        line = line + "measure q[" + to_string(_pCircuit->finalMapping(i)) + "]->c[" + to_string(i) + "];\n";
    }
    FILE* fout = fopen(fileName.c_str(), "w");
    fprintf(fout, "%s", line.c_str());
    fclose(fout);
}


string OLSQ::outputQASMStr(){
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
        vvTimeGate[gate.executionTime()].emplace_back(i+_pCircuit->nGate());
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
    line = line + "\n// measurement\n";
    for (i = 0; i < _pCircuit->nProgramQubit(); ++i){
        line = line + "measure q[" + to_string(_pCircuit->finalMapping(i)) + "]->c[" + to_string(i) + "];\n";
    }
    return line.c_str();
}

OLSQ_NAMESPACE_CPP_END