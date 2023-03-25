/***********************************************************************
  File        [ cir.cpp ]
  System      [ OLSQ: optimal quantum layout synthesis tool]
  Package     [ cir ]
  Synopsis    [ Circuit class implementation ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 22, Nov., 2022 ]
***********************************************************************/
#include "cir/circuit.hpp"

OLSQ_NAMESPACE_CPP_START

void Circuit::addGate( string const & gateName, vector<unsigned_t> const & vTargetQubit, unsigned_t duration){
    unsigned_t idx = _vGate.size();
    _vGate.emplace_back(idx, gateName, duration);
    for(unsigned_t q: vTargetQubit){
        assert(isValidQubitIdx(q));
        _vGate[idx].addTargetProgramQubit(q);
    }
}

void Circuit::addSwapGate(vector<unsigned_t> const & vTargetQubit, unsigned_t duration){
    unsigned_t idx = _vSwapGate.size(), i;
    _vSwapGate.emplace_back(_vGate.size() + _vSwapGate.size(), "swap", duration);
    for(i = 0; i < vTargetQubit.size(); ++i){
        // assert(isValidQubitIdx(vTargetQubit[i]));
        _vSwapGate[idx].addTargetProgramQubit(vTargetQubit[i]);
        _vSwapGate[idx].addTargetPhysicalQubit(i, vTargetQubit[i]);
    }
}

void Circuit::setInitialMapping(vector<unsigned_t> const &  vInitialMapping){
    assert(vInitialMapping.size() == _nProgramQubit);
    for(unsigned_t i = 0; i < _nProgramQubit; i++){
        assert(isValidQubitIdx(vInitialMapping[i]));
        _vInitialMapping[i] = vInitialMapping[i];
    }
}

void Circuit::printCircuit(){
    unsigned_t gateId = 0, i;
    fprintf(stdout, "[Info] Circuit Info                              \n");
    fprintf(stdout, "       ------------------------------------------\n");
    fprintf(stdout, "       - name             : %s\n", _name.c_str());
    fprintf(stdout, "       - #Program Qubit   : %d\n", _nProgramQubit);
    fprintf(stdout, "       ------------------------------------------\n");
    fprintf(stdout, "       - Gate list\n");
    for (Gate &gate: _vGate){
        fprintf(stdout, "        - Gate %d, name: %s, duration: %d, target qubit:", gateId, gate.name().c_str(), gate.duration());
        for ( i = 0; i < gate.nTargetQubit();  ++i )
            fprintf(stdout, " %d", gate.targetProgramQubit(i));
        fprintf(stdout, "\n");
        ++gateId;
    }
    fprintf(stdout, "\n");
}


void Circuit::printCircuitLayout(){
    unsigned_t gateId, qId, i;
    Gate gate;
    fprintf(stdout, "[Info] Compiled Circuit Info                              \n");
    fprintf(stdout, "       ------------------------------------------\n");
    fprintf(stdout, "       - name             : %s\n", _name.c_str());
    fprintf(stdout, "       - #Program Qubit   : %d\n", _nProgramQubit);
    fprintf(stdout, "       - #swap            : %lu\n", _vSwapGate.size());
    fprintf(stdout, "       ------------------------------------------\n");
    fprintf(stdout, "       - Initial mapping\n");
    for (qId = 0; qId < _nProgramQubit; ++qId){
        fprintf(stdout, "        - Program qubit %d is mapped to physical qubit: %d\n", qId, _vInitialMapping[qId]);
    }
    fprintf(stdout, "       ------------------------------------------\n");
    fprintf(stdout, "       - Final mapping\n");
    for (qId = 0; qId < _nProgramQubit; ++qId){
        fprintf(stdout, "        - Program qubit %d is mapped to physical qubit: %d\n", qId, _vFinalMapping[qId]);
    }
    fprintf(stdout, "       ------------------------------------------\n");
    fprintf(stdout, "       - SWAP Gate list\n");
    gateId = 0;
    for (Gate &gate: _vSwapGate){
        fprintf(stdout, "        - SWAP Gate %d, duration: %d, time: %d, target qubit: %d %d\n", gateId, gate.duration(), gate.executionTime(), gate.targetPhysicalQubit(0), gate.targetPhysicalQubit(1));
        ++gateId;
    }
    fprintf(stdout, "       ------------------------------------------\n");
    fprintf(stdout, "       - Gate list\n");
    gateId = 0;
    for (Gate &gate: _vGate){
        fprintf(stdout, "        - Gate %d, name: %s, duration: %d, time: %d, target program qubit:", gateId, gate.name().c_str(), gate.duration(), gate.executionTime());
        for ( i = 0; i < gate.nTargetQubit();  ++i ){
            fprintf(stdout, " %d", gate.targetProgramQubit(i));
        }
        fprintf(stdout, ", target physical qubit: ");
        for ( i = 0; i < gate.nTargetQubit();  ++i ){
            fprintf(stdout, "%d ", gate.targetPhysicalQubit(i));
        }
        fprintf(stdout, "\n");
        ++gateId;
    }
    fprintf(stdout, "\n");
}

OLSQ_NAMESPACE_CPP_END