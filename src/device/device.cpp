/***********************************************************************
  File        [ device.cpp ]
  System      [ OLSQ: optimal quantum layout synthesis tool]
  Package     [ device ]
  Synopsis    [ Device class implementation ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 22, Nov., 2022 ]
***********************************************************************/
#include "device/device.hpp"

OLSQ_NAMESPACE_CPP_START

void Device::constructQubit(){
    _vQubit.clear();
    unsigned_t i;
    _vQubit.reserve(_nQubit);
    for (i = 0; i < _nQubit; ++i){
        _vQubit.emplace_back(i);
    }
}

void Device::setEdge(vector<pair<unsigned_t, unsigned_t> >& vEdge){
    for (unsigned_t i = 0; i < vEdge.size(); ++i){
        assert(isValidQubitIdx(vEdge[i].first));
        assert(isValidQubitIdx(vEdge[i].second));
        _vEdge.emplace_back(i, vEdge[i].first, vEdge[i].second);
        _vQubit[vEdge[i].first].vSpanEdge.emplace_back(i);
        _vQubit[vEdge[i].second].vSpanEdge.emplace_back(i);
    }
}

void Device::addEdge(unsigned_t q0, unsigned_t q1){
    unsigned_t idx = _vEdge.size();
    assert(isValidQubitIdx(q0));
    assert(isValidQubitIdx(q1));
    _vEdge.emplace_back(idx, q0, q1);
    _vQubit[q0].vSpanEdge.emplace_back(idx);
    _vQubit[q1].vSpanEdge.emplace_back(idx);
}

void Device::printDevice(){
    unsigned_t edgeId, i, nSpanEdge, j;
    Edge edge;
    fprintf(stdout, "[Info] Device Info                              \n");
    fprintf(stdout, "       ------------------------------------------\n");
    fprintf(stdout, "       - name             : %s\n", _name.c_str());
    fprintf(stdout, "       - #Physical Qubit  : %d\n", _nQubit);
    fprintf(stdout, "       - #Edge            : %d\n", _nEdge);
    fprintf(stdout, "       ------------------------------------------\n");
    fprintf(stdout, "       - Edge list\n");
    for ( i = 0; i < _nEdge;  ++i ){
        fprintf(stdout, "        - Edge %d (%d,%d)\n", i, _vEdge[i].qubitId1(), _vEdge[i].qubitId2());
    }
    fprintf(stdout, "\n");
    fprintf(stdout, "       ------------------------------------------\n");
    fprintf(stdout, "       - Qubit span edge list\n");
    for (i = 0; i < _nQubit; ++i){
        Qubit& qubit = _vQubit[i];
        fprintf(stdout, "        - Span edge for qubit %d: ", i);
        nSpanEdge = qubit.vSpanEdge.size();
        for (j = 0; j < nSpanEdge; ++j){
            fprintf(stdout, " (%d,%d)", _vEdge[qubit.vSpanEdge[j]].qubitId1(), _vEdge[qubit.vSpanEdge[j]].qubitId2());
        }
        fprintf(stdout, "\n");
    }
}



OLSQ_NAMESPACE_CPP_END