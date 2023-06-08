/***********************************************************************
  File        [ device.cpp ]
  System      [ mOLSQ: multilevel quantum layout synthesis tool]
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

void Device::constructQubit(vector<pair<unsigned_t, unsigned_t> > vQubitXY){
    _vQubit.clear();
    unsigned_t i;
    _vQubit.reserve(_nQubit);
    for (i = 0; i < _nQubit; ++i){
        _vQubit.emplace_back(i, vQubitXY[i].first, vQubitXY[i].second);
        _maxX = max(_maxX, vQubitXY[i].first);
        _maxY = max(_maxY, vQubitXY[i].second);
        _umXY2QubitId[vQubitXY[i]] = i;
    }
}

void Device::setEdge(vector<pair<unsigned_t, unsigned_t> >& vEdge){
    for (unsigned_t i = 0; i < vEdge.size(); ++i){
        addEdge(vEdge[i].first, vEdge[i].second);
    }
}

void Device::addEdge(unsigned_t q0, unsigned_t q1){
    unsigned_t idx = _vEdge.size();
    assert(isValidQubitIdx(q0));
    assert(isValidQubitIdx(q1));
    _vEdge.emplace_back(idx, q0, q1);
    _vQubit[q0].vSpanEdge.emplace_back(idx);
    _vQubit[q1].vSpanEdge.emplace_back(idx);
    _vvAdjacentMatrix[q0][q1] = true;
    _vvAdjacentMatrix[q1][q0] = true;
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

void Device::calAPSP(){
    _vvASAP.clear();
    _vvASAP.resize(_nQubit, vector<unsigned_t>(_nQubit, MAX_UNSIGNED));
    // unsigned_t i, j, k;/
    // for(i = 0; i < _nEdge; ++i){
    //     _vvASAP[_vEdge[i].qubitId1()][_vEdge[i].qubitId2()] = 1;
    // }
    // for(i = 0; i < _nQubit; ++i){
    //     _vvASAP[i][i]= 0;
    // }
    // for(k = 0; k < _nQubit; ++k){ 
    //     for(i = 0; i < _nQubit; ++i){ 
    //         for(j = 0; j < _nQubit; ++j){ 
    //             if (_vvASAP[i][j] > _vvASAP[i][k] + _vvASAP[k][j]){
    //                 _vvASAP[i][j] = _vvASAP[i][k] + _vvASAP[k][j];
    //             }
    //         }
    //     }
    // }
    unsigned_t i, j, q;
    deque<unsigned_t> queue;
    unsigned_t cost, nSpanEdge;
    unordered_set<unsigned_t> sTraversedQubit;
    for(i = 0; i < _nQubit; ++i){ 
        // run bfs to all other nodes
        sTraversedQubit.clear();
        queue.clear();
        queue.emplace_back(i);
        _vvASAP[i][i] = 0;
        while(queue.size() > 0){
            Qubit& qubit = _vQubit[queue[0]];
            sTraversedQubit.insert(queue[0]);
            queue.pop_front();
            nSpanEdge = qubit.vSpanEdge.size();
            for (j = 0; j < nSpanEdge; ++j){
                q = _vEdge[qubit.vSpanEdge[j]].qubitId1();
                if (qubit.idx != q){
                    if (sTraversedQubit.find(q) == sTraversedQubit.end()){
                        sTraversedQubit.insert(q);
                        _vvASAP[i][q] = _vvASAP[i][qubit.idx] + 1; 
                        queue.emplace_back(q);
                        // cout << "phy " << i << " to phy " << q << ": dis " << cost << endl;
                    }

                }
                else{
                    q = _vEdge[qubit.vSpanEdge[j]].qubitId2();
                    if (sTraversedQubit.find(q) == sTraversedQubit.end()){
                        sTraversedQubit.insert(q);
                        queue.emplace_back(q);
                        _vvASAP[i][q] = _vvASAP[i][qubit.idx] + 1;
                        // cout << "phy " << i << " to phy " << q << ": dis " << cost << endl;
                    }
                }
            }
        }
    }

    // print distance
    // fprintf(stdout, "[Info] Device Qubit Distance Info                              \n");
    // for(i = 0; i < _nQubit; ++i){ 
    //     for(j = i; j < _nQubit; ++j){ 
    //         assert(_vvASAP[i][j] == _vvASAP[j][i]);
    //         fprintf(stdout, "        - (%d,%d): dis %d\n", i, j, _vvASAP[i][j]);
    //     }
    // }
}



OLSQ_NAMESPACE_CPP_END