/***********************************************************************
  File        [ device.hpp ]
  System      [ mOLSQ: multilevel quantum layout synthesis tool]
  Package     [ device ]
  Synopsis    [ device class header ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 22, Nov., 2022 ]
***********************************************************************/
#ifndef DEVICE_HPP
#define DEVICE_HPP

#include "misc/global.hpp"
#include <unordered_map>
#include <unordered_set>
#include <deque>

MOLSQ_NAMESPACE_HPP_START


////////////////////////////
// Struct Qubit
////////////////////////////
struct Qubit {
    Qubit():
        idx(0), vSpanEdge(0), x(0), y(0)
        {}
    Qubit(unsigned_t i):
        idx(i), vSpanEdge(0), x(0), y(0)
        {}
    Qubit(unsigned_t i, unsigned_t x, unsigned_t y):
        idx(i), vSpanEdge(0), x(x), y(y)
        {}
    ~Qubit() {}
    unsigned_t            idx; 
    vector<unsigned_t>    vSpanEdge;
    unsigned_t            x; // x coordinate 
    unsigned_t            y; // y coordinate 
};

////////////////////////////
// Struct Edge
////////////////////////////
struct Edge {
    Edge():
        idx(0), prEndPoints(0,0)
        {}
    Edge(unsigned_t i, unsigned_t p0, unsigned_t p1):
        idx(i), prEndPoints(p0,p1)
        {}
    ~Edge() {}
    unsigned_t  qubitId1()   const { return prEndPoints.first; }
    unsigned_t  qubitId2()   const { return prEndPoints.second; }
    unsigned_t            idx;
    pair<unsigned_t, unsigned_t>   prEndPoints;
};

class Device{
    public:
        Device(const string& name) :   
            _name(name)
            {
                _umXY2QubitId.clear();
                constructQubit();
                _vEdge.clear();
            }
        Device(const string& name, unsigned_t nQubit, unsigned_t nEdge) :   
            _name(name), _nQubit(nQubit), _nEdge(nEdge)
            {
                _umXY2QubitId.clear();
                constructQubit();
                _vEdge.clear();
                _vEdge.reserve(nEdge);
                _vvAdjacentMatrix.resize(nQubit, vector<bool>(nQubit, false));
            }
        Device(const string& name, unsigned_t nQubit, unsigned_t nEdge, vector<pair<unsigned_t, unsigned_t> > vQubitXY) :   
            _name(name), _nQubit(nQubit), _nEdge(nEdge)
            {
                _umXY2QubitId.clear();
                constructQubit(vQubitXY);
                _vEdge.clear();
                _vEdge.reserve(nEdge);
                _vvAdjacentMatrix.resize(nQubit, vector<bool>(nQubit, false));
            }
        ~Device() {}
        
        string             name()                               const { return _name; }
        unsigned_t         nQubit()                             const { return _nQubit; }
        unsigned_t         nEdge()                              const { return _nEdge; }
        unsigned_t         maxX()                               const { return _maxX; }
        unsigned_t         maxY()                               const { return _maxY; }
        
        bool               isValidQubitIdx(unsigned_t idx)      const { return 0 <= idx && idx < _nQubit;}
        bool               isValidEdgeIdx(unsigned_t idx)       const { return 0 <= idx && idx < _nEdge;}


        Qubit&             qubit(unsigned_t idx)  { assert(isValidQubitIdx(idx)); return _vQubit[idx]; }
        Edge&              edge(unsigned_t idx)   { assert(isValidEdgeIdx(idx)); return _vEdge[idx]; }
        unsigned_t         qubitId(unsigned_t x, unsigned_t y)  { assert(x <= _maxX && y <= _maxY); return _umXY2QubitId[make_pair(x,y)]; }

        void setQubit(unsigned_t nQubit) { _nQubit = nQubit; constructQubit(); _vvAdjacentMatrix.resize(nQubit, vector<bool>(nQubit, false)); };
        void setXYCoord(vector<pair<unsigned_t, unsigned_t> >& vXY){
            unsigned_t i;
            for (i = 0; i < _nQubit; ++i){
                _vQubit[i].x = vXY[i].first;
                _vQubit[i].y = vXY[i].second;
                _umXY2QubitId[vXY[i]] = i;
            }
        }
        void setEdge(vector<pair<unsigned_t, unsigned_t> >& vEdge);
        void addEdge(unsigned_t q0, unsigned_t q1);
        void printDevice();

        void calAPSP();
        unsigned_t getDistance(unsigned_t i, unsigned_t j){
            return _vvASAP[i][j];
        }
        bool isAdjacent(unsigned_t i, unsigned_t j) { return _vvAdjacentMatrix[i][j]; }

    private:
        ////////////////////////////
        // Private member
        ////////////////////////////
        string          _name;
        unsigned_t      _nQubit;
        unsigned_t      _nEdge;
        unsigned_t      _maxX;
        unsigned_t      _maxY;
        vector<Qubit>   _vQubit;
        vector<Edge>    _vEdge;
        unordered_map<pair<unsigned_t, unsigned_t>, unsigned_t, pair_hash> _umXY2QubitId;

        vector<vector<unsigned_t> > _vvASAP;
        vector<vector<bool> >       _vvAdjacentMatrix;

        ////////////////////////////
        // Private functions
        ////////////////////////////
        void constructQubit();
        void constructQubit(vector<pair<unsigned_t, unsigned_t> > vQubitXY);
};


MOLSQ_NAMESPACE_HPP_END

#endif //DEVICE_HPP