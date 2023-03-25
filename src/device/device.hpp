/***********************************************************************
  File        [ device.hpp ]
  System      [ OLSQ: optimal quantum layout synthesis tool]
  Package     [ device ]
  Synopsis    [ device class header ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 22, Nov., 2022 ]
***********************************************************************/
#ifndef DEVICE_HPP
#define DEVICE_HPP

#include "misc/global.hpp"

OLSQ_NAMESPACE_HPP_START

////////////////////////////
// Struct Qubit
////////////////////////////
struct Qubit {
    Qubit():
        idx(0), vSpanEdge(0)
        {}
    Qubit(unsigned_t i):
        idx(i), vSpanEdge(0)
        {}
    ~Qubit() {}
    unsigned_t            idx; // f = g + h
    vector<unsigned_t>    vSpanEdge;
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
        Device() {}
        Device(const string& name, unsigned_t nQubit, unsigned_t nEdge) :   
            _name(name), _nQubit(nQubit), _nEdge(nEdge)
            {
                constructQubit();
                _vEdge.clear();
                _vEdge.reserve(nEdge);
            }
        ~Device() {}
        
        string             name()                               const { return _name; }
        unsigned_t         nQubit()                             const { return _nQubit; }
        unsigned_t         nEdge()                              const { return _nEdge; }
        
        bool               isValidQubitIdx(unsigned_t idx)      const { return 0 <= idx && idx < _nQubit;}
        bool               isValidEdgeIdx(unsigned_t idx)       const { return 0 <= idx && idx < _nEdge;}


        Qubit&              qubit(unsigned_t idx)  { assert(isValidQubitIdx(idx)); return _vQubit[idx]; }
        Edge&               edge(unsigned_t idx)   { assert(isValidEdgeIdx(idx)); return _vEdge[idx]; }

        void setEdge(vector<pair<unsigned_t, unsigned_t> >& vEdge);
        void addEdge(unsigned_t q0, unsigned_t q1);
        void printDevice();

    private:
        ////////////////////////////
        // Private member
        ////////////////////////////
        string          _name;
        unsigned_t      _nQubit;
        unsigned_t      _nEdge;
        vector<Qubit>   _vQubit;
        vector<Edge>    _vEdge;

        ////////////////////////////
        // Private functions
        ////////////////////////////
        void constructQubit();
};


OLSQ_NAMESPACE_HPP_END

#endif //DEVICE_HPP