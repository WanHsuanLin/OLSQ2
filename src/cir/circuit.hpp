/***********************************************************************
  File        [ circuit.hpp ]
  System      [ OLSQ: optimal quantum layout synthesis tool]
  Package     [ cir ]
  Synopsis    [ circuit class header ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 22, Nov., 2022 ]
***********************************************************************/
#ifndef CIRCUIT_HPP
#define CIRCUIT_HPP

#include "misc/global.hpp"
#include "cir/gate.hpp"


OLSQ_NAMESPACE_HPP_START

class Circuit
{
    public:
        Circuit() :
            _name(""), _nProgramQubit(0)
            {
                _vGate.clear();
                _vSwapGate.clear();
                _vInitialMapping.clear();
                _vFinalMapping.clear();
            }
        Circuit(const string& cirName, unsigned_t nProgramQubit, unsigned_t nGate):
            _name(cirName), _nProgramQubit(nProgramQubit)
            {
                _vGate.clear();
                _vGate.reserve(nGate);
                _vSwapGate.clear();
                _vInitialMapping.clear();
                _vInitialMapping.resize(nProgramQubit, -1);
                _vFinalMapping.clear();
                _vFinalMapping.resize(nProgramQubit, -1);
            }
        ~Circuit() {}

        
        // get function
        string             name()                                 const { return _name; }
        unsigned_t         nProgramQubit()                        const { return _nProgramQubit; }
        unsigned_t         nGate()                                const { return _vGate.size(); }
        unsigned_t         nSwapGate()                            const { return _vSwapGate.size(); }
        unsigned_t         circuitDepth()                         const { return _circuitDepth; }
        bool               isValidQubitIdx(unsigned_t idx)        const { return 0 <= idx && idx < _nProgramQubit;}
        bool               isValidGateIdx(unsigned_t idx)         const { return 0 <= idx && idx < _vGate.size();}
        bool               isValidSwapGateIdx(unsigned_t idx)     const { return 0 <= idx && idx < _vSwapGate.size();}
        Gate&              gate(unsigned_t idx)                         { assert(isValidGateIdx(idx)); return _vGate[idx]; }
        Gate&              swapGate(unsigned_t idx)                     { assert(isValidSwapGateIdx(idx)); return _vSwapGate[idx]; }
        int_t              initialMapping(unsigned_t idx)         const { assert(isValidQubitIdx(idx)); return _vInitialMapping[idx]; }
        int_t              finalMapping(unsigned_t idx)           const { assert(isValidQubitIdx(idx)); return _vFinalMapping[idx]; }
        // set function
        void addGate( string const & gateName, vector<unsigned_t> const & vTargetQubit, unsigned_t duration = 1);
        // void addGate( string const & gateName, vector<int_t> const & vTargetQubit, unsigned_t duration = 1);
        void addSwapGate(vector<unsigned_t> const & vTargetQubit, unsigned_t duration = 1);
        void setInitialMapping(int_t pro_q, int_t phy_q)        { _vInitialMapping[pro_q] = phy_q; }
        void setInitialMapping(vector<unsigned_t> const &  vInitialMapping);
        void setFinalMapping(int_t pro_q, int_t phy_q)          { _vFinalMapping[pro_q] = phy_q; }
        void setCircuitDepth(int_t t)                           { _circuitDepth = t; }
        void clearSwap()                                        { _vSwapGate.clear(); }

        void printCircuit();
        void printCircuitLayout();

    private:
        string                       _name;
        unsigned_t                   _nProgramQubit;
        vector<Gate>                 _vGate;

        // store the layout synthesis result
        unsigned_t                   _circuitDepth;
        vector<Gate>                 _vSwapGate;
        vector<int_t>                _vInitialMapping;
        vector<unsigned_t>           _vFinalMapping;
};

OLSQ_NAMESPACE_HPP_END

#endif // CIRCUIT_HPP
