/***********************************************************************
  File        [ circuit.hpp ]
  System      [ mOLSQ: multilevel quantum layout synthesis tool]
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


MOLSQ_NAMESPACE_HPP_START

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
                _vpGateDependency.clear();
            }
        Circuit(const string& cirName) :
            _name(cirName), _nProgramQubit(0)
            {
                _vGate.clear();
                _vSwapGate.clear();
                _vInitialMapping.clear();
                _vFinalMapping.clear();
                _vpGateDependency.clear();
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
                _vpGateDependency.clear();
                _vsQubitRegion.resize(nProgramQubit);
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
        void               qasm(string const & fileName);
        string             qasmstr();
        // set function
        void setQubitNum(unsigned_t n) { _nProgramQubit = n; _vsQubitRegion.clear(); _vsQubitRegion.resize(_nProgramQubit); }
        void constructDependency();
        void addDependency(Gate& g1, Gate& g2){
            _vpGateDependency.emplace_back(make_pair(g1.idx(), g2.idx()));
        }
        void addDependency(unsigned_t g1, unsigned_t g2){
            _vpGateDependency.emplace_back(make_pair(g1, g2));
        }
        void setDependency(vector<pair<unsigned_t, unsigned_t> >& vDependencies){
            for (pair<unsigned_t, unsigned_t>& p : vDependencies){
                assert(isValidGateIdx(p.first));
                assert(isValidGateIdx(p.second));
                addDependency(p.first, p.second);
            }
        }
        void addGate( Gate & gate ) { _vGate.emplace_back(gate); };
        void addGate( string const & gateName, vector<unsigned_t> const & vTargetQubit, unsigned_t duration = 1);
        // void addGate( string const & gateName, vector<int_t> const & vTargetQubit, unsigned_t duration = 1);
        void addSwapGate(unsigned_t idx, vector<unsigned_t> const & vTargetQubit, unsigned_t duration = 1);
        void addQubitRegion(unsigned_t idx, int_t q)            { _vsQubitRegion[idx].insert(q); }
        void addQubitRegion(unsigned_t idx, set<unsigned>& s)   { _vsQubitRegion[idx].merge(s); }
        void setInitialMapping(int_t pro_q, int_t phy_q)        { _vInitialMapping[pro_q] = phy_q; }
        void setInitialMapping(vector<unsigned_t> const &  vInitialMapping);
        void setFinalMapping(int_t pro_q, int_t phy_q)          { _vFinalMapping[pro_q] = phy_q; }
        void setFinalMapping(vector<unsigned_t> const &  vFinalMapping);
        void resetsetFinalMapping(){
            for(unsigned_t i = 0; i < _vFinalMapping.size(); ++i){
                _vFinalMapping[i] = -1;
            }
        }
        void resetsetQubitRegion()                              { _vsQubitRegion.clear(); _vsQubitRegion.resize(_nProgramQubit); }
        void resetsetQubitRegion(unsigned_t i)                  { _vsQubitRegion[i].clear(); }
        void setQubitRegion(vector<set<unsigned_t> >& vsQubitRegion){
            _vsQubitRegion.clear();
            for (set<unsigned_t>& s : vsQubitRegion){
                _vsQubitRegion.emplace_back(s);
            }
        }
        void setCircuitDepth(int_t t)                           { _circuitDepth = t; }

        vector<pair<unsigned_t, unsigned_t> >* pvpGateDependency()       { return &_vpGateDependency; };
        vector<set<int_t>>*                    pQubitRegion()             { return &_vsQubitRegion; };
        set<int_t>&                            qubitRegion(unsigned_t i) { return _vsQubitRegion[i]; };

        void clearSwap()                                        { _vSwapGate.clear(); }

        void printCircuit();
        void printCircuitLayout();
        void printCircuitInitialMapping();
        void printCircuitFinalMapping();
        void printDependency();

    private:
        string                       _name;
        unsigned_t                   _nProgramQubit;
        vector<Gate>                 _vGate;

        // store the layout synthesis result
        unsigned_t                   _circuitDepth;
        vector<Gate>                 _vSwapGate;
        vector<int_t>                _vInitialMapping;
        vector<unsigned_t>           _vFinalMapping;

        vector<pair<unsigned_t, unsigned_t> >   _vpGateDependency;
        vector<set<int_t>>           _vsQubitRegion;
};

MOLSQ_NAMESPACE_HPP_END

#endif // CIRCUIT_HPP
