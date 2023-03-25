/***********************************************************************
  File        [ gate.hpp ]
  System      [ OLSQ: optimal quantum layout synthesis tool]
  Package     [ cir ]
  Synopsis    [ gate class header ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 22, Nov., 2022 ]
***********************************************************************/
#ifndef GATE_HPP
#define GATE_HPP

#include "misc/global.hpp"

OLSQ_NAMESPACE_HPP_START

class Gate
{
    public:
        Gate() :
            _idx(0), _name(""), _duration(0)
            {
                _vTargetProgramQubit.clear();
                _vTargetPhysicalQubit.clear();
            }
        Gate(unsigned_t idx, const string& gateName, unsigned_t duration = 1):
            _idx(idx), _name(gateName), _duration(duration)
            {
                _vTargetProgramQubit.clear();
                _vTargetPhysicalQubit.clear();
            }
        ~Gate() {}

        // set function
        void       addTargetProgramQubit(const unsigned_t q)                                { _vTargetProgramQubit.emplace_back(q); _vTargetPhysicalQubit.emplace_back(0);}
        void       addTargetPhysicalQubit(const unsigned_t idx, const unsigned_t phy_q)     { assert(isValidQubitIdx(idx)); _vTargetPhysicalQubit[idx] = phy_q; }
        void       setvProgramTargetQubit(vector<unsigned_t>& vTargetQubit)                 { _vTargetProgramQubit = move(vTargetQubit); }
        void       setDuration(const unsigned_t d)                                          { _duration = d; }
        void       setExecutionTime(const unsigned_t t)                                     { _executionTime = t; }
        
        // get function
        unsigned_t idx()                               const { return _idx;  }
        string     name()                              const { return _name; }
        unsigned_t duration()                          const { return _duration; }
        unsigned_t executionTime()                     const { return _executionTime; }
        unsigned_t nTargetQubit()                      const { return _vTargetProgramQubit.size(); }
        bool       isValidQubitIdx(unsigned_t idx)     const { return 0 <= idx && idx < _vTargetProgramQubit.size();}
        unsigned_t targetProgramQubit(unsigned_t idx)  const { assert(isValidQubitIdx(idx)); return _vTargetProgramQubit[idx]; }
        unsigned_t targetPhysicalQubit(unsigned_t idx) const { assert(isValidQubitIdx(idx)); return _vTargetPhysicalQubit[idx]; }


    private:
        unsigned_t         _idx;
        string             _name;
        vector<unsigned_t> _vTargetProgramQubit;
        unsigned_t         _duration;
        // store the layout synthesis result
        vector<unsigned_t> _vTargetPhysicalQubit;
        unsigned_t         _executionTime;

};

OLSQ_NAMESPACE_HPP_END

#endif // GATE_HPP