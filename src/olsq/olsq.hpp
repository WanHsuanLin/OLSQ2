/***********************************************************************
  File        [ olsq.hpp ]
  System      [ OLSQ: optimal quantum layout synthesis tool]
  Package     [ olsq ]
  Synopsis    [ OLSQ class header ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 22, Nov., 2022 ]
***********************************************************************/
#ifndef OLSQ_HPP
#define OLSQ_HPP

#include "misc/global.hpp"
#include "cir/circuit.hpp"
#include "device/device.hpp"
#include "misc/timeUsage.hpp"
#include <algorithm>
#include <bitwuzla/bitwuzla.h>
#include <pb2cnf.h>
#include <set>
#include <map>

OLSQ_NAMESPACE_HPP_START

class OLSQ{
    public:
        OLSQ(Circuit& cir, Device& device)
        : _pCircuit(&cir), _device(device), _swapIdx(cir.nGate()), _verbose(2), _outputQASM("out.qasm"){
            // cout << "+++++++++++++++++++++=" << &cir << endl;
            // cout << "+++++++++++++++++++++=" << _pCircuit << endl;
            _vpGateDependency.clear();
        }
        ~OLSQ() {}
        bool                          isValidGateIdx(unsigned_t idx)       const { return 0 <= idx && idx < _pCircuit->nGate();}
        bool                          isValidDependencyIdx(unsigned_t idx) const { return 0 <= idx && idx < _vpGateDependency.size();}
        pair<unsigned_t, unsigned_t>& dependency(unsigned_t idx)                 { assert(isValidDependencyIdx(idx)); return _vpGateDependency[idx]; }        
        
        void reset()                                            { _smt.reset(_olsqParam.timeout, 0);  
                                                                    if(_olsqParam.is_transition)
                                                                        initializeTransitionMode(); 
                                                                    else
                                                                        initializeNormalMode();
                                                                    if(_olsqParam.is_optimize_swap)
                                                                        setOptimizeForSwap();}
        void setVerbose(unsigned_t verbose)                     { _verbose = verbose; }
        void setCircuit(Circuit & cir)                          { _pCircuit = &cir; }
        void useCircuitInitalMapping()                          { _olsqParam.is_given_initial_mapping = true; }                 
        void setSwapDuration(unsigned_t d)                      { _olsqParam.swap_duration = d; }                 
        void setSabreForSwap(bool s, unsigned_t bound)          { _olsqParam.is_use_SABRE_for_swap = s; _olsqParam.sabre_swap_bound = bound; }                 
        void setOptimizeForSwap(){ 
            _olsqParam.is_optimize_swap = true; 
            _olsqParam.min_depth = _olsqParam.max_depth - 1;
        }                 
        void initializeTransitionMode(unsigned_t max_depth = 5, unsigned_t min_depth = 1){
            _olsqParam.is_transition = true;
            _olsqParam.is_given_depth = true;
            _olsqParam.max_depth = max_depth;
            _olsqParam.min_depth = min_depth;
            _olsqParam.min_depth_expand_step = 1;
            _olsqParam.max_depth_expand_factor = 2;
        }
        void initializeNormalMode(unsigned_t max_depth = 0, unsigned_t min_depth = 0){
            _olsqParam.is_transition = false;
            if (min_depth == 0){
                _olsqParam.is_given_depth = false;
            }
            else{
                _olsqParam.min_depth = min_depth;
            }
            if (max_depth == 0){
                _olsqParam.max_depth = 2 * _olsqParam.min_depth;
            }
            else{
                _olsqParam.max_depth = max_depth;
            }
            _olsqParam.max_depth = max_depth;
            _olsqParam.min_depth_expand_step = 5;
            _olsqParam.max_depth_expand_factor = 2;
        }
        void run(string const & fileName);
        void outputQASM(string const & fileName);
        void setDependency(vector<pair<unsigned_t, unsigned_t> > & vDependencies);
        void printDependency();
    
    ////////////////////////////
    // Struct OLSQParam
    ////////////////////////////
    private:
        struct OLSQParam {
            bool         is_transition                 = true;
            bool         is_optimize_swap              = false;
            bool         is_use_SABRE_for_swap         = false;
            bool         is_given_dependency           = false;
            bool         is_given_depth                = false;
            bool         is_given_initial_mapping      = false;
            unsigned_t   max_depth                     = 5;
            unsigned_t   min_depth                     = 1;
            unsigned_t   max_depth_expand_factor       = 2;
            unsigned_t   min_depth_expand_step         = 1;
            unsigned_t   sabre_swap_bound              = 0;
            unsigned_t   timeout                       = 86400;
            unsigned_t   swap_duration                  = 1;
        } _olsqParam;

    ////////////////////////////
    // Struct OLSQParam
    ////////////////////////////
     struct smt {
            smt(){
                // z3::context c;
                // solver = z3::solver(c);
                bitwuzla_set_option(pSolver, BITWUZLA_OPT_PRODUCE_MODELS, 1);
                bitwuzla_set_option(pSolver, BITWUZLA_OPT_INCREMENTAL, 1);
                vvPi.clear();
                vTg.clear();
                vvSigma.clear();
            }
            ~smt(){
                bitwuzla_delete(pSolver);
            }

            static int_t isTimeout(void* data){
                pair<TimeState, unsigned_t>* pState = ( pair<TimeState, unsigned_t>* ) data;
                TimeState& fullStart = pState->first;
                unsigned_t timeout = pState->second;
                TimeState curSt;
                curSt.checkUsage();
                TimeState dur = diff(fullStart, curSt);
                return (timeout < dur.userTime);
            }

            void reset(unsigned_t timeout, unsigned_t iter){
                if(iter == 0){
                    state.first.checkUsage();
                    state.second = timeout;
                }
                bitwuzla_reset(pSolver);
                bitwuzla_set_option(pSolver, BITWUZLA_OPT_PRODUCE_MODELS, 1);
                bitwuzla_set_option(pSolver, BITWUZLA_OPT_INCREMENTAL, 1);
                bitwuzla_set_termination_callback(pSolver, isTimeout, &state);
                vvPi.clear();
                vTg.clear();
                vvSigma.clear();
            }
            Bitwuzla *                              pSolver = bitwuzla_new();
            vector<vector<const BitwuzlaTerm*> >    vvPi;         // t->qId
            vector<const BitwuzlaTerm*>             vTg;
            vector<vector<const BitwuzlaTerm*> >    vvSigma;      // t->qId
            pair<TimeState, unsigned_t>             state;
        } _smt;
    ////////////////////////////
    // Private member
    ////////////////////////////
    private:
        TimeUsage                               _timer;
        Circuit*                                _pCircuit;
        Device&                                 _device;
        vector<pair<unsigned_t, unsigned_t> >   _vpGateDependency;
        unsigned_t                              _swapIdx;
        unsigned_t                              _iter;
        unsigned_t                              _verbose;
        string                                  _outputQASM;
    ////////////////////////////
    // Private functions
    ////////////////////////////
    private:
        void runSMT();

        void generateFormulationZ3();
        void constructVariableZ3();
        void addInjectiveMappingConstraintsZ3();
        void addValidTwoQubitGateConstraintsZ3();
        void addDependencyConstraintsZ3();
        void addSwapConstraintsZ3();
        void addTransformationConstraintsZ3();
        void addDepthConstraintsZ3();
        void addSwapCountConstraintsZ3(unsigned_t bound);
        bool checkModel();

        bool optimize();
        bool optimizeDepth();
        bool optimizeSwap();
        bool optimizeSwapForDepth(unsigned_t lower_swap_bound, unsigned_t upper_swap_bound, bool firstRun);

        void extractModel();
        void asapScheduling();

        void increase_depth_bound();
        void constructDependency();
        void addDependency(Gate& g1, Gate& g2){
            _vpGateDependency.emplace_back(make_pair(g1.idx(), g2.idx()));
        }
        void addDependency(unsigned_t g1, unsigned_t g2){
            _vpGateDependency.emplace_back(make_pair(g1, g2));
        }
        unsigned_t extract_longest_chain();
};
OLSQ_NAMESPACE_HPP_END

#endif // OLSQ_HPP
