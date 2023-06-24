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
#include <queue>

OLSQ_NAMESPACE_HPP_START

class OLSQ{
    public:
        OLSQ(Circuit& cir, Device& device)
        : _pCircuit(&cir), _device(device), _swapIdx(cir.nGate()), _verbose(2), _outputQASMFile("out.qasm"){
            // cout << "+++++++++++++++++++++=" << &cir << endl;
            // cout << "+++++++++++++++++++++=" << _pCircuit << endl;
            _vpGateDependency.clear();
            _vpGateTimeWindow.clear();
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
        void setSwapDuration(unsigned_t d)                      { _olsqParam.swap_duration = d; }                 
        void setSabreForSwap(bool s, unsigned_t bound)          { _olsqParam.is_use_SABRE_for_swap = s; _olsqParam.sabre_swap_bound = bound; }                 
        void enableGateTimeWindow()                             { _olsqParam.use_window_range_for_gate = 1; }                 
        void setOptimizeForSwap(){ 
            _olsqParam.is_optimize_swap = true; 
        }                 
        void initializeTransitionMode(unsigned_t min_depth = 1){
            _olsqParam.is_transition = true;
            _olsqParam.is_given_depth = true;
            _olsqParam.min_depth = min_depth;
            _olsqParam.max_depth_bit = ceil(log2(min_depth + 1));
            _olsqParam.max_depth = pow(2, _olsqParam.max_depth_bit);
        }
        void initializeNormalMode(unsigned_t min_depth = 0){
            _olsqParam.is_transition = false;
            if (min_depth == 0){
                _olsqParam.is_given_depth = false;
            }
            else{
                _olsqParam.min_depth = min_depth;
                _olsqParam.max_depth_bit = ceil(log2(min_depth + 1));
                _olsqParam.max_depth = pow(2, _olsqParam.max_depth_bit);
            }
        }
        void run(string const & fileName);
        void dump();
        void setDependency(vector<pair<unsigned_t, unsigned_t> > & vDependencies);
        void printDependency();

        void outputQASM(string const & fileName);
        string outputQASMStr();
    
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
            bool         use_window_range_for_gate     = false;
            unsigned_t   max_depth                     = 8;  //  always (power of 2), for bit length 
            unsigned_t   max_depth_bit                 = 3;  
            unsigned_t   min_depth                     = 1;
            unsigned_t   max_depth_expand_factor       = 1;
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
                bitwuzla_set_option(pSolver, BITWUZLA_OPT_SEED, 0);
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
                bitwuzla_set_option(pSolver, BITWUZLA_OPT_SEED, 0);
                bitwuzla_set_termination_callback(pSolver, isTimeout, &state);
                vvPi.clear();
                vTg.clear();
                vvSigma.clear();
            }
            Bitwuzla *                              pSolver = bitwuzla_new();
            vector<vector<const BitwuzlaTerm*> >    vvPi;         // t->qId
            vector<const BitwuzlaTerm*>             vTg;
            vector<vector<const BitwuzlaTerm*> >   vvSigma;      // t->qId
            pair<TimeState, unsigned_t>             state;
        } _smt;

    struct Node {
            Node(int_t p, int_t idx, int_t qIdx, int_t eIdx, int_t dis): parentIdx(p), parentEIdx(eIdx), idx(idx), qIdx(qIdx), dis(dis){};
            ~Node() {};
            int_t parentIdx;
            int_t parentEIdx;
            int_t idx;
            int_t qIdx;
            int_t dis;
        };
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
        vector<pair<unsigned_t, unsigned_t> >   _vpGateTimeWindow; // pair<start time, end time>
        string                                  _outputQASMFile;
    ////////////////////////////
    // Private functions
    ////////////////////////////
    private:
        void runSMT();

        void generateFormulationZ3();
        void constructVariableZ3();
        void addInjectiveMappingConstraintsZ3(unsigned_t boundOffset = 0);
        void addValidTwoQubitGateConstraintsZ3(unsigned_t boundOffset = 0); // if boundOffset > 0, it means we increase min_depth so we need to add constraints for partial t
        void addDependencyConstraintsZ3();
        void addSwapConstraintsZ3(unsigned_t boundOffset = 0);
        void addTransformationConstraintsZ3(unsigned_t boundOffset = 0);
        void addDepthConstraintsZ3();
        void addSwapCountConstraintsZ3(unsigned_t bound);
        bool checkModel();

        bool optimize();
        bool optimizeDepth();
        bool optimizeSwap();
        bool optimizeSwapForDepth(unsigned_t lower_swap_bound, unsigned_t upper_swap_bound, bool firstRun);

        void extractModel();
        void asapScheduling();

        void increaseDepthBound();
        void constructDependency();
        void addDependency(Gate& g1, Gate& g2){
            _vpGateDependency.emplace_back(make_pair(g1.idx(), g2.idx()));
        }
        void addDependency(unsigned_t g1, unsigned_t g2){
            _vpGateDependency.emplace_back(make_pair(g1, g2));
        }
        unsigned_t extract_longest_chain();
        void updateSMT(unsigned_t d);

        // construct overconstrained problem
        void constructGateTimeWindow();
        void updateGateTimeWindow(unsigned_t d);
        void printGateTimeWindow();
};
OLSQ_NAMESPACE_HPP_END

#endif // OLSQ_HPP
