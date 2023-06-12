#include "misc/global.hpp"
#include "olsq/olsq.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>


OLSQ_NAMESPACE_CPP_START

namespace py = pybind11;


PYBIND11_MODULE(olsqpyb, m) {
    // optional module docstring
    m.doc() = "pybind11 olsq plugin";

    // bindings to Circuit class
    py::class_<Circuit>(m, "Circuit")
        .def(py::init<const string &, unsigned_t, unsigned_t>())
        // .def("addGate", [](Circuit& self, const string& gateName, pybind11::list pyTargetQubit, unsigned_t duration){
        //         // convert pyTargetQubit into vTargetQubit
        //         vector<unsigned_t> vTargetQubit;
        //         for (py::handle item : pyTargetQubit){
        //             vTargetQubit.emplace_back(item.cast<unsigned_t>());
        //         }
        //         return self.addGate(gateName, vTargetQubit, duration); 
        //     }, "Add a gate to circuit", py::arg("gateName"), py::arg("pyTargetQubit"), py::arg("duration") = 1)
        //     vector<unsigned_t> vInitialMapping;
        //     for (py::handle item : pyInitialMapping){
        //         vInitialMapping.emplace_back(item.cast<unsigned_t>());
        //     }
        //     return self.setInitialMapping(vInitialMapping); 
        // }, "Set initial mapping")
        .def("printCircuit", &Circuit::printCircuit, "Print circuit information")
        .def("printCircuitLayout", &Circuit::printCircuitLayout, "Print compiled circuit information");
    m.def("setMappingRegion", [](Circuit& self, pybind11::list pyMappingRegion){
        vector<set<int_t> > vsMappingRegion;
        for (py::handle item : pyMappingRegion){
            vsMappingRegion.emplace_back(item.cast<set<int_t>>());
        }
        vsMappingRegion.resize(self.nProgramQubit());
        return self.setQubitRegion(vsMappingRegion); 
    }, "Set mapping region for circuit");

    m.def("addGate", [](Circuit& self, const string& gateName, pybind11::list pyTargetQubit, unsigned_t duration){
        // convert pyTargetQubit into vTargetQubit
        vector<unsigned_t> vTargetQubit;
        for (py::handle item : pyTargetQubit){
            vTargetQubit.emplace_back(item.cast<unsigned_t>());
        }
        return self.addGate(gateName, vTargetQubit, duration); 
    }, "Add a gate to circuit", py::arg("self"), py::arg("gateName"), py::arg("pyTargetQubit"), py::arg("duration") = 1);
    m.def("setInitialMapping", [](Circuit& self, pybind11::list pyInitialMapping){
        vector<unsigned_t> vInitialMapping;
        for (py::handle item : pyInitialMapping){
            vInitialMapping.emplace_back(item.cast<unsigned_t>());
        }
        return self.setInitialMapping(vInitialMapping); 
    }, "Set initial mapping");

    // bindings to Device class
    py::class_<Device>(m, "Device")
        .def(py::init<const string &, unsigned_t, unsigned_t>())
        .def("nQubit", &Device::nQubit, "get number of physical qubits")
        .def("addEdge", &Device::addEdge, "Add an edge to device")
        // .def("setEdge", [](Device& self, pybind11::list pyEdge){
        //     // convert pyEdge into vEdge
        //     vector<pair<unsigned_t, unsigned_t> > vEdge;
        //     for (py::handle item : pyEdge){
        //         vEdge.emplace_back(item.cast<pair<unsigned_t, unsigned_t> >());
        //     }
        //     return self.setEdge(vEdge); 
        // }, "Set device edges")
        .def("printDevice", &Device::printDevice, "Print device information");
    m.def("setEdge", [](Device& self, pybind11::list pyEdge){
        // convert pyEdge into vEdge
        vector<pair<unsigned_t, unsigned_t> > vEdge;
        for (py::handle item : pyEdge){
            vEdge.emplace_back(item.cast<pair<unsigned_t, unsigned_t> >());
        }
        return self.setEdge(vEdge); 
    }, "Set device edges");

    // bindings to OLSQ class
    py::class_<OLSQ>(m, "OLSQ")
        .def(py::init<Circuit&, Device&>())
        .def("setSwapDuration", &OLSQ::setSwapDuration, "use given initial mapping (default = 1)")
        .def("setSabreForSwap", &OLSQ::setSabreForSwap, "set SABRE swap count as SWAP upper bound")
        .def("initializeTransitionMode", &OLSQ::initializeTransitionMode, "use transition based mode (default)", py::arg("min_depth") = 1)
        .def("initializeNormalMode", &OLSQ::initializeNormalMode, "use normal based mode", py::arg("min_depth") = 0)
        .def("setOptimizeForSwap", &OLSQ::setOptimizeForSwap, "Set optimization object to swap count")
        .def("run", &OLSQ::run, "run quantum layout synthesis")
        // .def("setDependency", [](OLSQ& self, pybind11::list pyDependencies){
        //     // convert pyDependencies into vDependencies
        //     vector<pair<unsigned_t, unsigned_t> > vDependencies;
        //     for (py::handle item : pyDependencies){
        //         vDependencies.emplace_back(item.cast<pair<unsigned_t, unsigned_t> >());
        //     }
        //     return self.setDependency(vDependencies); 
        // }, "Set dependency for circuit")
        .def("printDependency", &OLSQ::printDependency, "Print dependency information");
    m.def("setDependency", [](OLSQ& self, pybind11::list pyDependencies){
        // convert pyDependencies into vDependencies
        vector<pair<unsigned_t, unsigned_t> > vDependencies;
        for (py::handle item : pyDependencies){
            vDependencies.emplace_back(item.cast<pair<unsigned_t, unsigned_t> >());
        }
        return self.setDependency(vDependencies); 
    }, "Set dependency for circuit");
    
}

OLSQ_NAMESPACE_CPP_END
