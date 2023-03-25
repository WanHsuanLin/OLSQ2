#include <cstdio>
#include "cir/circuit.hpp"
#include "device/device.hpp"
#include "olsq/olsq.hpp"
#include "misc/global.hpp"

using namespace OLSQ_NAMESPACE;

int main(int argc, char *argv[]) {
  printf("hello world!\n");
  Circuit cir("test", 4, 4);
  vector<unsigned_t> vTargetQubit(2,0);
  vTargetQubit[0] = 0;
  vTargetQubit[1] = 1;
  cir.addGate("cx", vTargetQubit);
  vTargetQubit[0] = 0;
  vTargetQubit[1] = 2;
  cir.addGate("cx", vTargetQubit);
  vTargetQubit[0] = 0;
  vTargetQubit[1] = 3;
  cir.addGate("cx", vTargetQubit);
  vTargetQubit[0] = 1;
  vTargetQubit[1] = 2;
  cir.addGate("cx", vTargetQubit);
  cir.printCircuit();
  vector<pair<unsigned_t, unsigned_t> > vEdge;
  vEdge.emplace_back(make_pair(0,1));
  vEdge.emplace_back(make_pair(1,2));
  vEdge.emplace_back(make_pair(2,3));
  vEdge.emplace_back(make_pair(0,3));
  Device device("test", 4, 4);
  device.setEdge(vEdge);
  device.printDevice();
  OLSQ olsq(cir, device);
  olsq.run("test.qasm");
  return 0;
}
