OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
creg c[16];

// moment 0
rzz(pi/4) q[14], q[10];

// moment 1
SWAP q[10], q[11];

// moment 2
rzz(pi/4) q[14], q[10];

// moment 3
rzz(pi/4) q[14], q[15];

// moment 4
rzz(pi/4) q[11], q[15];

// moment 5
rzz(pi/4) q[11], q[10];
SWAP q[14], q[15];

// moment 6
rzz(pi/4) q[14], q[10];

// measurement
measure q[15]->c[0];
measure q[11]->c[1];
measure q[14]->c[2];
measure q[10]->c[3];

// The check result is True
// Swap objective is swap
// Mode is transition
// Device is grid
// We used sabre
// Swap duarion is 1
// And the output mode is default