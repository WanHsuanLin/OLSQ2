OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
creg c[16];

// moment 0
x q[7];
x q[12];
x q[15];
x q[1];
x q[4];
x q[9];
x q[5];
cx q[11], q[10];
cx q[3], q[2];
x q[6];
x q[14];
x q[0];

// moment 1
x q[7];
x q[15];
cx q[8], q[12];
cx q[3], q[2];
cx q[9], q[5];
cx q[13], q[14];

// moment 2
x q[10];
x q[11];
x q[15];
x q[2];
cx q[3], q[7];
cx q[5], q[4];
cx q[9], q[8];
x q[0];

// moment 3
cx q[11], q[10];
x q[8];
cx q[3], q[2];
cx q[9], q[5];
cx q[15], q[14];

// moment 4
x q[15];
x q[5];
x q[10];
x q[4];
cx q[3], q[2];
cx q[9], q[8];

// measurement
measure q[4]->c[0];
measure q[15]->c[1];
measure q[7]->c[2];
measure q[11]->c[3];
measure q[14]->c[4];
measure q[13]->c[5];
measure q[5]->c[6];
measure q[6]->c[7];
measure q[3]->c[8];
measure q[10]->c[9];
measure q[8]->c[10];
measure q[12]->c[11];
measure q[2]->c[12];
measure q[1]->c[13];
measure q[9]->c[14];
measure q[0]->c[15];

// Swap objective is swap
// Mode is normal
// Device is grid
// We used sabre
// Swap duarion is 1
// And the output mode is default
// The check result is False
// Our resulting depth matches the known ground truth