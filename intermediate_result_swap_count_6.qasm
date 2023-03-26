OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
creg c[16];

// moment 0
rzz(pi/4) q[4], q[0];
rzz(pi/4) q[4], q[8];
rzz(pi/4) q[15], q[14];
rzz(pi/4) q[14], q[13];
rzz(pi/4) q[12], q[8];

// moment 1
rzz(pi/4) q[5], q[1];
rzz(pi/4) q[0], q[1];
rzz(pi/4) q[0], q[4];
rzz(pi/4) q[15], q[11];
rzz(pi/4) q[14], q[13];
rzz(pi/4) q[8], q[4];
rzz(pi/4) q[12], q[13];
rzz(pi/4) q[13], q[9];
rzz(pi/4) q[1], q[2];
rzz(pi/4) q[11], q[7];

// moment 2
rzz(pi/4) q[14], q[10];
rzz(pi/4) q[9], q[10];
rzz(pi/4) q[8], q[4];
rzz(pi/4) q[11], q[7];
rzz(pi/4) q[3], q[7];
rzz(pi/4) q[3], q[2];
rzz(pi/4) q[7], q[6];
rzz(pi/4) q[6], q[10];
rzz(pi/4) q[6], q[2];
SWAP q[4], q[5];
SWAP q[9], q[13];

// moment 3
SWAP q[8], q[12];
SWAP q[3], q[7];

// moment 4
SWAP q[14], q[15];

// moment 5

// moment 6

// moment 7
SWAP q[8], q[9];

// moment 8

// moment 9

// moment 10

// measurement
measure q[14]->c[0];
measure q[13]->c[1];
measure q[9]->c[2];
measure q[6]->c[3];
measure q[5]->c[4];
measure q[15]->c[5];
measure q[7]->c[6];
measure q[4]->c[7];
measure q[2]->c[8];
measure q[8]->c[9];
measure q[1]->c[10];
measure q[12]->c[11];
measure q[10]->c[12];
measure q[11]->c[13];
measure q[3]->c[14];
measure q[0]->c[15];
