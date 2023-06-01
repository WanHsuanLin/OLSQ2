// parameters: output_mode:IR objective:swap mode:transition device:grid sabre:true swap_duration:1 device_size:4 verify_result:true

OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
rzz(pi/4) q[0],q[1];
rzz(pi/4) q[0],q[3];
rzz(pi/4) q[0],q[2];
rzz(pi/4) q[1],q[2];
rzz(pi/4) q[1],q[3];
rzz(pi/4) q[2],q[3];
