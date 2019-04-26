clc;clear;
Am=[1 1;0 1];
Bm=[0.5;1];
Cm=[1 0];
[A,B,C]=format_change(Am,Bm,Cm);

Q=C'*C;
a=0;
R=0;
N=4;
Np=10;
[Omega,Psi]=dmpc(A,B,a,N,Np,Q,R);

