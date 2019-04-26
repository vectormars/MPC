clc;clear;
Ap=[1.4 -0.48 1;1 0 0;0 0 0];
Bp=[0;0;1];
Cp=[1 0 0];
[A,B,C]=format_change(Ap,Bp,Cp);

Q=C'*C;
R=0.1;
Np=46;N=8;
a=0.5;
alpha=1;

[K,P,Z]=dlqr(A,B,Q,R);
gamma=1/alpha;
g=gamma*gamma;
Qa=g*Q+(1-g)*P;
Ra=g*R;
[Omega,Psi]=dmpc(A/alpha,B/alpha,a,N,Np,Qa,Ra);
cond(Omega)