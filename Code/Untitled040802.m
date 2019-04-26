clc;clear;
Am=[1 1;0 1];
Bm=[0.5;1];
Cm=[1 0];
[A,B,C]=format_change(Am,Bm,Cm);

Q=C'*C;
a=0;
R=1;
N=4;
Np=150;
alpha=1;
[K,P,Z]=dlqr(A,B,Q,R);

gamma=1/alpha;
g=gamma*gamma;
Qa=g*Q+(1-g)*P;
Ra=g*R;

[Omega,Psi]=dmpc(A/alpha,B/alpha,a,N,Np,Qa,Ra);
G=Omega\Psi;
cond(Omega)