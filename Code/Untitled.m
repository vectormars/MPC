clc;clear;
A=[2.6 -2.24 0.64 -0.5;1 0 0 0;0 1 0 0;0 0 0 0];
B=[1 0 0 1]';
C=[1 0 0 0];

Q=C'*C;R=1;
[K,P,Z]=dlqr(A,B,Q,R);

Np=46;N=3;a=0.3;

alpha=1.2;
gamma=1/alpha;
g=gamma*gamma;
Qa=g*Q+(1-g)*P;
Ra=g*R;
[Omega,Psi]=dmpc(A/alpha,B/alpha,a,N,Np,Qa,Ra);
[Al,L0]=lagd(a,N);
Kmpc=L0'*(Omega\Psi);


