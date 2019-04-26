clc;clear;
num=[-5.7980 19.5128 -21.6452 7.9547];
den=[1 -3.0228 3.8630 -2.6426 0.8084];
[Ap,Bp,Cp,Dp]=tf2ss(num,den);
[A,B,C]=format_change(Ap,Bp,Cp);

Q=C'*C;
R=0.3;
[K,P,Z]=dlqr(A,B,Q,R);
X0=[0.1 0.2 0.3 0.4 0.5]';

Np=46;N=10;
a=0.4;

alpha=1.2;
gamma=1/alpha;
g=gamma*gamma;
Qa=g*Q+(1-g)*P;
Ra=g*R;
[Omega,Psi]=dmpc(A/alpha,B/alpha,a,N,Np,Qa,Ra);
[Al,L0]=lagd(a,N);
Kmpc=L0'*(Omega\Psi);

