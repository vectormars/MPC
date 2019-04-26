clc;clear;
n11=12.8*conv([-1 4],[-1 4]);
d11=conv([16.7 1],conv([1 4],[1 4]));
n12=-1.89;
d12=[21 1];
n21=0;
d21=1;
n22=-19.4*conv([-3 4],[-3 4]);
d22=conv([14.4 1],conv([3 4],[3 4]));
h=1; 
Gs=tf({n11 n12;n21 n22},{d11 d12;d21 d22},h);
Gsmin=ss(Gs,'min');
[Ap,Bp,Cp,Dp]=ssdata(Gsmin);
[n1,n_in]=size(Bp);
[A,B,C]=format_change(Ap,Bp,Cp);
Q=C'*C;
R=eye(n_in,n_in);
alpha=1.2;
a=[0 0];
N=[6 6];
Np=140;
[K,P,E]=dlqr(A,B,Q,R);

gamma=1/alpha;
Qa=gamma^2*Q+(1-gamma^2)*P;
Ra=gamma^2*R;

A_hat=A/alpha;
B_hat=B/alpha;
A_hat_t=A_hat';
[Omega,Psi]=dmpc(A/alpha,B/alpha,a,N,Np,Qa,Ra);
cond(Omega)
