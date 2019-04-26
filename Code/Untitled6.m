clc;clear;
n11=12.8*conv([-1 4],[-1 4]);
d11=conv([16.7 1],conv([1 4],[1 4]));
n12=-1.89*conv([-3 4],[-3 4]);
d12=conv([21 1],conv([3 4],[3 4]));
n21=1.28*conv([-7 4],[-7 4]);
d21=conv([10.9 1],conv([7 4],[7 4]));
n22=-19.4*conv([-3 4],[-3 4]);
d22=conv([14.4 1],conv([3 4],[3 4]));

h=1; %sampling interval
Gs=tf({n11 n12;n21 n22},{d11 d12;d21 d22});
Gsmin=ss(Gs,'min');
[Ac,Bc,Cc,Dc]=ssdata(Gsmin);
[Ap,Bp,Cp,Dp]=c2dm(Ac,Bc,Cc,Dc,h);
[A,B,C]=format_change(Ap,Bp,Cp);
[n1,n_in]=size(Bp);
a=[0 0];N=[5 5];
Np=140;
R=eye(n_in,n_in);Q=C'*C;

alpha=1.2;
[K,P,Z]=dlqr(A,B,Q,R);
gamma=1/alpha;
g=gamma*gamma;
Qa=g*Q+(1-g)*P;
Ra=g*R;
A_hat=A/alpha;
B_hat=B/alpha;
[Omega,Psi]=dmpc(A_hat,B_hat,a,N,Np,Qa,Ra);
cond(Omega)