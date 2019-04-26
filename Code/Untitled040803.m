clc;clear;
num11=12.8*conv([-1 4],[-1 4]);
den11=conv([16.7 1],conv([1 4],[1 4]));
num12=-1.89*conv([-3 4],[-3 4]);
den12=conv([21 1],conv([3 4],[3 4]));
num21=1.28*conv([-7 4],[-7 4]);
den21=conv([10.9 1],conv([7 4],[7 4]));
num22=-19.4*conv([-3 4],[-3 4]);
den22=conv([14.4 1],conv([3 4],[3 4]));
Gs=tf({num11 num12; num21 num22},{den11 den12;den21 den22});
Gsmin=ss(Gs,'min');
[Ac,Bc,Cc,Dc]=ssdata(Gs);
Delta_t=1;
[Ap,Bp,Cp,Dp]=c2dm(Ac,Bc,Cc,Dc,Delta_t,'zoh');
[A,B,C]=format_change(Ap,Bp,Cp);

Q=C'*C;
a1=0;
a2=0;
a=[a1 a2];
R=eye(2);
N1=5;
N2=5;
N=[N1 N2];
Np=140;
alpha=1.2;
[K,P,Z]=dlqr(A,B,Q,R);
gamma=1/alpha;
g=gamma*gamma;
Qa=g*Q+(1-g)*P;
Ra=g*R;
[Omega,Psi]=dmpc(A/alpha,B/alpha,a,N,Np,Qa,Ra);
cond(Omega)


