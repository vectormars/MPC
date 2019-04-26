clc;clear;
num=1;
den=[1 2*0.1*3 3*3];
[Ac,Bc,Cc,Dc]=tf2ss(num,den);
Delta_t=0.1;
[Ad,Bd,Cd,Dd]=c2dm(Ac,Bc,Cc,Dc,Delta_t);
[A,B,C]=format_change(Ad,Bd,Cd);
N=8;Np=46;Nc=15;
a=0.7;
R=0.3;
Q=C'*C;
X0=[0.1;0.2;0.3];

[m1,n1]=size(Cd);
[n1,n_in]=size(Bd);
y=zeros(m1,1);
u=zeros(n_in,1);
u(:,1)=6;
[Al,L0]=lagd(a,N);
xm=zeros(n1,1);
xm(:,1)=[0.1;0.2];
N_sim=Np;
sp=ones(1,N_sim+10);
[M,Lzerot]=Mdu(a,N,1,1);
[Omega,Psi]=dmpc(A,B,a,N,Np,Q,R);
[u1,y1,deltau1,k]=simuuc(xm,u,y,sp,Ad,Bd,Cd,N_sim,Omega,Psi,Lzerot);
t=0:1:Np-1;
stairs(t,deltau1)






