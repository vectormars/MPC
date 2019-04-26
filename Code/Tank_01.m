clc;clear;
A=1;
B=1;
C=0.1;

Nc=3;Np=20;
[Phi,F]=mpcgain(A,B,C,Nc,Np);
P_P=Phi'*Phi;
P_F=Phi'*F;
P_Rs=P_F(:,end);
H=2*(P_P+0.01*eye(Nc));
Pole=0;
K_ob=acker(A',C',Pole)';

N_sim=61;
x=zeros(Nc,N_sim);
y=zeros(1,N_sim);
Delta_U=zeros(1,N_sim);
u=zeros(1,N_sim);
U=zeros(Nc,N_sim);

%%%%%%%%%%%%%%%% with constraint 2.18 %%%%%%%%%%%%%%%%
x=zeros(Nc,N_sim);
y=zeros(1,N_sim);
Delta_U=zeros(1,N_sim);
u=zeros(1,N_sim);
U=zeros(Nc,N_sim);
A_cons=[1 0 0;-1 0 0;1 0 0;-1 0 0];
for kk=1:N_sim-1
    Delta_u=(P_P+0.01*eye(Nc))/(P_Rs-P_F*x(:,kk));
    Tem=Delta_u(1,1);
    U(:,kk+1)=U(:,kk)+Delta_u;
    Temp1=1;
    for i=1:Nc;
        Temp1=Temp1*((U(i,kk+1)>=0)&&(U(i,kk+1)<=10));
    end
    Temp2=((Tem>=-1.5)&&(Tem<=3));
    Temp=Temp1*Temp2;
    if Temp==1;
        Delta_U(:,kk)=Tem;
        DeltaU(:,kk)=Delta_u;
    else
        f=-2*(P_Rs-P_F*x(:,kk));
        b=[3;1.5;6-u(kk);3+u(kk)];
        eta=QPhild(H,f,A_cons,b);
        Delta_U(:,kk)=eta(1,1);
        DeltaU(:,kk)=eta;
    end  
    x(:,kk+1)=A*x(:,kk)+B*Delta_U(:,kk)+K_ob*(y(kk)-C*x(:,kk));
    y(kk+1)=C*x(:,kk+1);
    u(kk+1)=u(kk)+Delta_U(:,kk);
end
