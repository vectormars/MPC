clc;clear;
num=10;
den=[1 0.1 3];
[Ac,Bc,Cc,Dc]=tf2ss(num,den);
 
Delta_t=0.01;
[Ad,Bd,Cd,Dd]=c2dm(Ac,Bc,Cc,Dc,Delta_t);
[A,B,C]=format_change(Ad,Bd,Cd);
Nc=3;Np=20;
[Phi,F]=mpcgain(A,B,C,Nc,Np);
P_P=Phi'*Phi;
P_F=Phi'*F;
P_Rs=P_F(:,3);
H=2*(P_P+0.01*eye(Nc));
Pole=[0 0 0];
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
    Delta_u=inv(P_P+0.01*eye(Nc))*(P_Rs-P_F*x(:,kk));
    Tem=Delta_u(1,1);
    U(:,kk+1)=U(:,kk)+Delta_u;
    Temp1=1;
    for i=1:Nc;
        Temp1=Temp1*((U(i,kk+1)>=-3)&&(U(i,kk+1)<=6));
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
%%%%%%%%%%%%%%%% with constraint 2.19 %%%%%%%%%%%%%%%%
x1=zeros(Nc,N_sim);
y1=zeros(1,N_sim);
Delta_U1=zeros(1,N_sim);
u1=zeros(1,N_sim);
U1=zeros(Nc,N_sim);

A_cons1=[eye(3);-eye(3);1 0 0;1 1 0;1 1 1;-1 0 0;-1 -1 0;-1 -1 -1];
for kk=1:N_sim-1
    Delta_u1=inv(P_P+0.01*eye(Nc))*(P_Rs-P_F*x1(:,kk));
    Tem=Delta_u1(1,1);
    U1(:,kk+1)=U1(:,kk)+Delta_u1;
    Temp1=1;
    for i=1:Nc;
        Temp1=Temp1*((U1(i,kk+1)>=-3)&&(U1(i,kk+1)<=6));
    end
    Temp2=1;
    for i=1:Nc;
        Temp2=Temp2*((Delta_u1(i,:))>=-1.5&&(Delta_u1(i,:)<=3));
    end
    Temp=Temp1*Temp2;
    if Temp==1;
        Delta_U1(:,kk)=Tem;
        DeltaU1(:,kk)=Delta_u1;
    else
        f=-2*(P_Rs-P_F*x1(:,kk));
        b1=[3;3;3;1.5;1.5;1.5;6-u1(kk);6-u1(kk);6-u1(kk);3+u1(kk);3+u1(kk);3+u1(kk)];
        eta=QPhild(H,f,A_cons1,b1);
        Delta_U1(:,kk)=eta(1,1);
        DeltaU1(:,kk)=eta;
    end  
    x1(:,kk+1)=A*x1(:,kk)+B*Delta_U1(:,kk)+K_ob*(y1(kk)-C*x1(:,kk));
    y1(kk+1)=C*x1(:,kk+1);
    u1(kk+1)=u1(kk)+Delta_U1(:,kk);
end

%%%%%%%%%%%%%%%%% Print %%%%%%%%%%%%%%%%%%%%%%%
t=0:1:N_sim-1;

figure(1)
subplot(211)
plot(t,y);hold on
plot(t,y1,'r')
xlabel('Sampling Instant')
ylabel('Output y')
legend('constraint 2.18','constraint 2.19',2)
subplot(212)
stairs(t,u);hold on
stairs(t,u1,'r')
xlabel('Sampling Instant')
ylabel('Control signal u')
legend('constraint 2.18','constraint 2.19')
 
figure(2)
plot(t,Delta_U);hold on
plot(t,Delta_U1,'r')
xlabel('Sampling Instant')
ylabel('Delta u')
legend('constraint 2.18','constraint 2.19')

