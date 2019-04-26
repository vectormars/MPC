clc;clear;
num=10;
den=[1 0.1 3];
[Ac,Bc,Cc,Dc]=tf2ss(num,den);

Delta_t=0.01;
[Ad,Bd,Cd,Dd]=c2dm(Ac,Bc,Cc,Dc,Delta_t);
Nc=3;Np=20;
[Phi_Phi,Phi_F,Phi_R,A,B,C]=mpcgain(Ad,Bd,Cd,Nc,Np);
Phi_Rs=Phi_F(:,3);
x0=[0;0;0];
A_cons=[1 0 0;0 1 0;0 0 1;-1 0 0;0 -1 0;0 0 -1];
H=2*(Phi_Phi+0.01*eye(Nc));
Pole=[0 0 0];
K_ob=acker(A',C',Pole)';

N_sim=61;
x=zeros(Nc,N_sim);
y=zeros(1,N_sim);
Delta_U=zeros(1,N_sim);
u=zeros(1,N_sim);
Delta_u_storage=zeros(Nc,N_sim);
%%%%%%%%%%%%%%%% with constraint u1,u2,u3 %%%%%%%%%%%%%%%%
for kk=1:N_sim-1
    Delta_u=(Phi_Phi+0.01*eye(Nc))\(Phi_Rs-Phi_F*x(:,kk));
    Tem=Delta_u;
    Temp=1;
    for i=1:Nc;
        Temp=Temp*((Delta_u(i,1)>=-1.5)&&(Delta_u(i,1)<=3));
    end
    if Temp==1;
        Delta_U(:,kk)=Tem(1,1);
        Delta_u_storage(:,kk)=Tem;
    else
        f=-2*(Phi_Rs-Phi_F*x(:,kk));
        b=[3;3;3;1.5;1.5;1.5];
        eta=QPhild(H,f,A_cons,b);
        Delta_U(:,kk)=eta(1,1);
        Delta_u_storage(:,kk)=eta;
    end
    x(:,kk+1)=A*x(:,kk)+B*Delta_U(:,kk)+K_ob*(y(kk)-C*x(:,kk));
    y(kk+1)=C*x(:,kk+1);
    u(kk+1)=u(kk)+Delta_U(:,kk);
end

figure(1)
t=0:60;
plot(t,Delta_u_storage(1,:),t,Delta_u_storage(2,:),t,Delta_u_storage(3,:));hold on
plot([0,60],[-1.5,-1.5],'--','color','black')
plot([0,60],[3,3],'--','color','black')
legend('\Deltau_k','\Deltau_{k+1}','\Deltau_{k+2}')
savefig('Eg21501.fig')
print('-dpng','Eg21501.png')

%%%%%%%%%%%%%%%% with constraint u1 only %%%%%%%%%%%%%%%%
b1=[3;1.5];
A_cons1=[1 0 0;-1 0 0];

x1=zeros(Nc,N_sim);
y1=zeros(1,N_sim);
Delta_U1=zeros(1,N_sim);
u1=zeros(1,N_sim);

for kk=1:N_sim-1
    Delta_u=(Phi_Phi+0.01*eye(Nc))\(Phi_Rs-Phi_F*x1(:,kk));
    Tem=Delta_u(1,1);
    if ((Tem>=-1.5)&&(Tem<=3))
        Delta_U1(:,kk)=Tem;
    else
        f=-2*(Phi_Rs-Phi_F*x1(:,kk));
        eta=QPhild(H,f,A_cons1,b1);
        Delta_U1(:,kk)=eta(1,1);
    end
    x1(:,kk+1)=A*x1(:,kk)+B*Delta_U1(:,kk)+K_ob*(y1(kk)-C*x1(:,kk));
    y1(kk+1)=C*x1(:,kk+1);
    u1(kk+1)=u1(kk)+Delta_U1(:,kk);
end

figure(2)
subplot(211)
plot(t,y)
hold on
plot(t,y1,'r')
xlabel('Sampling Instant')
ylabel('Output y')
legend('with constraint \Deltau1,\Deltau2,\Deltau3','with constraint \Deltau1 only',4)
subplot(212)
stairs(t,u)
hold on
stairs(t,u1,'r')
xlabel('Sampling Instant')
ylabel('Control signal u')
legend('with constraint \Deltau1,\Deltau2,\Deltau3','with constraint \Deltau1 only')
savefig('Eg21502.fig')
print('-dpng','Eg21502.png')

figure(3)
plot(t,Delta_U)
hold on
plot(t,Delta_U1,'r')
plot([0,60],[-1.5,-1.5],'--','color','black')
plot([0,60],[3,3],'--','color','black')
xlabel('Sampling Instant')
ylabel('Delta u')
legend('with constraint \Deltau1,\Deltau2,\Deltau3','with constraint \Deltau1 only')
savefig('Eg21503.fig')
print('-dpng','Eg21503.png')
