clc;clear;
num=10;
den=[1 0.1 3];
[Ac,Bc,Cc,Dc]=tf2ss(num,den);
 
Delta_t=0.01;
[Ad,Bd,Cd,Dd]=c2dm(Ac,Bc,Cc,Dc,Delta_t);
Nc=3;Np=20;
[Phi,F,Phi_Phi,Phi_F,Phi_R,A,B,C]=mpcgain_D(Ad,Bd,Cd,Nc,Np);
Phi_Rs=Phi_F(:,3);
x0=[0;0;0];
A_cons=[-Phi;Phi];
Ymin=0;
Ymax=1;


H=2*(Phi_Phi+0.01*eye(Nc));
Pole=[0 0 0];
K_ob=acker(A',C',Pole)';

N_sim=400;
x=zeros(Nc,N_sim);
y=zeros(1,N_sim);
Delta_U=zeros(1,N_sim);
u=zeros(1,N_sim);
U=zeros(Nc,N_sim);

for kk=1:N_sim-1
    Delta_u=inv(Phi_Phi+0.01*eye(Nc))*(Phi_Rs-Phi_F*x(:,kk));
    Tem=Delta_u(1,1);
    U(:,kk+1)=U(:,kk)+Delta_u;
    Y=C*x(:,kk);
    Temp1=1;
    Temp1=Temp1*((Y>=0)&&(Y<=1));
%     for i=1:Np;
%         Temp1=Temp1*((Y(i)>=0)&&(Y(i)<=1));
%     end

    if Temp1==1;
        Delta_U(:,kk)=Tem;
        DeltaU(:,kk)=Delta_u;
    else
        f=-2*(Phi_Rs-Phi_F*x(:,kk));
        b1=-ones(Np,1)*Ymin+F*x(:,kk);
        b2=ones(Np,1)*Ymax-F*x(:,kk);
        b=[b1;b2];
        eta=QPhild(H,f,A_cons,b);
        Delta_U(:,kk)=eta(1,1);
        DeltaU(:,kk)=eta;
    end  
    x(:,kk+1)=A*x(:,kk)+B*Delta_U(:,kk)+K_ob*(y(kk)-C*x(:,kk));
    y(kk+1)=C*x(:,kk+1);
    u(kk+1)=u(kk)+Delta_U(:,kk);
end

t=0:1:N_sim-1;

figure(1)
plot(t,y);hold on
plot([0,N_sim-1],[Ymin,Ymin],'--','color','black')
plot([0,N_sim-1],[Ymax,Ymax],'--','color','black')
plot([0,N_sim-1],[1,1],'--','color','red')
xlabel('Sampling Instant')
ylabel('Output y')
title('0 to 1')
savefig('Eg22002.fig')
print('-dpng','Eg22002.png')