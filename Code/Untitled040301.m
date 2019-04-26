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
x0=[0;0;0];
b=[3;1.5];
A_cons=[1 0 0;1 1 0;1 1 1;-1 0 0;-1 -1 0;-1 -1 -1];
H=2*(P_P+0.01*eye(Nc));
Pole=[0 0 0];
K_ob=acker(A',C',Pole)';

N_sim=61;
x=zeros(Nc,N_sim);
y=zeros(1,N_sim);
Delta_U=zeros(1,N_sim);
u=zeros(1,N_sim);

Delta_u1=zeros(1,N_sim);
x1=zeros(Nc,N_sim);
y1=zeros(1,N_sim);
u1=zeros(1,N_sim);
for kk=1:N_sim-1
    Delta_U1=inv(P_P+0.01*eye(Nc))*(P_Rs-P_F*x1(:,kk));
    Delta_u1(:,kk)=Delta_U1(1,1);
    x1(:,kk+1)=A*x1(:,kk)+B*Delta_u1(:,kk)+K_ob*(y1(kk)-C*x1(:,kk));
    y1(kk+1)=C*x1(:,kk+1);
    u_tem=u1(kk)+Delta_u1(kk);
    if ((u_tem>=-3)&&(u_tem<=6))
        u1(kk+1)=u_tem;
    else
        f=-2*(P_Rs-P_F*x1(:,kk));
        b=[6-u1(kk);6-u1(kk);6-u1(kk);3+u1(kk);3+u1(kk);3+u1(kk)];
        eta=QPhild(H,f,A_cons,b);
        Delta_u1(:,kk)=eta(1,1);
        u1(kk+1)=Delta_u1(:,kk)+u1(kk);
    end
        
end

t=0:1:N_sim-1;
figure
subplot(211)
plot(t,y1)
xlabel('Sampling Instant')
ylabel('Output y')
subplot(212)
stairs(t,u1)
xlabel('Sampling Instant')
ylabel('Control signal u')

figure
plot(t,Delta_u1)
xlabel('Sampling Instant')
ylabel('Delta u')


