clc;clear
K=1; % gain of the motor
T=1; % time constant of the motor
Delta_t=0.1; % sampling interval
num=K;
den=conv([1 0],[T 1]);
[numd,dend]=c2dm(num,den,Delta_t);
Ap=zeros(3,3);
Ap(1,:)=[-dend(1,2:3) numd(1,3)];Ap(2,:)=[1 0 0];Ap(3,:)=[0 0 0];
Bp(1,1)=numd(1,2);Bp(2,1)=0;Bp(3,1)=1;
Cp(1,1)=1;Cp(1,2)=0;Cp(1,3)=0;
Dp(1,1)=0;
[A,B,C]=format_change(Ap,Bp,Cp);

n1=3;
m1=1;
n_in=1;

Q=C'*C;
R=0.1;
a=0.7;
N=10; %you can increase N
Np=46;
N_sim=250;
[Omega,Psi]=dmpc(A,B,a,N,Np,Q,R);
[A1,L0]=lagd(a,N);
%%%%%%%%%%% With constraints on y  %%%%%%%%%
xm=zeros(n1,1);
y=0;
y_delta_k=0;
y_delta_k_m1=0;
u_delta_k_m1=0;
u=0;
up=0;
deltau1=zeros(1,N_sim);
u1=zeros(1,N_sim);
y1=zeros(1,N_sim);
r=[ones(1,120) -zeros(1,120) ones(1,200)]*0;
d=[ones(1,120) -zeros(1,120) ones(1,200)];
Xf=[y_delta_k; y_delta_k_m1;u_delta_k_m1;y];
u_min=-1.5;
u_max=1.5;
deltau_min=-0.4;
deltau_max=0.4;
y_min=-0.14;
y_max=0.14;
M_act=C*B*L0';
E=C*A;

for kk=1:N_sim;
    eta=-(Omega\Psi)*Xf;
    deltau=L0'*eta;
    beta=E*Xf+M_act*eta;
    if (beta<y_min)
        M_act1=-M_act;
        lambda=-(M_act1*((Omega)\M_act1'))\(-y_min+E*Xf+M_act1*((Omega)\Psi)*Xf);
        eta=-(Omega)\(Psi*Xf+M_act1'*lambda);
        deltau=L0'*eta;
    end
    if (beta>y_max)
        lambda=-(M_act*((Omega)\M_act'))\(y_max-E*Xf+M_act*((Omega)\Psi)*Xf);
        eta=-(Omega)\(Psi*Xf+M_act'*lambda);
        deltau=L0'*eta; 
    end
    if(deltau>deltau_max) 
        deltau=deltau_max;
    end
    if (deltau<deltau_min)
        deltau=deltau_min;
    end
    u=up+deltau;
    if (u>u_max) 
        deltau=u_max-up; 
        u=u_max;
    end
    if (u<u_min) 
        deltau=u_min-up;
        u=u_min; 
    end 
    deltau1(1,kk)=deltau;
    u1(1,kk)=u;
    y1(1,kk)=y;
%%%%
%plant simulation
%%%%%%
yp=y;
xm_old=xm;
xm=Ap*xm+Bp*(u+d(kk)); % calculate xm(k+1)
y=Cp*xm; %calculate y(k+1)
%updating feedback state variable Xf
y_delta_k_m1=y_delta_k;
y_delta_k=y-yp;
u_delta_k_m1=deltau;
Xf=[y_delta_k; y_delta_k_m1;u_delta_k_m1;y];
up=u;
end

%%%%%%%% Print %%%%%%%%%%%%%%%%%%%%%%%%%%%
k=0:(N_sim-1);
figure(1)
subplot(211)
plot(k,y1);hold on
%plot(k,y2,'r')
ylabel('y')
xlabel('Sampling Instant')
axis([0 N_sim-1 -0.5 0.5])
line([0,N_sim],[0.14,0.14],'linestyle',':')
line([0,N_sim],[-0.14,-0.14],'linestyle',':')
subplot(212)
plot(k,u1);hold on
%plot(k,u2,'r')
ylabel('u')
xlabel('Sampling Instant')
line([0,N_sim],[1.5,1.5],'linestyle',':')
line([0,N_sim],[-1.5,-1.5],'linestyle',':')
axis([0 N_sim-1 -2 2])
figure(2)
plot(k,deltau1);hold on
%plot(k,deltau2,'r')
ylabel('Delta u')
xlabel('Sampling Instant')
line([0,N_sim],[-0.4,-0.4],'linestyle',':')
line([0,N_sim],[0.4,0.4],'linestyle',':')
axis([0 N_sim-1 -0.51 0.51])
