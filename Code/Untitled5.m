clc;clear;
num=[-5.7980 19.5128 -21.6452 7.9547];
den=[1 -3.0228 3.8630 -2.6426 0.8084];
[Ap,Bp,Cp,Dp]=tf2ss(num,den);
[A,B,C]=format_change(Ap,Bp,Cp);

Q=C'*C;
R=0.3;

X0=[0.1 0.2 0.3 0.4 0.5]';
Np=46;N=10;a=0.4;
alpha=1.2;
[K,P,Z]=dlqr(A,B,Q,R);
N_sim=60;
%%% closed-loop predictive control with alpha=1.2 without constraint %%%
gamma=1/alpha;
g=gamma*gamma;
Qa=g*Q+(1-g)*P;
Ra=g*R;
[Omega,Psi]=dmpc(A/alpha,B/alpha,a,N,Np,Qa,Ra);
[Al,L0]=lagd(a,N);
n1=4;
xm=zeros(n1,1);
y=0;
x1_delta=0;
x2_delta=0;
x3_delta=0;
x4_delta=0;
r=ones(1,120);
Xf=[x1_delta;x2_delta;x3_delta;x4_delta;y-r(1,1)];
up=0;
deltau2=zeros(1,N_sim);
u2=zeros(1,N_sim);
y2=zeros(1,N_sim);
for kk=1:N_sim;
    eta=-(Omega\Psi)*Xf;
    deltau=L0'*eta;
    u=up+deltau;
    deltau2(1,kk)=deltau;    %save data
    u2(1,kk)=u;              %save data
    y2(1,kk)=y;              %save data
    %%%%
    %plant simulation
    %%%%%%
    xm_old=xm;
    xm=Ap*xm+Bp*u; % calculate xm(k+1)
    y=Cp*xm; %calculate y(k+1)
    %updating feedback state variable Xf
    Xf=[xm-xm_old;y-r(1,kk+1)];
    up=u;
end
%%% closed-loop predictive control with alpha=1.2 with constraint %%%
xm=zeros(n1,1);
y=0;
x1_delta=0;
x2_delta=0;
x3_delta=0;
x4_delta=0;
r=ones(1,120);
Xf=[x1_delta;x2_delta;x3_delta;x4_delta;y-r(1,1)];
up=0;

u_min=0.1;
u_max=0.4;
deltau_min=-0.01;
deltau_max=0.05;

deltau1=zeros(1,N_sim);
u1=zeros(1,N_sim);
y1=zeros(1,N_sim);
for kk=1:N_sim;
    eta=-(Omega\Psi)*Xf;
    deltau=L0'*eta;
    if (deltau>deltau_max)
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
    deltau1(1,kk)=deltau;    %save data
    u1(1,kk)=u;              %save data
    y1(1,kk)=y;              %save data
    %%%%
    %plant simulation
    %%%%%%
    xm_old=xm;
    xm=Ap*xm+Bp*u; % calculate xm(k+1)
    y=Cp*xm; %calculate y(k+1)
    %updating feedback state variable Xf
    Xf=[xm-xm_old;y-r(1,kk+1)];
    up=u;
end
%%%%%%%%%%%% Print %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=0:1:N_sim-1;
figure(1)
plot(t,deltau2);hold on
plot(t,deltau1,'r')
xlabel('Sampling Instant')
ylabel('Delta u')
legend('without constraint','constraint on the first sample')
line([0,N_sim],[-0.01,-0.01],'linestyle',':')
line([0,N_sim],[0.05,0.05],'linestyle',':')
axis([0 N_sim-1 -0.1 0.15])
figure(2)
plot(t,u2);hold on
plot(t,u1,'r')
xlabel('Sampling Instant')
ylabel('u')
legend('without constraint','constraint on the first sample')
line([0,N_sim],[0.1,0.1],'linestyle',':')
line([0,N_sim],[0.4,0.4],'linestyle',':')
axis([0 N_sim-1 -0.1 0.5])
figure(3)
plot(t,y2);hold on
plot(t,y1,'r')
xlabel('Sampling Instant')
ylabel('y')
legend('without constraint','constraint on the first sample')
axis([0 N_sim-1 -1 1.5])