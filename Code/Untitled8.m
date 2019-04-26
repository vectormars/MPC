clc;clear;
Ts=1;
G=tf(0.1,[1 -0.8],Ts,'Inputdelay',6);
Gsmin=ss(G);
[Ad,Bd,Cd,Dd]=ssdata(Gsmin);
[A,B,C]=format_change(Ad,Bd,Cd);
Q=C'*C;
R=0.1;
alpha=1.2;
N=1;
Np=46;
[K,P,Z]=dlqr(A,B,Q,R);
N_sim=50;
n1=1;
gamma=1/alpha;
g=gamma*gamma;
Qa=g*Q+(1-g)*P;
Ra=g*R;
%%%%%%%%%%% a=0 %%%%%%%%%%%%
a=0;
[Omega,Psi]=dmpc(A/alpha,B/alpha,a,N,Np,Qa,Ra);
[Al,L0]=lagd(a,N);
xm=zeros(n1,1);
y=0;
x_delta=0;
r=ones(1,120);
Xf=[x_delta;y-r(1,1)];
up=0;
deltau1=zeros(1,N_sim);
u1=zeros(1,N_sim);
y1=zeros(1,N_sim);
for kk=1:N_sim;
    eta=-(Omega\Psi)*Xf;
    deltau=L0'*eta;
    u=up+deltau;
    deltau1(1,kk)=deltau;    %save data
    u1(1,kk)=u;              %save data
    y1(1,kk)=y;              %save data
    %%%%
    %plant simulation
    %%%%%%
    xm_old=xm;
    xm=Ad*xm+Bd*u; % calculate xm(k+1)
    y=Cd*xm; %calculate y(k+1)
    %updating feedback state variable Xf
    Xf=[xm-xm_old;y-r(1,kk+1)];
    up=u;
end
y1=[zeros(1,6) y1];
y1=y1(1:50);
%%%%%%%%%%% a=0.3 %%%%%%%%%%%%
a=0.3;
[Omega,Psi]=dmpc(A/alpha,B/alpha,a,N,Np,Qa,Ra);
[Al,L0]=lagd(a,N);
xm=zeros(n1,1);
y=0;
x_delta=0;
r=ones(1,120);
Xf=[x_delta;y-r(1,1)];
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
    xm=Ad*xm+Bd*u; % calculate xm(k+1)
    y=Cd*xm; %calculate y(k+1)
    %updating feedback state variable Xf
    Xf=[xm-xm_old;y-r(1,kk+1)];
    up=u;
end
y2=[zeros(1,6) y2];
y2=y2(1:50);
%%%%%%%%%%% a=0.6 %%%%%%%%%%%%
a=0.6;
[Omega,Psi]=dmpc(A/alpha,B/alpha,a,N,Np,Qa,Ra);
[Al,L0]=lagd(a,N);
xm=zeros(n1,1);
y=0;
x_delta=0;
r=ones(1,120);
Xf=[x_delta;y-r(1,1)];
up=0;
deltau3=zeros(1,N_sim);
u3=zeros(1,N_sim);
y3=zeros(1,N_sim);
for kk=1:N_sim;
    eta=-(Omega\Psi)*Xf;
    deltau=L0'*eta;
    u=up+deltau;
    deltau3(1,kk)=deltau;    %save data
    u3(1,kk)=u;              %save data
    y3(1,kk)=y;              %save data
    %%%%
    %plant simulation
    %%%%%%
    xm_old=xm;
    xm=Ad*xm+Bd*u; % calculate xm(k+1)
    y=Cd*xm; %calculate y(k+1)
    %updating feedback state variable Xf
    Xf=[xm-xm_old;y-r(1,kk+1)];
    up=u;
end
y3=[zeros(1,6) y3];
y3=y3(1:50);
%%%%%%%%%%% a=0.9 %%%%%%%%%%%%
a=0.9;
[Omega,Psi]=dmpc(A/alpha,B/alpha,a,N,Np,Qa,Ra);
[Al,L0]=lagd(a,N);
xm=zeros(n1,1);
y=0;
x_delta=0;
r=ones(1,120);
Xf=[x_delta;y-r(1,1)];
up=0;
deltau4=zeros(1,N_sim);
u4=zeros(1,N_sim);
y4=zeros(1,N_sim);
for kk=1:N_sim;
    eta=-(Omega\Psi)*Xf;
    deltau=L0'*eta;
    u=up+deltau;
    deltau4(1,kk)=deltau;    %save data
    u4(1,kk)=u;              %save data
    y4(1,kk)=y;              %save data
    %%%%
    %plant simulation
    %%%%%%
    xm_old=xm;
    xm=Ad*xm+Bd*u; % calculate xm(k+1)
    y=Cd*xm; %calculate y(k+1)
    %updating feedback state variable Xf
    Xf=[xm-xm_old;y-r(1,kk+1)];
    up=u;
end
y4=[zeros(1,6) y4];
y4=y4(1:50);
%%%%%%%%%%% closed-loop DLQR %%%%%%%%%%%%
xm=zeros(n1,1);
y=0;
x_delta=0;
r=ones(1,120);
Xf=[x_delta;y-r(1,1)];
up=0;
deltau5=zeros(1,N_sim);
u5=zeros(1,N_sim);
y5=zeros(1,N_sim);
for kk=1:N_sim;
    deltau=-K*Xf;
    u=up+deltau;
    deltau5(1,kk)=deltau;    %save data
    u5(1,kk)=u;              %save data
    y5(1,kk)=y;              %save data
    %%%%
    %plant simulation
    %%%%%%
    xm_old=xm;
    xm=Ad*xm+Bd*u; % calculate xm(k+1)
    y=Cd*xm; %calculate y(k+1)
    %updating feedback state variable Xf
    Xf=[xm-xm_old;y-r(1,kk+1)];
    up=u;
end
y5=[zeros(1,6) y5];
y5=y5(1:50);
%%%%%%%%%%%%%%%%% Print %%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=0:1:N_sim-1;
figure(1)
plot(t,deltau1);hold on
plot(t,deltau2,'r');hold on
plot(t,deltau3,'g');hold on
plot(t,deltau4,'k');hold on
plot(t,deltau5,'--');
axis([-4 40 -1 3])
xlabel('Sampling Instant')
ylabel('Delta u')
legend('a=0','a=0.3','a=0.6','a=0.9','DLQR')
figure(2)
plot(t,u1);hold on
plot(t,u2,'r');hold on
plot(t,u3,'g');hold on
plot(t,u4,'k');hold on
plot(t,u5,'--')
xlabel('Sampling Instant')
ylabel('u')
axis([-4 40 0 4])
legend('a=0','a=0.3','a=0.6','a=0.9','DLQR')
figure(3)
plot(t,y1);hold on
plot(t,y2,'r');hold on
plot(t,y3,'g');hold on
plot(t,y4,'k');hold on
plot(t,y5,'--')
xlabel('Sampling Instant')
ylabel('Output y')
axis([-4 40 -0.1 1.1])
legend('a=0','a=0.3','a=0.6','a=0.9','DLQR')
