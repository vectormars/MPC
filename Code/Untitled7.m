clc;clear;
Ts=0.1;
num=1;
den=conv([1 -0.8],[1 -0.6]);
G=tf(num,den,Ts,'Inputdelay',16);
Gsmin=ss(G);
[Ad,Bd,Cd,Dd]=ssdata(Gsmin);
[A,B,C]=format_change(Ad,Bd,Cd);
Q=C'*C;R=0.1;a=0.5;
N=8;Np=46;
[K,P,Z]=dlqr(A,B,Q,R);
N_sim=50;
n1=2;
alpha=1.2;
gamma=1/alpha;
g=gamma*gamma;
Qa=g*Q+(1-g)*P;
Ra=g*R;

%%%%%%%%%%% a=0.5 %%%%%%%%%%%%
a=0.5;
[Omega,Psi]=dmpc(A/alpha,B/alpha,a,N,Np,Qa,Ra);
[Al,L0]=lagd(a,N);
xm=zeros(n1,1);
y=0;
x1_delta=0;
x2_delta=0;
r=ones(1,120);
Xf=[x1_delta;x2_delta;y-r(1,1)];
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
y1=[zeros(1,16) y1];
y1=y1(1:50);
plot(y1)