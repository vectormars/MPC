clc;clear;
num=1;
den=[1 2*0.5*3 3*3];
[Ac,Bc,Cc,Dc]=tf2ss(num,den);
Delta_t=0.1;
[Ap,Bp,Cp,Dp]=c2dm(Ac,Bc,Cc,Dc,Delta_t);
[A,B,C]=format_change(Ap,Bp,Cp);

n1=2;
m1=1;
n_in=1;

Q=C'*C; 
R=0.3;
Np=46;
a=0.7;
N=18;
alpha=1.2;
[K,P,Z]=dlqr(A,B,Q,R);
N_sim=50;

%%%%%%%%%%% Without constraints %%%%%%%%%%%%%%%%%%
gamma=1/alpha;
g=gamma*gamma;
Qa=g*Q+(1-g)*P;
Ra=g*R;
[Omega,Psi]=dmpc(A/alpha,B/alpha,a,N,Np,Qa,Ra);
[Al,L0]=lagd(a,N);

xm=zeros(n1,1);
y=0;
x1_delta=0;
x2_delta=0;
r=ones(1,120);
Xf=[x1_delta;x2_delta;y-r(1,1)];
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
    u_delta_k_m1=deltau;
    Xf=[xm-xm_old;y-r(1,kk+1)];
    up=u;
end

t=0:1:N_sim-1;
figure(1)
plot(t,deltau2)
axis([0 40 -0.1 2])
figure(2)
plot(t,u2)
axis([0 40 0 10])
figure(3)
plot(t,y2)
axis([0 40 -0.1 1.1])


