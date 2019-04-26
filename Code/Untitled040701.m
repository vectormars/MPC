clc;clear;
num=1;
den=[1 2*0.1*3 3*3];
[Ac,Bc,Cc,Dc]=tf2ss(num,den);
Delta_t=0.1;
[Ad,Bd,Cd,Dd]=c2dm(Ac,Bc,Cc,Dc,Delta_t);
[A,B,C]=format_change(Ad,Bd,Cd);
N=8;Np=46;Nc=15;
a=0.7;
R=0.3;
Q=C'*C;

%%%%%%%%%%%%%%% with constraint on U %%%%%%%%%%%%%%%%%%%%
[Omega,Psi]=dmpc(A,B,a,N,Np,Q,R);
X0=[0.1;0.2;0.3];
up=6;
u_min=1.8;
u_max=4;
M0=Mu(a,N,1,Np);
M=[M0;-M0];
gamma=[(u_max-up)*ones(Np,1);(-u_min+up)*ones(Np,1)];
eta=QPhild(Omega,Psi*X0,M,gamma);
[Al,L0]=lagd(a,N);
L(:,1)=L0;

for kk=2:Np;
    L(:,kk)=Al*L(:,kk-1);   % total L sequence
end
 
delta_U=zeros(1,Np);
for i=1:N
    delta_U=delta_U+eta(i,:)*L(i,:);
end

x=X0;
for i=1:Np
    x(:,i+1)=A*x(:,i)+B*delta_U(i);
    y(:,i)=C*x(:,i);
end

U=zeros(1,Np);
U(1,1)=up;
for i=1:Np
    U(i+1)=U(i)+delta_U(i);
end
U=U(2:Np+1);

%%%%%%%%%%%% without constraint %%%%%%%%%%%%%%%%%%%%%%%%
[Omega,Psi]=dmpc(A,B,a,N,Np,Q,R);
eta1=-(Omega\Psi)*X0;
[Al,L0]=lagd(a,N);
L(:,1)=L0;
for kk=2:Np;
    L(:,kk)=Al*L(:,kk-1);
end

delta_U1=zeros(1,Np);
for i=1:N
    delta_U1=delta_U1+eta1(i,1)*L(i,:);
end

x1=X0;
for i=1:Np
    x1(:,i+1)=A*x1(:,i)+B*delta_U1(i);
    y1(:,i)=C*x1(:,i);
end

U1=zeros(1,Np+1);
U1(1,1)=up;
for i=1:Np
    U1(i+1)=U1(i)+delta_U1(i);
end
U1=U1(2:Np+1);

%%%%%%%%%%%% Print %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=0:1:Np-1;
figure
stairs(t,delta_U);hold on
stairs(t,delta_U1,'r')
xlabel('Sampling Instant')
ylabel('Delta u')
legend('with constraint on u','without constraint')
figure
stairs(t,U);hold on
stairs(t,U1,'r')
xlabel('Sampling Instant')
ylabel('u')
legend('with constraint on u','without constraint')
axis([0 45 1 5])
line([0,45],[4,4],'linestyle',':')
line([0,45],[1.8,1.8],'linestyle',':')