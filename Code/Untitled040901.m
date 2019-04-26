clc;clear;
num=1;
den=[1 2*0.5*3 3*3];
[Ac,Bc,Cc,Dc]=tf2ss(num,den);
Delta_t=0.1;
[Ap,Bp,Cp,Dp]=c2dm(Ac,Bc,Cc,Dc,Delta_t);
[A,B,C]=format_change(Ap,Bp,Cp);
Q=C'*C;

R=0.3;
Np=46;
a=0.7;
N=18;
X0=[0.1 0.2 0.3]';
alpha=1.2;
[K,P,Z]=dlqr(A,B,Q,R);

gamma=1/alpha;
g=gamma*gamma;
Qa=g*Q+(1-g)*P;
Ra=g*R;
[Omega,Psi]=dmpc(A/alpha,B/alpha,a,N,Np,Qa,Ra);
[Al,L0]=lagd(a,N);
Kmpc=L0'*(Omega\Psi);

%%%%%%%%%%%%%%% DMPC %%%%%%%%%%%%%%%%%%%%
x=X0;
delta_U=zeros(1,Np);
y=zeros(1,Np);
U=zeros(1,Np);
for i=1:Np
    x(:,i+1)=(A-B*Kmpc)*x(:,i);
    delta_U(i)=-Kmpc*x(:,i);
end

for i=1:Np
    U(i)=-K*((A/alpha)-(B/alpha)*K)^(i-1)*X0;
end
%%%%%%%%%%%% DLQR %%%%%%%%%%%%%%%%%%%%%%%
delta_U1=zeros(1,Np);
y1=zeros(1,Np);
U1=zeros(1,Np);

for i=1:Np
    delta_U1(i)=-K*(A-B*K)^(i-1)*X0;
end

x1=X0;
for i=1:Np
    x1(:,i+1)=A*x1(:,i)+B*delta_U1(i);
end

for i=1:Np
    U1(i+1)=U(i)+delta_U(i);
end
%%%%%%%%%%%% Print %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=0:1:Np-1;
figure(1)
stairs(t,delta_U);hold on
stairs(t,delta_U1,'r')
xlabel('Sampling Instant')
ylabel('Optimal Control')
legend('DMPC','DLQR')
axis([0 20 -2.5 0.5])
figure(2)
%stairs(t,x(3,:));hold on
%stairs(t,x1(3,:),'r')
xlabel('Sampling Instant')
ylabel('Optimal Control')
%legend('DMPC','DLQR')
axis([0 20 -2.5 0.5])
figure(3)
x=x(:,1:Np);
x1=x1(:,1:Np);
plot(t,x(3,:));hold on
plot(t,x1(3,:),'r')
xlabel('Sampling Instant')
ylabel('Optimal Output')
legend('DMPC','DLQR')
axis([0 20 -0.1 1])




