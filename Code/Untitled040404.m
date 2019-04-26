clc;clear;
num=1;
den=[1 2*0.1*3 3*3];
[Ac,Bc,Cc,Dc]=tf2ss(num,den);
Delta_t=0.1;
[Ad,Bd,Cd,Dd]=c2dm(Ac,Bc,Cc,Dc,Delta_t);
[A,B,C]=format_change(Ad,Bd,Cd);
N=8;Np=46;
a=0.7;
R=0.3;
Q=C'*C;
X0=[0.1;0.2;0.3];
%%%%%%%%%%%%%%% DMPC with N=5 %%%%%%%%%%%%%%%%%%%%
[Al,L0]=lagd(a,N);
[Omega,Psi]=dmpc(A,B,a,N,Np,Q,R);
Kmpc=L0'*(Omega\Psi);
x=X0;
for i=1:Np
    x(:,i+1)=(A-B*Kmpc)*x(:,i);
    delta_U(i)=-Kmpc*x(:,i);
    y(:,i)=C*x(:,i);
end

%%%%%%%%%%%% DLQR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[K,P,Z]=dlqr(A,B,Q,R);  % DLQR
for i=1:Np
    delta_U1(i)=-K*(A-B*K)^(i-1)*X0;
end

x1=X0;
for i=1:Np
    x1(:,i+1)=A*x1(:,i)+B*delta_U1(i);
    y1(:,i)=C*x1(:,i);
end

%%%%%%%%%%%% Print %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=0:1:Np-1;
figure
plot(t,y);hold on
plot(t,y1,'r')
xlabel('Sampling Instant')
ylabel('Optimal Output')
legend('DMPC','DLQR')
axis([0 46 -0.5 1])
figure
stairs(t,delta_U);hold on
stairs(t,delta_U1,'r')
xlabel('Sampling Instant')
ylabel('Optimal Control')
legend('DMPC','DLQR')
axis([0 46 -1.5 0.5])