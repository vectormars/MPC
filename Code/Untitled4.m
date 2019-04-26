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
%%%%%%%%%%%% DLQR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[K,P,Z]=dlqr(A,B,Q,R);
delta_u1=zeros(1,Np);
for i=1:Np
    delta_u1(i)=-(K*(((A/alpha)-(B/alpha)*K)^(i-1)))*X0;
end
 
x1=X0;
y1=zeros(1,Np);
for i=1:Np
    x1(:,i+1)=(A/alpha)*x1(:,i)+(B/alpha)*delta_u1(i);
    y1(:,i)=C*x1(:,i);
end
%%%%%%%%%%%% a=0.4,Np=46 and N=10 %%%%%%%%%%%%%%%%%%%
gamma=1/alpha;
g=gamma*gamma;
Qa=g*Q+(1-g)*P;
Ra=g*R;
[Omega,Psi]=dmpc(A/alpha,B/alpha,a,N,Np,Qa,Ra);
[Al,L0]=lagd(a,N);

eta=-Omega\Psi*X0;
L(:,1)=L0;
for i=2:Np;
    L(:,i)=Al*L(:,i-1);
end
delta_u2=zeros(size(L(1,:)));
for i=1:N
    delta_u2=delta_u2+eta(i,1)*L(i,:);
end
 
x2=X0;
y2=zeros(1,Np);
for i=1:Np
    x2(:,i+1)=(A/alpha)*x2(:,i)+(B/alpha)*delta_u2(i);
    y2(:,i)=C*x2(:,i);
end
%%%%%%%%%%%% Print %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=0:1:Np-1;
delta_u1=delta_u1(1:46);
y1=y1(1:46);
figure(1)
stairs(t,delta_u1);hold on
stairs(t,delta_u2,'r');
axis([0 40 -0.1 0.1])
xlabel('Sampling Instant')
ylabel('Delta u')
legend('LQR','a=0.4,N=10,alpha=1.2')
figure(2)
plot(t,y1);hold on
plot(t,y2,'r');
xlabel('Sampling Instant')
ylabel('Optimal Output')
axis([0 40 -0.5 1])
legend('LQR','a=0.4,N=10,alpha=1.2')

