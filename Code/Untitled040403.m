clc;clear;
A=[0.8 0;0.8 1];
B=[0.6;0.6];
C=[0 1];
x0=[0.1;0.2];
Q=C'*C;

R=0.1;
a=0.6;N=5;
Np=16;

[mA,mB]=size(A);
[Al,L0]=lagd(a,N);
Sc_1=B*L0';
[m,n]=size(Sc_1);   %mA by N
Sc(:,1:n)=Sc_1;   
for i=1:Np-1
    Sc(:,n*i+1:n*i+n)=A*Sc(:,n*i-n+1:n*i)+Sc_1*(Al^(i))';
end
 
for i=1:Np
    Omega_1(:,n*i-n+1:n*i)=Sc(:,n*i-n+1:n*i)'*Q*Sc(:,n*i-n+1:n*i);
end
 
E=zeros(N,N);
for i=1:Np
    E(:,:)=E(:,:)+Omega_1(:,n*i-n+1:n*i);
end
E=E+R*eye(N,N);
 
 
for i=1:Np
    Cha_0(:,mA*i-1:mA*i)=Sc(:,n*i-n+1:n*i)'*Q*(A^(i));
end
 
H(:,:)=zeros(N,mA);
for i=1:Np
    H(:,:)=H(:,:)+Cha_0(:,mA*i-1:mA*i);
end
 
eta=-(E\H)*x0;

%%%%%%%%%%%%%%% Laguerre N=5 %%%%%%%%%%%%%%%%%%%%
L(:,1)=L0;
for kk=2:Np;
    L(:,kk)=Al*L(:,kk-1);   % total L sequence
end
 
delta_U=zeros(1,Np);
for i=1:N
    delta_U=delta_U+eta(i,1)*L(i,:);
end

x=x0;
for i=1:Np
    x(:,i+1)=A*x(:,i)+B*delta_U(i);
    y(:,i)=C*x(:,i);
end
%%%%%%%%%%%% DLQR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[K,P,Z]=dlqr(A,B,Q,R);  % DLQR
for i=1:Np
    delta_U1(i)=-K*(A-B*K)^(i-1)*x0;
end

x1=x0;
for i=1:Np
    x1(:,i+1)=A*x1(:,i)+B*delta_U1(i);
    y1(:,i)=C*x1(:,i);
end

%%%%%%%%%%%% Print %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=0:1:Np-1;
figure
stairs(t,delta_U);hold on
stairs(t,delta_U1,'r')
xlabel('Sampling Instant')
ylabel('Optimal Control')
legend('Laguerre (N=5)','DLQR')

figure
plot(t,y);hold on
plot(t,y1,'r')
xlabel('Sampling Instant')
ylabel('Optimal Output')
legend('Laguerre (N=5)','DLQR')
