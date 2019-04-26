clear;clc;
Am=[0 1;-4 0];
Bm=[1;0];
Cm=[0 1];
Dm=0;
delta_t=0.1;
Np=10;
Nc=3;
[Ad,Bd,Cd,Dd]=c2dm(Am,Bm,Cm,Dm,delta_t);
 
[m1,n1]=size(Cd);
[n1,n_in]=size(Bd);
A_e=eye(n1+m1,n1+m1);
A_e(1:n1,1:n1)=Ad;
A_e(n1+1:n1+m1,1:n1)=Cd*Ad;
B_e=zeros(n1+m1,n_in);
B_e(1:n1,:)=Bd;
B_e(n1+1:n1+m1,:)=Cd*Bd;
C_e=zeros(m1,n1+m1);
C_e(:,n1+1:n1+m1)=eye(m1,m1);
 
A=A_e;
B=B_e;
C=C_e;
 
Phi=zeros(Np,Nc);
 
for i=1:Np
    F(i,:)=C*A^(i);
end
 
 
for i=1:Np
    Phi(i,1)=C*(A^(i-1))*B;
end
 
v=Phi(:,1);
 
for i=2:Nc
    Phi(:,i)=[zeros(i-1,1);v(1:Np-i+1,1)];
end
 
K_mpc=inv(Phi'*Phi)*(Phi'*F);
K_mpc=K_mpc(1,:);
Rs=ones(Np,1);
K_y=inv(Phi'*Phi)*(Phi'*Rs);
K_y=K_y(1,1);

Pole=[0 0 0];
K_ob=acker(A',C',Pole)';

N=15;
delta_u=0;
x=zeros(3,N+1);
delta_U=zeros(1,N+1);
y=zeros(1,N+1);
r=ones(1,N+1);
u=zeros(1,N+1);


for kk=1:N
    x(1:3,kk+1)=A*x(1:3,kk)+B*delta_U(kk)+K_ob*(y(kk)-C*x(1:3,kk));
    delta_U(kk+1)=K_y*r(kk+1)-K_mpc*x(1:3,kk+1);
    u(kk+1)=u(kk)+delta_U(kk+1);
    if -25<=u(kk+1)&&u(kk+1)<=25;
        u(kk+1)=u(kk+1);
    elseif u(kk+1)>25
        u(kk+1)=25;
        delta_U(kk+1)=25-u(kk);
    else u(kk+1)<-25
        u(kk+1)=-25
        delta_U(kk+1)=-25-u(kk);
    end
    y(kk+1)=C*x(1:3,kk+1);  
end


subplot(211)
y1=y(2:N+1);
t=0:1:N-1;
plot(t,y1)
axis([0 8 0 1.5])
xlabel('Sampling Instant')
ylabel('Output y')
subplot(212)
u1=u(1:N);
stairs(t,u1)
axis([0 8 -40 40])
xlabel('Sampling Instant')
ylabel('Input u')




