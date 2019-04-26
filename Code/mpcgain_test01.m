clc;clear;
Am=0.8;
Bm=0.1;
Cm=1;
Np=10;Nc=4;
[Phi_Phi,Phi_F,Phi_R,A,B,C]=mpcgain(Am,Bm,Cm,Nc,Np);
x0=[0.1;0.2];
Delta_U1=inv(Phi_Phi)*(Phi_R-Phi_F*x0);
rw=10;
Delta_U2=inv(Phi_Phi+rw*eye(4,4))*(Phi_R-Phi_F*x0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h(1,:)=C;
F(1,:)=C*A;
for kk=2:Np
    h(kk,:)=h(kk-1,:)*A;
    F(kk,:)= F(kk-1,:)*A;
end
v=h*B;
Phi=zeros(Np,Nc); %declare the dimension of Phi
Phi(:,1)=v; % first column of Phi
for i=2:Nc
    Phi(:,i)=[zeros(i-1,1);v(1:Np-i+1,1)]; %Toeplitz matrix
end

% Control I. with rw=0;
N=10;
delta_u1=zeros(Nc,N+1);
x=zeros(2,N+1);
x(:,1)=x0;
y=zeros(1,N);
u=zeros(1,N+1);

for i=1:N
    Delta_U1=inv(Phi_Phi)*(Phi_R-Phi_F*x(:,i));
    y(i)=C*x(:,i);
    x(:,i+1)=A*x(:,i)+B*Delta_U1(1);
end
t=1:N;
figure(1)
plot(t,y,'-.',t,x(1,1:end-1));
xlabel('Sample instant')
ylabel('Response')
legend('y','\Deltax_m')
savefig('rw_0.fig')
print('-dpng','rw_0.png')

% Control II. with rw=10;
for i=1:N
    Delta_U2=inv(Phi_Phi+rw*eye(4,4))*(Phi_R-Phi_F*x0);
    y(i)=C*x(:,i);
    x(:,i+1)=A*x(:,i)+B*Delta_U2(1);
end
t=1:N;
figure(2)
plot(t,y,'-.',t,x(1,1:end-1));
xlabel('Sample instant')
ylabel('Response')
legend('y','\Deltax_m',2)
savefig('rw_10.fig')
print('-dpng','rw_10.png')