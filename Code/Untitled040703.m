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
[Omega,Psi]=dmpc(A,B,a,N,Np,Q,R);
Deltau_min=-1;
Deltau_max=0.25;
Nc=15;
X0=[0.1;0.2;0.3];
[M0,Lzerot]=Mdu(a,N,1,Nc);
M=[M0;-M0];
gamma=[Deltau_max*ones(Nc,1);-Deltau_min*ones(Nc,1)];
eta=QPhild(Omega,Psi*X0,M,gamma);

[Al,L0]=lagd(a,N);
L(:,1)=L0;
for kk=2:Np;
    L(:,kk)=Al*L(:,kk-1);   % total L sequence
end
 
delta_U=zeros(1,Np);
for i=1:N
    delta_U=delta_U+eta(i,1)*L(i,:);
end

x=X0;
for i=1:Np
    x(:,i+1)=A*x(:,i)+B*delta_U(i);
    y(:,i)=C*x(:,i);
end

t=0:1:Np-1;
figure
stairs(t,delta_U);
xlabel('Sampling Instant')
ylabel('Optimal Control')
line([0,45],[-1,-1],'linestyle',':')
line([0,45],[0.25,0.25],'linestyle',':')
axis([0 46 -1.5 0.5])
figure
plot(t,y)
xlabel('Sampling Instant')
ylabel('Optimal Output')
axis([0 46 -0.5 1])
