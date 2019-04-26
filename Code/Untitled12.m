clc;clear;
A=[0.8 0;0.8 1];
B=[0.6;0.6];
Q=[0 0;0 1];
R=1;
X=[0.1 0.2]';
Np=16;
 
a=0.6;
N=5;
 
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
E=E+eye(N,N);
 
 
for i=1:Np
    Cha_0(:,mA*i-1:mA*i)=Sc(:,n*i-n+1:n*i)'*Q*(A^(i));
end
 
H(:,:)=zeros(N,mA);
for i=1:Np
    H(:,:)=H(:,:)+Cha_0(:,mA*i-1:mA*i);
end
 
eta=-(inv(E)*H)*X;
L(:,1)=L0;
for kk=2:Np;
    L(:,kk)=Al*L(:,kk-1);
end
 
delta_U=zeros(size(L(1,:)));
for i=1:N
    delta_U=delta_U+eta(i,1)*L(i,:);
end
 
x=X;
for i=1:Np
    x(:,i+1)=A*x(:,i)+B*delta_U(i);
end

x5=x;
save('x5')






