clc;clear;
num=[1 -3];
den=[1 2*0.3*3 9];
[Am Bm Cm Dm]=tf2ss(num,den);
delta_t=0.1;
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

Np=36;
Q=C'*C;
R=0.3;

a=0.4;
N=8;
 
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