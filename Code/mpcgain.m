function [Phi,F]=mpcgain(A,B,C,Nc,Np);
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