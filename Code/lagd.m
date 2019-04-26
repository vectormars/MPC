function [Al,L0]=lagd(a,N)
b=sqrt(1-a*a);
L0=zeros(N,1);
for i=1:N;
    L0(i,1)=(-1)^(i-1)*a^(i-1);
    L0(i,1)=L0(i,1)*b;
end

Al=zeros(N,N);
A1=zeros(N,1);
A1(1,1)=a;
for i=2:N;
    A1(i,1)=(-a)^(i-2)*(1-a*a);
end
Al(1:N,1)=A1;
for i=2:N;
    Al(1:N,i)=[zeros(i-1,1);A1(1:N-i+1,1)];
end
