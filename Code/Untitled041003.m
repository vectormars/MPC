clc;clear;
delta_t=0.01;
Tm=48;
N_sample=Tm/delta_t;  
t=0:delta_t:(N_sample-1)*delta_t;

num=[-1 1];
den=conv([1 1.1],conv([1 0.7],[1 0.1]));
h_t=impulse(num,den,t)'; 

a1=1.1;a2=0.7;a3=0.1;
Ak=[-a1 0 0;-2*a2 -a2 0;-2*a1 -a2 -a3];
Bk=ones(3,1);
Ck=[sqrt(2*a1) 0 0;0 sqrt(2*a2) 0;0 0 sqrt(2*a3)];