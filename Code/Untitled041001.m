clc;clear
p=10;
N=3;
[Ap,L0]=lagc(p,N);
delta_t=0.01;
Tm=8;
N_sample=Tm/delta_t;
t=0:delta_t:(N_sample-1)*delta_t;
%solution of the differential equation
for i=1:N_sample;
    L(:,i)=expm(Ap*t(i))*L0;
end
figure
plot(t,L(1,:),t,L(2,:),'--',t,L(3,:),'r')
legend('L1','L2','L3')
xlabel('Time (sec)')