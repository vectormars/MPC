clc;clear;
p=1;
N=4;
delta_t=0.01;
Tm=8;            % simulation time
N_sample=Tm/delta_t;  % sample number
t=0:delta_t:(N_sample-1)*delta_t;
[Ap,L0]=lagc(p,N);
L=zeros(N,N_sample);
for i=1:N_sample;
    L(:,i)=expm(Ap*t(i))*L0;
end

num=1;
den=conv([1 1],[1 2]);
h_t=impulse(num,den,t)';  % impluse response
hL_1=L(1,:).*h_t;    
c1=INTE(hL_1,delta_t);
hL_2=L(2,:).*h_t;
c2=INTE(hL_2,delta_t);
hL_3=L(3,:).*h_t;
c3=INTE(hL_3,delta_t);
hL_4=L(4,:).*h_t;
c4=INTE(hL_4,delta_t);
%%%%%%%%%  Print  %%%%%%%%%%%%%%%%%%%%%%%
h_model2=c1*L(1,:)+c2*L(2,:);
h_model4=c1*L(1,:)+c2*L(2,:)+c3*L(3,:)+c4*L(4,:);
figure 
subplot(211)
plot(t,h_t,t,h_model2,'--');
legend('impulse','c1,c2')
xlabel('Time(sec)')
subplot(212)
plot(t,h_t,t,h_model4,'--');
legend('impulse','c1,c2,c3,c4')
xlabel('Time (sec)')

