clc;clear;
Ac=[0 1;-4 0];
Bc=[0;1];
Cc=[0 1];
Dc=0;
Delta_t=0.1;
[Ad,Bd,Cd,Dd]=c2dm(Ac,Bc,Cc,Dc,Delta_t);

N=60;
x0=[1;0];
x0_hat=[0.3;0];
X=zeros(2,N);
X(:,1)=x0;
for i=1:60;
    X(:,i+1)=Ad*X(:,i);
end
X1=X(1,1:end-1);
X2=X(2,1:end-1);

% Without observor
XI=zeros(2,N);
XI(:,1)=x0_hat;
for i=1:60;
    XI(:,i+1)=Ad*XI(:,i);
end
XI1=XI(1,1:end-1);
XI2=XI(2,1:end-1);  

t=1:N;
figure(1)
plot(t,X1,t,XI1,'--')
legend('x_1','xhat_1')
savefig('Est_no_ob_1.fig')
print('-dpng','Est_no_ob_1.png')

figure(2)
plot(t,X2,t,XI2,'--')
legend('x_2','xhat_2',2)
savefig('Est_no_ob_2.fig')
print('-dpng','Est_no_ob_2.png')

% With observor
Pole=[0.1 0.2];
K_ob=place(Ad',Cd',Pole)';
XII=zeros(2,N);
XII(:,1)=x0_hat;
for i=1:60;
    X(:,i+1)=Ad*X(:,i);
    XII(:,i+1)=Ad*XII(:,i)+K_ob*(X(2,i)-XII(2,i));
end
XII1=XII(1,1:end-1);
XII2=XII(2,1:end-1);  

figure(3)
plot(t,X1,t,XII1,'--')
legend('x_1','xhat_1')
savefig('Est_ob_1.fig')
print('-dpng','Est_ob_1.png')

figure(4)
plot(t,X2,t,XII2,'--')
legend('x_2','xhat_2',2)
savefig('Est_ob_2.fig')
print('-dpng','Est_ob_2.png')
