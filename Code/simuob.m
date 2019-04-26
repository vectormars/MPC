function [u1,y1,deltau1,k]=simuob(xm,u,y,sp,Ap,Bp,Cp,A,B,C,N_sim,Omega,Psi,K_ob,Lzerot)
%closed-loop simulation without constraints
[ny,n1]=size(Cp);
[n1,nu]=size(Bp);
X_hat=zeros(n1,1);
for kk=1:N_sim;
Xsp=[zeros(n1-ny,1);sp(:,kk)];
eta=-(Omega\Psi)*(X_hat-Xsp);
deltau=Lzerot*eta;
u=u+deltau; %update u
deltau1(:,kk)=deltau;
u1(1:nu,kk)=u; %keep u
y1(1:ny,kk)=y; %keep y
X_hat=A*X_hat+K_ob*(y-C*X_hat)+B*deltau;
%u and y to generate X_hat(k+1)
%%%%
%plant simulation
%%%%%%
xm=Ap*xm+Bp*u; % calculate xm(k+1)
y=Cp*xm; %calculate y(k+1)
end
k=0:(N_sim-1);