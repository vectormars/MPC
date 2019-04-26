% Plant transfer function: g=5.72exp(-14s)/(60s+1)
% Disturbance transfer function: gd=1.52exp(-15s)/ (25s+1)

% Build the step response models for a sampling period of 7.
delt1=0;          % Sampling time.  Can be 0 (for continuous-time system)
delay1=14;        % Pure time delay
num1=5.72;
den1=[60 1];
G = tf(num1,den1,'InputDelay',14);
step(G)

g=poly2tfd(num1,den1,delt1,delay1);   % Create transfer functions
tfinal=245;       % total step=243
delt2=7;          % sampling period=7
nout1=1;          % indicating a stable output
plant=tfd2step(tfinal,delt2,nout1,g);  
% Determines the step response model of a transfer function model
plotstep(plant)


delay2=15;
num2=1.52;
den2=[25 1];
gd=poly2tfd(num2,den2,delt1,delay2);
delt2=7;
nout2=1;
dplant=tfd2step(tfinal,delt2,nout2,gd);

% Calculate the MPC controller gain matrix for   (不考虑控制增量权重!)

model=plant;   % Assume no plant/model mismatch
ywt=1;         % Output Weight=1
uwt=0;         % Input Weight=0
M=5;           % Input Horizon=5
P=20;          % Output Horizon=20

Kmpc1=mpccon(model,ywt,uwt,M,P);  
% Calculate MPC controller gains for unconstrained case

% Simulate and plot response for unmeasured and measured
% step disturbance through dplant.
tend=245;
r=[ ]; usat=[ ]; tfilter=[ ];
dmodel=[ ];   % no measured disturbance
dstep=1;
[y1,u1]=mpcsim(plant,model,Kmpc1,tend,r,usat,tfilter,...
dplant,dmodel,dstep);

dmodel=dplant; % add measured disturbance
[y2,u2]=mpcsim(plant,model,Kmpc1,tend,r,usat,tfilter,...
dplant,dmodel,dstep);
plotall([y1,y2],[u1,u2],delt2);
% Perfect rejection for measured disturbance case.

% Simulate and plot response for unmeasured step
% disturbance through dplant with and without
% input constraints.
% No plant/model mismatch,
% Output Weight = 1, Input Weight = 0
% Input Horizon = 5, Output Horizon = 20
% Minimum Constraint on Input = -0.4
% Maximum Constraint on Input = inf
% Delta Constraint on Input = 0.1
model = plant;
ywt = 1; uwt = 0;
M = 5; P = 20;
tend = 245;
r = 0;
ulim = [ ];
ylim = [ ]; tfilter = [ ]; dmodel = [ ];
dstep = 1;
[y9,u9] = cmpc(plant,model,ywt,uwt,M,P,tend,r,...
ulim,ylim,tfilter,dplant,dmodel,dstep);
ulim = [ -0.4, inf, 0.1]; % impose constraints
[y10,u10] = cmpc(plant,model,ywt,uwt,M,P,tend,r,...
ulim,ylim,tfilter,dplant,dmodel,dstep);
plotall([y9,y10],[u9,u10],delt2);

% Simulate and plot response for unmeasured
% step disturbance through dplant with uwt=0,
% with and without noise filtering.
tend=245;
r=[ ]; usat=[ ]; dmodel=[ ];
tfilter=[ ];
dstep=1;
[y5,u5]=mpcsim(plant,model,Kmpc1,tend,r,usat,tfilter,...
dplant,dmodel,dstep);
tfilter=20; % noise filtering time constant=20
[y6,u6]=mpcsim(plant,model,Kmpc1,tend,r,usat,tfilter,...
dplant,dmodel,dstep);
plotall([y5,y6],[u5,u6],delt2);

