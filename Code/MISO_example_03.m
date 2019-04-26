clc;
clear;
sys=ss(tf({1,1,1},{[1 .5 1],[1 1],[.7 .5 1]}),'min');
Ts=.2;             % sampling time
model=c2d(sys,Ts); % discrete model

model=setmpcsignals(model,'MV',1,'MD',2,'UD',3);


%  Define the structure of models used by the MPC controller. %
clear Model
Model.Plant=model;  % Predictive model
Model.Disturbance=tf(sqrt(1000),[1 0]); 
% Disturbance model: Integrator driven by white noise with variance = 1000

%  Define prediction and control horizons. %
p=[];       % prediction horizon (take default one)
m=3;        % control horizon

MPCobj=mpc(Model,Ts,p,m);    % define objective function
MPCobj.MV=struct('Min',0,'Max',1,'RateMin',-10,'RateMax',10);
% define constraints

Tstop=30;                               % simulation time
Tf=round(Tstop/Ts);                     % number of simulation steps=150
r=ones(Tf,1);                           % reference trajectory
v=[zeros(Tf/3,1);ones(2*Tf/3,1)];       % measured disturbance trajectory
% First 50 steps=0, rest 100 steps=1  ==>in MD1

d=[zeros(2*Tf/3,1);-0.5*ones(Tf/3,1)]; % unmeasured disturbance trajectory
% first 100 steps=0,rest 50 steps=-0.5 
SimOptions=mpcsimopt(MPCobj);
SimOptions.Unmeas=d;                          % unmeasured input disturbance
SimOptions.OutputNoise=.001*(rand(Tf,1)-.5);  % output measurement noise
SimOptions.InputNoise=.05*(rand(Tf,1)-.5);    % noise on manipulated variables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

simModel=ss(tf({1,1,1},{[1 .8 1],[1 2],[.6 .6 1]}),'min');
simModel=setmpcsignals(simModel,'MV',1,'MD',2,'UD',3);
simModel=struct('Plant',simModel);
simModel.Nominal.Y=0.1; % The nominal value of the output of the true plant is 0.1
simModel.Nominal.X=-.1*[1 1 1 1 1];

SimOptions.Model=simModel;
SimOptions.plantinit=[0.1 0 -0.1 0 .05]; % Initial state of the true plant
SimOptions.OutputNoise=[];  % remove output measurement noise
SimOptions.InputNoise=[];   % remove noise on manipulated variables

close all
sim(MPCobj,Tf,r,v,SimOptions);