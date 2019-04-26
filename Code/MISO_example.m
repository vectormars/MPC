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

close all
sim(MPCobj,Tf,r,v);


