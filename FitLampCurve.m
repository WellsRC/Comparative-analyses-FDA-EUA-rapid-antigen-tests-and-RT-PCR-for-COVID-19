 clear;
rng default
%Clinical Evaluation of Self-Collected Saliva by Quantitative Reverse
%Transcription-PCR (RT-qPCR), Direct RT-qPCR, Reverse Transcription–Loop-Mediated Isothermal Amplification, and a Rapid Antigen Test To Diagnose COVID-19
load('LAMP-Data.mat');
Time=LAMP(:,1);
Positive=LAMP(:,2);
TotalTest=LAMP(:,3);
LB=[-4 -4];
UB2=[log10(2) log10(2)];

tL=2.9;
ts=8.29;
options = optimoptions('ga','PlotFcn', @gaplotbestf,'MaxGenerations',5*10^3,'FunctionTolerance',10^(-16),'UseParallel',false);

[part,fvalt]=ga(@(x)LikelihoodLAMPCurve(x,Time+ts,Positive,TotalTest,tL),length(LB),[],[],[],[],LB,UB2,[],options);

opts= optimset('MaxIter',10^4,'MaxFunEvals',10^4,'TolFun',10^(-16),'TolX',10^(-16),'UseParallel',false,'Display','off');
[par,MLE]=fmincon(@(x)LikelihoodLAMPCurve(x,Time+ts,Positive,TotalTest,tL),part,[],[],[],[],LB,UB2,[],opts);
beta=10.^par;
save('MLE-Estimate-LAMP-Hill.mat');