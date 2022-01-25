clear;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Average detection during incubation period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

addpath([pwd '\Delta_Variant']);
addpath([pwd '\Delta_Variant\Results']);

[~,IncubationI,~,ts,~] = BaselineParameters;
ts1=ts;
load('MLE-Estimate-RTPCR.mat','beta');


S = integral(@(t)TestSensitivity(t,ts,[],beta),0,ts)./ts;

load('RTPCR_Parameter_Uncertainty.mat','betaRTPCRv');
SU=zeros(length(betaRTPCRv(:,1)),1);
for ii=1:length(betaRTPCRv(:,1))
   SU(ii)=  integral(@(t)TestSensitivity(t,ts,[],betaRTPCRv(ii,:)),0,ts)./ts;
end

[~,LBS,UBS]=Credible_Interval_High_Density(S,SU,0.95,'continuous',[0 1]);
fprintf('=================================================================================== \n'); 
fprintf('Baseline \n');
fprintf('=================================================================================== \n');
fprintf('Percent of infection after the incubation period: %3.1f %% \n',100.*[1-IncubationI]);
fprintf('Average sensitivity during incubation period: %4.3f (%4.3f - %4.3f) \n',[S LBS UBS]);
rmpath([pwd '\Delta_Variant']);
rmpath([pwd '\Delta_Variant\Results']);

addpath([pwd '\Non_Delta']);
addpath([pwd '\Non_Delta\Results']);

[~,IncubationI,~,ts,~] = BaselineParameters;
ts1=ts;
load('MLE-Estimate-RTPCR_Non_Delta.mat','beta');


S = integral(@(t)TestSensitivity(t,ts,[],beta),0,ts)./ts;

load('RTPCR_Parameter_Uncertainty.mat','betaRTPCRv');
SU=zeros(length(betaRTPCRv(:,1)),1);
for ii=1:length(betaRTPCRv(:,1))
   SU(ii)=  integral(@(t)TestSensitivity(t,ts,[],betaRTPCRv(ii,:)),0,ts)./ts;
end

[~,LBS,UBS]=Credible_Interval_High_Density(S,SU,0.95,'continuous',[0 1]);
fprintf('=================================================================================== \n'); 
fprintf('Different incubation period and infectivity curve \n');
fprintf('=================================================================================== \n');
fprintf('Percent of infection after the incubation period: %3.1f %% \n',100.*[1-IncubationI]);
fprintf('Average sensitivity during incubation period: %4.3f (%4.3f - %4.3f) \n',[S LBS UBS]);
rmpath([pwd '\Non_Delta']);
rmpath([pwd '\Non_Delta\Results']);