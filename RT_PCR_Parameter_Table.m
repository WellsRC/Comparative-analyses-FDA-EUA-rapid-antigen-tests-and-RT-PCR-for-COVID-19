clear;
close all;
t=linspace(0,20,100001);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Baseline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
addpath([pwd '\Delta_Variant']);
addpath([pwd '\Delta_Variant\Results']);
load('MLE-Estimate-RTPCR.mat','beta');
S = max(PCRSens(t,beta));
PeakPCR(1,1)=S;
mmm=exp(beta(1)-beta(2)^2);
C(1,1)=beta(3)./lognpdf(mmm,beta(1),beta(2));
K(1,1)=beta(1);
z(1,1)=beta(2);

load('RTPCR_Parameter_Uncertainty.mat','betaRTPCRv');
S=zeros(size(betaRTPCRv(:,1)));
for ii=1:length(S)
   S(ii)=max(PCRSens(t,betaRTPCRv(ii,:)));
end
% NOTE: Rounded this instance for the MLE as numerically the upper bound
% was below the MLE.
[~,PeakPCR(1,2),PeakPCR(1,3)]=Credible_Interval_High_Density(round(PeakPCR(1,1),4),S,0.95,'continuous',[0 1]);
mmm=exp(betaRTPCRv(:,1)-betaRTPCRv(:,2).^2);
[~,C(1,2),C(1,3)]=Credible_Interval_High_Density(C(1,1),betaRTPCRv(:,3)./lognpdf(mmm,betaRTPCRv(:,1),betaRTPCRv(:,2)),0.95,'continuous',[0 max(betaRTPCRv(:,3)./lognpdf(mmm,betaRTPCRv(:,1),betaRTPCRv(:,2))).*1.01]);

[~,K(1,2),K(1,3)]=Credible_Interval_High_Density(K(1,1),betaRTPCRv(:,1),0.95,'continuous',[min(betaRTPCRv(:,1)).*0.995 max(betaRTPCRv(:,1)).*1.005]);

[~,z(1,2),z(1,3)]=Credible_Interval_High_Density(z(1,1),betaRTPCRv(:,2),0.95,'continuous',[min(betaRTPCRv(:,2)).*0.995 max(betaRTPCRv(:,2)).*1.005]);


rmpath([pwd '\Delta_Variant']);
rmpath([pwd '\Delta_Variant\Results']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Non-Delta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55

addpath([pwd '\Non_Delta']);
addpath([pwd '\Non_Delta\Results']);
load('MLE-Estimate-RTPCR_Non_Delta.mat','beta');

S = max(PCRSens(t,beta));
PeakPCR(3,1)=S;

mmm=exp(beta(1)-beta(2)^2);
C(3,1)=beta(3)./lognpdf(mmm,beta(1),beta(2));
K(3,1)=beta(1);
z(3,1)=beta(2);

load('RTPCR_Parameter_Uncertainty.mat','betaRTPCRv');


for ii=1:length(S)
   S(ii)=max(PCRSens(t,betaRTPCRv(ii,:)));
end
[~,PeakPCR(3,2),PeakPCR(3,3)]=Credible_Interval_High_Density(PeakPCR(3,1),S,0.95,'continuous',[0 1]);

mmm=exp(betaRTPCRv(:,1)-betaRTPCRv(:,2).^2);

[~,C(3,2),C(3,3)]=Credible_Interval_High_Density(C(3,1),betaRTPCRv(:,3)./lognpdf(mmm,betaRTPCRv(:,1),betaRTPCRv(:,2)),0.95,'continuous',[0 max(betaRTPCRv(:,3)./lognpdf(mmm,betaRTPCRv(:,1),betaRTPCRv(:,2))).*1.01]);

[~,K(3,2),K(3,3)]=Credible_Interval_High_Density(K(3,1),betaRTPCRv(:,1),0.95,'continuous',[min(betaRTPCRv(:,1)).*0.995 max(betaRTPCRv(:,1)).*1.005]);

[~,z(3,2),z(3,3)]=Credible_Interval_High_Density(z(3,1),betaRTPCRv(:,2),0.95,'continuous',[min(betaRTPCRv(:,2)).*0.995 max(betaRTPCRv(:,2)).*1.005]);

rmpath([pwd '\Non_Delta']);
rmpath([pwd '\Non_Delta\Results']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Alternative curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
addpath([pwd '\Alternative_Curve_Delta_Variant']);
addpath([pwd '\Alternative_Curve_Delta_Variant\Results']);
load('MLE-Estimate-RTPCR_Alternative.mat','beta');

S = max(PCRSens(t,beta));
PeakPCR(2,1)=S;

mmm=1./beta(1);
C(2,1)=beta(3)./tpdf(log(mmm.*beta(1)),beta(2));
K(2,1)=beta(1);
z(2,1)=beta(2);

load('RTPCR_Parameter_Uncertainty.mat','betaRTPCRv');



for ii=1:length(S)
   S(ii)=max(PCRSens(t,betaRTPCRv(ii,:)));
end
[~,PeakPCR(2,2),PeakPCR(2,3)]=Credible_Interval_High_Density(PeakPCR(2,1),S,0.95,'continuous',[0 1]);

mmm=1./betaRTPCRv(:,1);


[~,C(2,2),C(2,3)]=Credible_Interval_High_Density(C(2,1),betaRTPCRv(:,3)./tpdf(log(mmm.*betaRTPCRv(:,1)),betaRTPCRv(:,2)),0.95,'continuous',[0 max(betaRTPCRv(:,3)./tpdf(log(mmm.*betaRTPCRv(:,1)),betaRTPCRv(:,2))).*1.01]);

% THERE IS NO VARIATION IN THIS PARAMETER HERE
K(2,2)=K(2,1);
K(2,3)=K(2,1);

[~,z(2,2),z(2,3)]=Credible_Interval_High_Density(z(2,1),betaRTPCRv(:,2),0.95,'continuous',[min(betaRTPCRv(:,2)).*0.995 max(betaRTPCRv(:,2)).*1.005]);

rmpath([pwd '\Alternative_Curve_Delta_Variant']);
rmpath([pwd '\Alternative_Curve_Delta_Variant\Results']);


Inc={'4.4'; '4.4'; '5.72'};
FuncForm={'log-Normal'; 'log-Student t'; 'log-Normal'};

K95=cell(3,1);
z95=cell(3,1);
C95=cell(3,1);
Peak95=cell(3,1);
for ii=1:3
    K95{ii}=[num2str(K(ii,1),'%4.3g') ' (' num2str(K(ii,2),'%4.3g')  char(8211) num2str(K(ii,3),'%4.3g') ')'];
    z95{ii}=[num2str(z(ii,1),'%4.3g') ' (' num2str(z(ii,2),'%4.3g')  char(8211) num2str(z(ii,3),'%4.3g') ')'];
    C95{ii}=[num2str(C(ii,1),'%4.3g') ' (' num2str(C(ii,2),'%4.3g')  char(8211) num2str(C(ii,3),'%4.3g') ')'];
    Peak95{ii}=[num2str(PeakPCR(ii,1),'%4.3g') ' (' num2str(PeakPCR(ii,2),'%4.3g')  char(8211) num2str(PeakPCR(ii,3),'%4.3g') ')'];
end

T=table(Inc,FuncForm,K95,z95,C95,Peak95);

writetable(T,'Table_RT-PCR_paramters.csv');