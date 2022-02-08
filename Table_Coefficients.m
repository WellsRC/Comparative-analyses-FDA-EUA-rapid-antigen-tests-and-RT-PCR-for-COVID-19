%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Produces the output for table 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
clear;
clc;
addpath([pwd '\Delta_Variant']);
addpath([pwd '\Delta_Variant\Results']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Produces the output for table 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

load('RAgTest_Name_Order.mat');

NSS=1000;
Ntest=length(testName);

beta0_Text=cell(Ntest,1);
beta1_Text=cell(Ntest,1);
b0=zeros(1,NSS);
b1=zeros(1,NSS);
for ii=1:Ntest
    [~,betaAg]=ParameterCOVIDTest(testName{ii},1);
    MLE_beta0=betaAg(1);
    MLE_beta1=betaAg(2);
    parfor jj=1:NSS        
        [~,betaAg]=ParameterCOVIDTest(testName{ii},0);
        b0(jj)=betaAg(1);
        b1(jj)=betaAg(2);        
    end    
    [~,b0l,b0u]=Credible_Interval_High_Density(MLE_beta0,b0,0.95,'continuous',[-10^(-12) max(b0)+10]);
    [~,b1l,b1u]=Credible_Interval_High_Density(MLE_beta1,b1,0.95,'continuous',[min(b1)-10 10^(-12)]);
    beta0_Text{ii}=[num2str(MLE_beta0,'%4.3g') ' (' num2str(b0l,'%4.3g')  char(8211) num2str(b0u,'%4.3g') ')'];
    beta1_Text{ii}=[num2str(MLE_beta1,'%4.3g') ' (' num2str(b1l,'%4.3g')  char(8211) num2str(b1u,'%4.3g') ')'];
end

T2=table(testName,beta0_Text,beta1_Text);


rmpath([pwd '\Delta_Variant']);
rmpath([pwd '\Delta_Variant\Results']);

writetable(T2,'Table_PPA_Coefficients.csv');
