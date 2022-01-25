% Determiens the sensitivity of adding an entry to to exit test
clear;
clc;
close all;

addpath([pwd '\Delta_Variant']);
addpath([pwd '\Delta_Variant\Results']);

[pA,~,R0,ts,td] = BaselineParameters;

load('RAgTest_Name.mat','testName');
NumTests=length(testName);
S_Symp_EX=zeros(NumTests,14);
S_Asymp_EX=zeros(NumTests,14);

S_Symp_X=zeros(NumTests,14);
S_Asymp_X=zeros(NumTests,14);
for TestN=1:NumTests
    
    [betaRTPCR,betaAg]=ParameterCOVIDTest(testName{TestN},1);
    testtype=cell(2,1);
    testtype{1}=betaAg;
    testtype{2}=betaAg;
    for qq=1:14
        S_Symp_EX(TestN,qq) = 1-(1./ts).*integral(@(t)(1-TestSensitivity(t,ts,testtype{1},betaRTPCR)).*(1-TestSensitivity(t+qq,ts,testtype{2},betaRTPCR)),0,ts); 
        S_Asymp_EX(TestN,qq) = 1-(1./td).*integral(@(t)(1-TestSensitivity(t,ts,testtype{1},betaRTPCR)).*(1-TestSensitivity(t+qq,ts,testtype{2},betaRTPCR)),0,td);
        
        
        S_Symp_X(TestN,qq) = 1-(1./ts).*integral(@(t)(1-TestSensitivity(t+qq,ts,testtype{2},betaRTPCR)),0,ts); 
        S_Asymp_X(TestN,qq) = 1-(1./td).*integral(@(t)(1-TestSensitivity(t+qq,ts,testtype{2},betaRTPCR)),0,td);
    end
end


rmpath([pwd '\Delta_Variant']);
rmpath([pwd '\Delta_Variant\Results']);
