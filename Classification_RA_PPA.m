
addpath('Delta_Variant/Results');
addpath('Delta_Variant');

load('RAgTest_Name.mat','testName');


StablePPA={};
GradualPPA={};
RapidPPA={};
RC=1;
SC=1;
GC=1;
t=[0 20 40];
for TestN=1:length(testName)
    [betaRTPCR,betaAg]=ParameterCOVIDTest(testName{TestN},1);
    temp=LR(t,betaAg);
    if(temp(2)<0.01)
        RapidPPA{RC}=testName{TestN};
        RC=RC+1;
    elseif(1-temp(3)./temp(1)<0.01)
        StablePPA{SC}=testName{TestN};
        SC=SC+1;
    else
        GradualPPA{GC}=testName{TestN};
        GC=GC+1;
    end
end
rmpath('Delta_Variant/Results');
rmpath('Delta_Variant');