addpath([pwd '\Delta_Variant']);
addpath([pwd '\Delta_Variant\Results']);

load('RAgTest_PlotOrder.mat');
NumTest=length(testName);

for ii=1:NumTest
    
[betaRTPCR,betaAg]=ParameterCOVIDTest(testName{ii},1);
[pA,~,~,ts,td] = BaselineParameters;
PPARAgA=zeros(15,1);
PPARAgS=zeros(15,1);
    for qq=1:15
        PPARAgA(qq)=(1./(td)).*integral(@(x)(TestSensitivity(x+(qq-1),ts,betaAg,betaRTPCR).*TestSensitivity(x+(qq-1),ts,[],betaRTPCR)),0,td);
        PPARAgS(qq)=(1./ts).*integral(@(x)(TestSensitivity(x+(qq-1),ts,betaAg,betaRTPCR).*(TestSensitivity(x+(qq-1),ts,[],betaRTPCR))),0,ts);
    end

    figure(ii)
    w=(ts-[0:14])./ts;
    w(w<0)=0;
    w=w';
    pAw=w.*pA+(1-w);
    plot([0:14],pAw.*PPARAgA+(1-pAw).*PPARAgS)
end