function [betaRTPCR,betaAg]=ParameterCOVIDTest(testName,MLEv)
if(MLEv~=0)
    load('MLE-Estimate-RTPCR','beta');
    betaRTPCR=beta;
else
    load('MLE-Estimate-RTPCR','beta','MLE');
    load('Uncertainty-BetaEstimate-RTPCR.mat','L','betaRTU');
    betaRTU=[beta; betaRTU];
    L=[-MLE;L];
    w=cumsum(exp(L)./sum(exp(L)));
    r=rand(1);
    findx=find(w>=r,1);
    betaRTPCR=betaRTU(findx,:);
end

betaAg=[];
if(~isempty(testName))
    if(MLEv~=0)
        load([testName '_LR_Parameters.mat'],'beta');
        betaAg=beta;
    else
        load([testName '_LR_Uncertainty.mat'],'L','betaS','beta');
        w=cumsum(exp(L)./sum(exp(L)));
        r=rand(1);
        findx=find(w>=r,1);
        betaAg=beta.*(1+betaS(findx,:));
    end
end

end