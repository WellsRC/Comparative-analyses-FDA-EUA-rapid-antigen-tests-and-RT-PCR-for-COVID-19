clear;
clc;
addpath([pwd '\Non_Delta']);
addpath([pwd '\Non_Delta\Results']);

P=[3 1 0];
N=[12 5 3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clopper Pearson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Entry test PPA: %3.1f%% (%3.1f%%-%3.1f%%) \n',100.*[P(1)./N(1) betainv(0.025,P(1),N(1)-P(1)+1) betainv(0.975,P(1)+1,N(1)-P(1))]);
fprintf('Day 3 test PPA: %3.1f%% (%3.1f%%-%3.1f%%) \n',100.*[P(2)./N(2) betainv(0.025,P(2),N(2)-P(2)+1) betainv(0.975,P(2)+1,N(2)-P(2))]);

% Day 4 is DIFFERENT B/C IT WAS ZERO!!!!!
fprintf('Day 4 test PPA: %3.1f%% (%3.1f%%-%3.1f%%) \n',100.*[P(3)./N(3) 0 1-(0.025)^(1./N(3))]);

[betaRTPCR,betaAg]=ParameterCOVIDTest('BD Veritor',1);
[pA,~,~,ts,td] = BaselineParameters;

test=[0 3 4];

PPARAgA=zeros(1,length(test));

PPARAgS=zeros(1,length(test));


    ii=1;
    PPARAgA(ii)=(1./(td)).*integral(@(x)(TestSensitivity(x+test(ii),ts,betaAg,betaRTPCR).*TestSensitivity(x+test(ii),ts,[],betaRTPCR)),0,td);

    PPARAgS(ii)=(1./ts).*integral(@(x)(TestSensitivity(x+test(ii),ts,betaAg,betaRTPCR).*(TestSensitivity(x+test(ii),ts,[],betaRTPCR))),0,ts);

    C=integral(@(x)(1-TestSensitivity(x+test(ii),ts,[],betaRTPCR)),0,td);
    
    CS=integral(@(x)(1-TestSensitivity(x+test(ii),ts,[],betaRTPCR)),0,ts-test(2));
    
    ii=2;

    % RT-PCR is when the individual is isolated
    PPARAgA(ii)=(1./C).*integral(@(x)((1-TestSensitivity(x+test(1),ts,[],betaRTPCR)).*TestSensitivity(x+test(ii),ts,betaAg,betaRTPCR).*TestSensitivity(x+test(ii),ts,[],betaRTPCR)),0,td);
    
    
    PPARAgS(ii)=(1./CS).*integral(@(x)((1-TestSensitivity(x+test(1),ts,[],betaRTPCR)).*TestSensitivity(x+test(ii),ts,betaAg,betaRTPCR).*TestSensitivity(x+test(ii),ts,[],betaRTPCR)),0,ts-test(ii));


    % RT-PCR is when the individual is isolated
    C=integral(@(x)(1-TestSensitivity(x+test(1),ts,[],betaRTPCR)).*(1-TestSensitivity(x+test(2),ts,[],betaRTPCR)),0,td);
    CS=integral(@(x)(1-TestSensitivity(x+test(1),ts,[],betaRTPCR)).*(1-TestSensitivity(x+test(2),ts,[],betaRTPCR)),0,ts-test(3));
    ii=3;


    PPARAgA(ii)=(1./C).*integral(@(x)((1-TestSensitivity(x+test(1),ts,[],betaRTPCR)).*(1-TestSensitivity(x+test(2),ts,[],betaRTPCR)).*TestSensitivity(x+test(ii),ts,betaAg,betaRTPCR).*TestSensitivity(x+test(ii),ts,[],betaRTPCR)),0,td);    
    
    PPARAgS(ii)=(1./CS).*integral(@(x)((1-TestSensitivity(x+test(1),ts,[],betaRTPCR)).*(1-TestSensitivity(x+test(2),ts,[],betaRTPCR)).*TestSensitivity(x+test(ii),ts,betaAg,betaRTPCR).*TestSensitivity(x+test(ii),ts,[],betaRTPCR)),0,ts-test(ii));
    
    
PPARAgA_UN=zeros(1000,length(test));

PPARAgS_UN=zeros(1000,length(test));

for jj=1:1000
    [betaRTPCR,betaAg]=ParameterCOVIDTest('BD Veritor',0);
    ii=1;
    PPARAgA_UN(jj,ii)=(1./(td)).*integral(@(x)(TestSensitivity(x+test(ii),ts,betaAg,betaRTPCR).*TestSensitivity(x+test(ii),ts,[],betaRTPCR)),0,td);

    PPARAgS_UN(jj,ii)=(1./ts).*integral(@(x)(TestSensitivity(x+test(ii),ts,betaAg,betaRTPCR).*(TestSensitivity(x+test(ii),ts,[],betaRTPCR))),0,ts);

    C=integral(@(x)(1-TestSensitivity(x+test(ii),ts,[],betaRTPCR)),0,td);
    
    CS=integral(@(x)(1-TestSensitivity(x+test(ii),ts,[],betaRTPCR)),0,ts-test(2));
    
    ii=2;

    % RT-PCR is when the individual is isolated
    PPARAgA_UN(jj,ii)=(1./C).*integral(@(x)((1-TestSensitivity(x+test(1),ts,[],betaRTPCR)).*TestSensitivity(x+test(ii),ts,betaAg,betaRTPCR).*TestSensitivity(x+test(ii),ts,[],betaRTPCR)),0,td);
    
    
    PPARAgS_UN(jj,ii)=(1./CS).*integral(@(x)((1-TestSensitivity(x+test(1),ts,[],betaRTPCR)).*TestSensitivity(x+test(ii),ts,betaAg,betaRTPCR).*TestSensitivity(x+test(ii),ts,[],betaRTPCR)),0,ts-test(ii));


    % RT-PCR is when the individual is isolated
    C=integral(@(x)(1-TestSensitivity(x+test(1),ts,[],betaRTPCR)).*(1-TestSensitivity(x+test(2),ts,[],betaRTPCR)),0,td);
    CS=integral(@(x)(1-TestSensitivity(x+test(1),ts,[],betaRTPCR)).*(1-TestSensitivity(x+test(2),ts,[],betaRTPCR)),0,ts-test(3));
    ii=3;


    PPARAgA_UN(jj,ii)=(1./C).*integral(@(x)((1-TestSensitivity(x+test(1),ts,[],betaRTPCR)).*(1-TestSensitivity(x+test(2),ts,[],betaRTPCR)).*TestSensitivity(x+test(ii),ts,betaAg,betaRTPCR).*TestSensitivity(x+test(ii),ts,[],betaRTPCR)),0,td);    
    
    PPARAgS_UN(jj,ii)=(1./CS).*integral(@(x)((1-TestSensitivity(x+test(1),ts,[],betaRTPCR)).*(1-TestSensitivity(x+test(2),ts,[],betaRTPCR)).*TestSensitivity(x+test(ii),ts,betaAg,betaRTPCR).*TestSensitivity(x+test(ii),ts,[],betaRTPCR)),0,ts-test(ii));
    
end

w=[ts ts-3 ts-4]./ts;
pAw=w.*pA+(1-w);
PPA=pAw.*PPARAgA+(1-pAw).*PPARAgS;
PPAU=repmat(pAw,1000,1).*PPARAgA_UN+(1-repmat(pAw,1000,1)).*PPARAgS_UN;

for ii=1:length(test)
    fprintf('Test day %d the PPA is %3.1f (%3.1f - %3.1f) \n', [test(ii) 100.*PPA(ii) 100.*prctile(PPAU(:,ii),[2.5 97.5])]);
end
rmpath([pwd '\Non_Delta']);
rmpath([pwd '\Non_Delta\Results']);