
load('RAgTest_PlotOrder.mat');


addpath([pwd '\Delta_Variant']);
addpath([pwd '\Delta_Variant\Results']);
Ntest=length(testName);


[pA,~,~,~,~] = BaselineParameters; % Proportion of infections being asymptomatic
    TestE=zeros(1000,Ntest);
    TestEE=zeros(1000,Ntest);
    
for ii=1:Ntest
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Determine the equivilant quarantine
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    Risk=1;   
    
    load('TestingonExit_RTPCR_24hrDelay_DeltaVOC_Uncertainty.mat','q','IDSLSv','IDSLAv')
    RRTPCR=Probability_Onward((1-pA).*IDSLSv(:,q==7)+pA.*IDSLAv(:,q==7),Risk);
    
    load(['TestingonExit_' testName{ii} '_NoDelay_DeltaVOC_Uncertainty.mat'],'q','IDSLSv','IDSLAv')
    RAg=Probability_Onward((1-pA).*IDSLSv+pA.*IDSLAv,Risk);
    
    load(['TestingEntryExit_' testName{ii} '_NoDelay_DeltaVOC_Uncertainty.mat'],'q','IDSLSv','IDSLAv')
    R2Ag=Probability_Onward((1-pA).*IDSLSv+pA.*IDSLAv,Risk);
    
    for jj=1:length(RRTPCR)
        TestE(jj,ii)=min(q(RAg(jj,:)<RRTPCR(jj)));
        TestEE(jj,ii)=min(q(R2Ag(jj,:)<RRTPCR(jj)));
    end
    
end


rmpath([pwd '\Delta_Variant']);
rmpath([pwd '\Delta_Variant\Results']);
