function [QDurE,QDurEE,FreqT,FPProb] = Quarantine_Testing_Freq(testName,pA,DurT,Scenario)

if(strcmp(Scenario,'Baseline'))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Determine the equivilant quarantine
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    Risk=1;
    load('TestingonExit_RTPCR_24hrDelay_DeltaVOC.mat','q','IDSLS','IDSLA');
    RRTPCR=Probability_Onward((1-pA).*IDSLS(q==7)+pA.*IDSLA(q==7),Risk);

    load(['TestingonExit_' testName '_NoDelay_DeltaVOC.mat'],'q','IDSLS','IDSLA')
    RAg=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    QDurE(1)=min(q(RAg<RRTPCR));

    load(['TestingonEntryExit_' testName '_NoDelay_DeltaVOC.mat'],'q','IDSLS','IDSLA')
    RAg=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    QDurEE(1)=min(q(RAg<RRTPCR));
    
    
    
    load('TestingonExit_RTPCR_24hrDelay_DeltaVOC_Uncertainty.mat','q','IDSLSv','IDSLAv')
    RRTPCR=Probability_Onward((1-pA).*IDSLSv(:,q==7)+pA.*IDSLAv(:,q==7),Risk);
    
    load(['TestingonExit_' testName '_NoDelay_DeltaVOC_Uncertainty.mat'],'q','IDSLSv','IDSLAv')
    RAg=Probability_Onward((1-pA).*IDSLSv+pA.*IDSLAv,Risk);
    
    load(['TestingEntryExit_' testName '_NoDelay_DeltaVOC_Uncertainty.mat'],'q','IDSLSv','IDSLAv')
    R2Ag=Probability_Onward((1-pA).*IDSLSv+pA.*IDSLAv,Risk);
    
    TestE=zeros(length(RRTPCR),1);
    TestEE=zeros(length(RRTPCR),1);
    for jj=1:length(RRTPCR)
        TestE(jj)=min(q(RAg(jj,:)<RRTPCR(jj)));
        TestEE(jj)=min(q(R2Ag(jj,:)<RRTPCR(jj)));
    end
    
    QDurE(2)=prctile(TestE,2.5);
    QDurE(3)=prctile(TestE,97.5);
    
    
    QDurEE(2)=prctile(TestEE,2.5);
    QDurEE(3)=prctile(TestEE,97.5);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Determine the minimal testing frequency
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    ft=[1:14];

    load(['Testing_Frequency_' testName '_DeltaVOC.mat'],'RTotA','RTotS')
    R=(1-pA).*RTotS+pA.*RTotA;
    FreqT=max(ft(R<1)); % take the maximum time between tests

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Risk of false positve
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    FPProb=CalcFalsePositive(testName,DurT,FreqT); % Calacute the probability of a 
elseif(strcmp(Scenario,'Non_Delta'))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Determine the equivilant quarantine
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    Risk=1;
    load('TestingonExit_RTPCR_24hrDelay_General.mat','q','IDSLS','IDSLA');
    RRTPCR=Probability_Onward((1-pA).*IDSLS(q==7)+pA.*IDSLA(q==7),Risk);

    load(['TestingonExit_' testName '_NoDelay_General.mat'],'q','IDSLS','IDSLA')
    RAg=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    QDurE(1)=min(q(RAg<RRTPCR));

    load(['TestingonEntryExit_' testName '_NoDelay_General.mat'],'q','IDSLS','IDSLA')
    RAg=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    QDurEE(1)=min(q(RAg<RRTPCR));
    
    
    
    load('TestingonExit_RTPCR_24hrDelay_General_Uncertainty.mat','q','IDSLSv','IDSLAv')
    RRTPCR=Probability_Onward((1-pA).*IDSLSv(:,q==7)+pA.*IDSLAv(:,q==7),Risk);
    
    load(['TestingonExit_' testName '_NoDelay_General_Uncertainty.mat'],'q','IDSLSv','IDSLAv')
    RAg=Probability_Onward((1-pA).*IDSLSv+pA.*IDSLAv,Risk);
    
    load(['TestingEntryExit_' testName '_NoDelay_General_Uncertainty.mat'],'q','IDSLSv','IDSLAv')
    R2Ag=Probability_Onward((1-pA).*IDSLSv+pA.*IDSLAv,Risk);
    
    TestE=zeros(length(RRTPCR),1);
    TestEE=zeros(length(RRTPCR),1);
    for jj=1:length(RRTPCR)
        TestE(jj)=min(q(RAg(jj,:)<RRTPCR(jj)));
        TestEE(jj)=min(q(R2Ag(jj,:)<RRTPCR(jj)));
    end
    
    QDurE(2)=prctile(TestE,2.5);
    QDurE(3)=prctile(TestE,97.5);
    
    
    QDurEE(2)=prctile(TestEE,2.5);
    QDurEE(3)=prctile(TestEE,97.5);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Determine the minimal testing frequency
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    ft=[1:14];

    load(['Testing_Frequency_' testName '_General.mat'],'RTotA','RTotS')
    R=(1-pA).*RTotS+pA.*RTotA;
    FreqT=max(ft(R<1)); % take the maximum time between tests

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Risk of false positve
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    FPProb=CalcFalsePositive(testName,DurT,FreqT); % Calacute the probability of a 
elseif(strcmp(Scenario,'Alternative_RTPCR'))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Determine the equivilant quarantine
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    Risk=1;
    load('TestingonExit_RTPCR_24hrDelay_DeltaVOC_Alternative_PCR.mat','q','IDSLS','IDSLA');
    RRTPCR=Probability_Onward((1-pA).*IDSLS(q==7)+pA.*IDSLA(q==7),Risk);

    load(['TestingonExit_' testName '_NoDelay_DeltaVOC_Alternative_PCR.mat'],'q','IDSLS','IDSLA')
    RAg=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    QDurE(1)=min(q(RAg<RRTPCR));

    load(['TestingonEntryExit_' testName '_NoDelay_DeltaVOC_Alternative_PCR.mat'],'q','IDSLS','IDSLA')
    RAg=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    QDurEE(1)=min(q(RAg<RRTPCR));
    
    
    
    load('TestingonExit_RTPCR_24hrDelay_DeltaVOC_Alternative_PCR_Uncertainty.mat','q','IDSLSv','IDSLAv')
    RRTPCR=Probability_Onward((1-pA).*IDSLSv(:,q==7)+pA.*IDSLAv(:,q==7),Risk);
    
    load(['TestingonExit_' testName '_NoDelay_DeltaVOC_Alternative_PCR_Uncertainty.mat'],'q','IDSLSv','IDSLAv')
    RAg=Probability_Onward((1-pA).*IDSLSv+pA.*IDSLAv,Risk);
    
    load(['TestingEntryExit_' testName '_NoDelay_DeltaVOC_Alternative_PCR_Uncertainty.mat'],'q','IDSLSv','IDSLAv')
    R2Ag=Probability_Onward((1-pA).*IDSLSv+pA.*IDSLAv,Risk);
    
    TestE=zeros(length(RRTPCR),1);
    TestEE=zeros(length(RRTPCR),1);
    for jj=1:length(RRTPCR)
        TestE(jj)=min(q(RAg(jj,:)<RRTPCR(jj)));
        TestEE(jj)=min(q(R2Ag(jj,:)<RRTPCR(jj)));
    end
    
    QDurE(2)=prctile(TestE,2.5);
    QDurE(3)=prctile(TestE,97.5);
    
    
    QDurEE(2)=prctile(TestEE,2.5);
    QDurEE(3)=prctile(TestEE,97.5);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Determine the minimal testing frequency
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    ft=[1:14];

    load(['Testing_Frequency_' testName '_DeltaVOC_Alternative_PCR.mat'],'RTotA','RTotS')
    R=(1-pA).*RTotS+pA.*RTotA;
    FreqT=max(ft(R<1)); % take the maximum time between tests

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Risk of false positve
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    FPProb=CalcFalsePositive(testName,DurT,FreqT); % Calacute the probability of a 
end
end


