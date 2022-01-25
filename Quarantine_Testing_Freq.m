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

    QDurE(1)=min(q(RAg<=RRTPCR));

    load(['TestingonEntryExit_' testName '_NoDelay_DeltaVOC.mat'],'q','IDSLS','IDSLA')
    RAg=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    QDurEE(1)=min(q(RAg<=RRTPCR));
    
    
    
    load('TestingonExit_RTPCR_24hrDelay_DeltaVOC_Uncertainty.mat','q','IDSLSv','IDSLAv')
    RRTPCR=Probability_Onward((1-pA).*IDSLSv(:,q==7)+pA.*IDSLAv(:,q==7),Risk);
    
    load(['TestingonExit_' testName '_NoDelay_DeltaVOC_Uncertainty.mat'],'q','IDSLSv','IDSLAv')
    RAg=Probability_Onward((1-pA).*IDSLSv+pA.*IDSLAv,Risk);
    
    load(['TestingEntryExit_' testName '_NoDelay_DeltaVOC_Uncertainty.mat'],'q','IDSLSv','IDSLAv')
    R2Ag=Probability_Onward((1-pA).*IDSLSv+pA.*IDSLAv,Risk);
    
    TestE=zeros(length(RRTPCR),1);
    TestEE=zeros(length(RRTPCR),1);
    for jj=1:length(RRTPCR)
        TestE(jj)=min(q(RAg(jj,:)<=RRTPCR(jj)));
        TestEE(jj)=min(q(R2Ag(jj,:)<=RRTPCR(jj)));
    end
    [~,QDurE(2),QDurE(3)]=Credible_Interval_High_Density(QDurE(1),TestE,0.95,'discrete',[1 14]);    
    [~,QDurEE(2),QDurEE(3)]=Credible_Interval_High_Density(QDurEE(1),TestEE,0.95,'discrete',[1 14]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Determine the minimal testing frequency
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    ft=[1:14];

    load(['Testing_Frequency_' testName '_DeltaVOC.mat'],'RTotA','RTotS')
    R=(1-pA).*RTotS+pA.*RTotA;
    FreqT(1)=max(ft(R<1)); % take the maximum time between tests
    
    load(['Testing_Frequency_' testName '_DeltaVOC_Uncertainty.mat'],'RTotAv','RTotSv')
    R=(1-pA).*RTotSv+pA.*RTotAv;
    TestF=zeros(length(R(1,:)),1);
    FPU=zeros(length(R(1,:)),250);
    for ii=1:length(TestF)
       TestF(ii)=max(ft(R(:,ii)<1)); 
       parfor jj=1:250
            FPU(ii,jj)=CalcFalsePositive(testName,DurT,TestF(ii),0);
       end
    end
    
    [~,FreqT(2),FreqT(3)]=Credible_Interval_High_Density(FreqT(1),TestF,0.95,'discrete',[1 14]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Risk of false positve
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    FPProb(1)=CalcFalsePositive(testName,DurT,FreqT(1),1); % Calacute the probability of a 
    [~,FPProb(2),FPProb(3)]=Credible_Interval_High_Density(FPProb(1),FPU(:),0.95,'continuous',[0 1]);
elseif(strcmp(Scenario,'Non_Delta'))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Determine the equivilant quarantine
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    Risk=1;
    load('TestingonExit_RTPCR_24hrDelay_General.mat','q','IDSLS','IDSLA');
    RRTPCR=Probability_Onward((1-pA).*IDSLS(q==7)+pA.*IDSLA(q==7),Risk);

    load(['TestingonExit_' testName '_NoDelay_General.mat'],'q','IDSLS','IDSLA')
    RAg=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    QDurE(1)=min(q(RAg<=RRTPCR));

    load(['TestingonEntryExit_' testName '_NoDelay_General.mat'],'q','IDSLS','IDSLA')
    RAg=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    QDurEE(1)=min(q(RAg<=RRTPCR));
    
    
    
    load('TestingonExit_RTPCR_24hrDelay_General_Uncertainty.mat','q','IDSLSv','IDSLAv')
    RRTPCR=Probability_Onward((1-pA).*IDSLSv(:,q==7)+pA.*IDSLAv(:,q==7),Risk);
    
    load(['TestingonExit_' testName '_NoDelay_General_Uncertainty.mat'],'q','IDSLSv','IDSLAv')
    RAg=Probability_Onward((1-pA).*IDSLSv+pA.*IDSLAv,Risk);
    
    load(['TestingEntryExit_' testName '_NoDelay_General_Uncertainty.mat'],'q','IDSLSv','IDSLAv')
    R2Ag=Probability_Onward((1-pA).*IDSLSv+pA.*IDSLAv,Risk);
    
    TestE=zeros(length(RRTPCR),1);
    TestEE=zeros(length(RRTPCR),1);
    for jj=1:length(RRTPCR)
        TestE(jj)=min(q(RAg(jj,:)<=RRTPCR(jj)));
        TestEE(jj)=min(q(R2Ag(jj,:)<=RRTPCR(jj)));
    end
    
    [~,QDurE(2),QDurE(3)]=Credible_Interval_High_Density(QDurE(1),TestE,0.95,'discrete',[1 14]);    
    [~,QDurEE(2),QDurEE(3)]=Credible_Interval_High_Density(QDurEE(1),TestEE,0.95,'discrete',[1 14]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Determine the minimal testing frequency
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    ft=[1:14];

    load(['Testing_Frequency_' testName '_General.mat'],'RTotA','RTotS')
    R=(1-pA).*RTotS+pA.*RTotA;
    FreqT(1)=max(ft(R<1)); % take the maximum time between tests
    
    load(['Testing_Frequency_' testName '_General_Uncertainty.mat'],'RTotAv','RTotSv')
    R=(1-pA).*RTotSv+pA.*RTotAv;
    TestF=zeros(length(R(1,:)),1);
    
    FPU=zeros(length(R(1,:)),250);
    for ii=1:length(R)
       TestF(ii)=max(ft(R(:,ii)<1)); 
       parfor jj=1:250
            FPU(ii,jj)=CalcFalsePositive(testName,DurT,TestF(ii),0);
       end
    end
    
    [~,FreqT(2),FreqT(3)]=Credible_Interval_High_Density(FreqT(1),TestF,0.95,'discrete',[1 14]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Risk of false positve
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    FPProb(1)=CalcFalsePositive(testName,DurT,FreqT(1),1); % Calacute the probability of a 
    [~,FPProb(2),FPProb(3)]=Credible_Interval_High_Density(FPProb(1),FPU(:),0.95,'continuous',[0 1]);
elseif(strcmp(Scenario,'Alternative_RTPCR'))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Determine the equivilant quarantine
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    Risk=1;
    load('TestingonExit_RTPCR_24hrDelay_DeltaVOC_Alternative_PCR.mat','q','IDSLS','IDSLA');
    RRTPCR=Probability_Onward((1-pA).*IDSLS(q==7)+pA.*IDSLA(q==7),Risk);

    load(['TestingonExit_' testName '_NoDelay_DeltaVOC_Alternative_PCR.mat'],'q','IDSLS','IDSLA')
    RAg=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    QDurE(1)=min(q(RAg<=RRTPCR));

    load(['TestingonEntryExit_' testName '_NoDelay_DeltaVOC_Alternative_PCR.mat'],'q','IDSLS','IDSLA')
    RAg=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    QDurEE(1)=min(q(RAg<=RRTPCR));
    
    
    
    load('TestingonExit_RTPCR_24hrDelay_DeltaVOC_Alternative_PCR_Uncertainty.mat','q','IDSLSv','IDSLAv')
    RRTPCR=Probability_Onward((1-pA).*IDSLSv(:,q==7)+pA.*IDSLAv(:,q==7),Risk);
    
    load(['TestingonExit_' testName '_NoDelay_DeltaVOC_Alternative_PCR_Uncertainty.mat'],'q','IDSLSv','IDSLAv')
    RAg=Probability_Onward((1-pA).*IDSLSv+pA.*IDSLAv,Risk);
    
    load(['TestingEntryExit_' testName '_NoDelay_DeltaVOC_Alternative_PCR_Uncertainty.mat'],'q','IDSLSv','IDSLAv')
    R2Ag=Probability_Onward((1-pA).*IDSLSv+pA.*IDSLAv,Risk);
    
    TestE=zeros(length(RRTPCR),1);
    TestEE=zeros(length(RRTPCR),1);
    for jj=1:length(RRTPCR)
        TestE(jj)=min(q(RAg(jj,:)<=RRTPCR(jj)));
        TestEE(jj)=min(q(R2Ag(jj,:)<=RRTPCR(jj)));
    end
    
    [~,QDurE(2),QDurE(3)]=Credible_Interval_High_Density(QDurE(1),TestE,0.95,'discrete',[1 14]);    
    [~,QDurEE(2),QDurEE(3)]=Credible_Interval_High_Density(QDurEE(1),TestEE,0.95,'discrete',[1 14]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Determine the minimal testing frequency
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    ft=[1:14];

    load(['Testing_Frequency_' testName '_DeltaVOC_Alternative_PCR.mat'],'RTotA','RTotS')
    R=(1-pA).*RTotS+pA.*RTotA;
    FreqT(1)=max(ft(R<1)); % take the maximum time between tests
    
    load(['Testing_Frequency_' testName '_DeltaVOC_Alternative_PCR_Uncertainty.mat'],'RTotAv','RTotSv')
    R=(1-pA).*RTotSv+pA.*RTotAv;
    TestF=zeros(length(R(1,:)),1);
    
    FPU=zeros(length(R(1,:)),250);
    for ii=1:length(R)
       TestF(ii)=max(ft(R(:,ii)<1)); 
       parfor jj=1:250
            FPU(ii,jj)=CalcFalsePositive(testName,DurT,TestF(ii),0);
       end
    end
    
    [~,FreqT(2),FreqT(3)]=Credible_Interval_High_Density(FreqT(1),TestF,0.95,'discrete',[1 14]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Risk of false positve
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    FPProb(1)=CalcFalsePositive(testName,DurT,FreqT(1),1); % Calacute the probability of a 
    [~,FPProb(2),FPProb(3)]=Credible_Interval_High_Density(FPProb(1),FPU(:),0.95,'continuous',[0 1]);
elseif(strcmp(Scenario,'Ag_Cut'))
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Determine the equivilant quarantine
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    Risk=1;
    load('TestingonExit_RTPCR_24hrDelay_AgCutoff=5.5934.mat','q','IDSLS','IDSLA');
    RRTPCR=Probability_Onward((1-pA).*IDSLS(q==7)+pA.*IDSLA(q==7),Risk);

    load(['TestingonExit_' testName '_NoDelay_AgCutoff=5.5934.mat'],'q','IDSLS','IDSLA')
    RAg=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    QDurE=min(q(RAg<=RRTPCR));

    load(['TestingonEntryExit_' testName '_NoDelay_AgCutoff=5.5934.mat'],'q','IDSLS','IDSLA')
    RAg=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    QDurEE=min(q(RAg<=RRTPCR));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Determine the minimal testing frequency
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    ft=[1:14];

    load(['Testing_Frequency_' testName '_AgCutoff=5.5934.mat'],'RTotA','RTotS')
    R=(1-pA).*RTotS+pA.*RTotA;
    FreqT=max(ft(R<1)); % take the maximum time between tests
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Risk of false positve
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    FPProb=CalcFalsePositive(testName,DurT,FreqT,1); % Calacute the probability of a 
elseif(strcmp(Scenario,'Baseline_Ag_Cut'))
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Determine the equivilant quarantine
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    Risk=1;
    load('TestingonExit_RTPCR_24hrDelay_DeltaVOC.mat','q','IDSLS','IDSLA');
    RRTPCR=Probability_Onward((1-pA).*IDSLS(q==7)+pA.*IDSLA(q==7),Risk);

    load(['TestingonExit_' testName '_NoDelay_DeltaVOC.mat'],'q','IDSLS','IDSLA')
    RAg=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    QDurE=min(q(RAg<=RRTPCR));

    load(['TestingonEntryExit_' testName '_NoDelay_DeltaVOC.mat'],'q','IDSLS','IDSLA')
    RAg=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    QDurEE=min(q(RAg<=RRTPCR));
    
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
    FPProb=CalcFalsePositive(testName,DurT,FreqT,1); % Calacute the probability of a 

elseif(strcmp(Scenario,'Ag_Cut_Long'))
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Determine the equivilant quarantine
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    Risk=1;
    load('TestingonExit_RTPCR_24hrDelay_AgCutoff=10.mat','q','IDSLS','IDSLA');
    RRTPCR=Probability_Onward((1-pA).*IDSLS(q==7)+pA.*IDSLA(q==7),Risk);

    load(['TestingonExit_' testName '_NoDelay_AgCutoff=10.mat'],'q','IDSLS','IDSLA')
    RAg=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    QDurE=min(q(RAg<=RRTPCR));

    load(['TestingonEntryExit_' testName '_NoDelay_AgCutoff=10.mat'],'q','IDSLS','IDSLA')
    RAg=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    QDurEE=min(q(RAg<=RRTPCR));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Determine the minimal testing frequency
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    ft=[1:14];

    load(['Testing_Frequency_' testName '_AgCutoff=10.mat'],'RTotA','RTotS')
    R=(1-pA).*RTotS+pA.*RTotA;
    FreqT=max(ft(R<1)); % take the maximum time between tests
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Risk of false positve
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    FPProb=CalcFalsePositive(testName,DurT,FreqT,1); % Calacute the probability of a 
end
end


