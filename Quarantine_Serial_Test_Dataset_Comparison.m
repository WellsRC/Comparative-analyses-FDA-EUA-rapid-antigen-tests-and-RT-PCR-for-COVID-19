clear;
clc;
close all;

addpath([pwd '\Delta_Variant']);
addpath([pwd '\Delta_Variant\Results']);


[pA,~,~,~,~] = BaselineParameters; % Proportion of infections being asymptomatic

load('RAgTest_Name_DatasetCompare.mat');
testNamev=testName;
for ii=1:length(testNamev(:,1))
    testNameFDA=testNamev{ii,1};
    testNameComm=testNamev{ii,2};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Determine the equivilant quarantine
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    Risk=1;
    load('TestingonExit_RTPCR_24hrDelay_DeltaVOC.mat','q','IDSLS','IDSLA');
    RRTPCR=Probability_Onward((1-pA).*IDSLS(q==7)+pA.*IDSLA(q==7),Risk);

    load(['TestingonExit_' testNameFDA '_NoDelay_DeltaVOC.mat'],'q','IDSLS','IDSLA')
    RAg=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    QDurE=min(q(RAg<=RRTPCR));

    load(['TestingonEntryExit_' testNameFDA '_NoDelay_DeltaVOC.mat'],'q','IDSLS','IDSLA')
    RAg=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    QDurEE=min(q(RAg<=RRTPCR));
    
    
    
    load('TestingonExit_RTPCR_24hrDelay_DeltaVOC_Uncertainty.mat','q','IDSLSv','IDSLAv')
    RRTPCR=Probability_Onward((1-pA).*IDSLSv(:,q==7)+pA.*IDSLAv(:,q==7),Risk);
    
    load(['TestingonExit_' testNameFDA '_NoDelay_DeltaVOC_Uncertainty.mat'],'q','IDSLSv','IDSLAv')
    RAg=Probability_Onward((1-pA).*IDSLSv+pA.*IDSLAv,Risk);
    
    load(['TestingEntryExit_' testNameFDA '_NoDelay_DeltaVOC_Uncertainty.mat'],'q','IDSLSv','IDSLAv')
    R2Ag=Probability_Onward((1-pA).*IDSLSv+pA.*IDSLAv,Risk);
    
    TestE=zeros(length(RRTPCR),1);
    TestEE=zeros(length(RRTPCR),1);
    for jj=1:length(RRTPCR)
        TestE(jj)=min(q(RAg(jj,:)<=RRTPCR(jj)));
        TestEE(jj)=min(q(R2Ag(jj,:)<=RRTPCR(jj)));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Determine the minimal testing frequency
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    ft=[1:14];

    load(['Testing_Frequency_' testNameFDA '_DeltaVOC.mat'],'RTotA','RTotS')
    R=(1-pA).*RTotS+pA.*RTotA;
    FreqT=max(ft(R<1)); % take the maximum time between tests
    
    load(['Testing_Frequency_' testNameFDA '_DeltaVOC_Uncertainty.mat'],'RTotAv','RTotSv')
    R=(1-pA).*RTotSv+pA.*RTotAv;
    TestF=zeros(length(R(1,:)),1);
    for jj=1:length(TestF)
       TestF(jj)=max(ft(R(:,jj)<1)); 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Determine the equivilant quarantine
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    Risk=1;
    load('TestingonExit_RTPCR_24hrDelay_DeltaVOC.mat','q','IDSLS','IDSLA');
    RRTPCR=Probability_Onward((1-pA).*IDSLS(q==7)+pA.*IDSLA(q==7),Risk);

    load(['TestingonExit_' testNameComm '_NoDelay_DeltaVOC.mat'],'q','IDSLS','IDSLA')
    RAg=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    QDurE_Comm=min(q(RAg<=RRTPCR));

    load(['TestingonEntryExit_' testNameComm '_NoDelay_DeltaVOC.mat'],'q','IDSLS','IDSLA')
    RAg=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    QDurEE_Comm=min(q(RAg<=RRTPCR));
    
    
    
    load('TestingonExit_RTPCR_24hrDelay_DeltaVOC_Uncertainty.mat','q','IDSLSv','IDSLAv')
    RRTPCR=Probability_Onward((1-pA).*IDSLSv(:,q==7)+pA.*IDSLAv(:,q==7),Risk);
    
    load(['TestingonExit_' testNameComm '_NoDelay_DeltaVOC_Uncertainty.mat'],'q','IDSLSv','IDSLAv')
    RAg=Probability_Onward((1-pA).*IDSLSv+pA.*IDSLAv,Risk);
    
    load(['TestingEntryExit_' testNameComm '_NoDelay_DeltaVOC_Uncertainty.mat'],'q','IDSLSv','IDSLAv')
    R2Ag=Probability_Onward((1-pA).*IDSLSv+pA.*IDSLAv,Risk);
    
    TestE_Comm=zeros(length(RRTPCR),1);
    TestEE_Comm=zeros(length(RRTPCR),1);
    for jj=1:length(RRTPCR)
        TestE_Comm(jj)=min(q(RAg(jj,:)<=RRTPCR(jj)));
        TestEE_Comm(jj)=min(q(R2Ag(jj,:)<=RRTPCR(jj)));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Determine the minimal testing frequency
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    ft=[1:14];

    load(['Testing_Frequency_' testNameComm '_DeltaVOC.mat'],'RTotA','RTotS')
    R=(1-pA).*RTotS+pA.*RTotA;
    FreqT_Comm=max(ft(R<1)); % take the maximum time between tests
    
    load(['Testing_Frequency_' testNameComm '_DeltaVOC_Uncertainty.mat'],'RTotAv','RTotSv')
    R=(1-pA).*RTotSv+pA.*RTotAv;
    TestF_Comm=zeros(length(R(1,:)),1);
    for jj=1:length(TestF_Comm)
       TestF_Comm(jj)=max(ft(R(:,jj)<1)); 
    end
    
    % Cannot to boot strapping as then they would have different RT-PCR
    % curves
    dQE=QDurE_Comm-QDurE;
    dQEu=TestE_Comm-TestE;
    [~,LBQE,UBQE]=Credible_Interval_High_Density(dQE,dQEu,0.95,'discrete',[-3 3]);
    
    dQEE=QDurEE_Comm-QDurEE;
    dQEEu=TestEE_Comm-TestEE;
    [~,LBQEE,UBQEE]=Credible_Interval_High_Density(dQEE,dQEEu,0.95,'discrete',[-3 3]);
    
    dST=FreqT_Comm-FreqT;
    dSTu=TestF_Comm-TestF;
    [~,LBST,UBST]=Credible_Interval_High_Density(dST,dSTu,0.95,'discrete',[-3 3]);
    
    fprintf('========================================================================================================== \n');
    fprintf([ testNameFDA '\n']);
    fprintf('========================================================================================================== \n');
    fprintf(['Change in the duration of the length of quarantine for test on exit for ' testNameFDA ': %2.0f (%2.0f - %2.0f) \n'],[dQE LBQE UBQE]);
    fprintf(['Change in the duration of the length of quarantine for test on entry and exit for ' testNameFDA ': %2.0f (%2.0f - %2.0f) \n'],[dQEE LBQEE UBQEE]);
    fprintf(['Change in the time between test during serial tesing for ' testNameFDA ': %2.0f (%2.0f - %2.0f) \n'],[dST LBST UBST]);
end

rmpath([pwd '\Delta_Variant']);
rmpath([pwd '\Delta_Variant\Results']);
