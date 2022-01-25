%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Produces the output for table 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

addpath([pwd '\Delta_Variant']);
addpath([pwd '\Delta_Variant\Results']);

load('RAgTest_Name_Order.mat');

Ntest=length(testName);

Quarantine_Duration_Exit=cell(Ntest,1);
Quarantine_Duration_Entry_Exit=cell(Ntest,1);
Serial_Testing_Freq=cell(Ntest,1);
False_Positive_Risk=cell(Ntest,1);



beta0=zeros(Ntest,1);
beta1=zeros(Ntest,1);

[pA,~,~,~,~] = BaselineParameters; % Proportion of infections being asymptomatic
DurT=14; % Duration to evalaute false positive

MLE_QEE=zeros(Ntest,1);
Unc_QEE=zeros(Ntest,1000);
Risk=1;
for ii=1:Ntest
    load('TestingonExit_RTPCR_24hrDelay_DeltaVOC.mat','q','IDSLS','IDSLA');
    RRTPCR=Probability_Onward((1-pA).*IDSLS(q==7)+pA.*IDSLA(q==7),Risk);

    load(['TestingonEntryExit_' testName{ii} '_NoDelay_DeltaVOC.mat'],'q','IDSLS','IDSLA')
    RAg=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    MLE_QEE(ii)=min(q(RAg<=RRTPCR));
    
    
    
    load('TestingonExit_RTPCR_24hrDelay_DeltaVOC_Uncertainty.mat','q','IDSLSv','IDSLAv')
    RRTPCR=Probability_Onward((1-pA).*IDSLSv(:,q==7)+pA.*IDSLAv(:,q==7),Risk);
        
    load(['TestingEntryExit_' testName{ii} '_NoDelay_DeltaVOC_Uncertainty.mat'],'q','IDSLSv','IDSLAv')
    R2Ag=Probability_Onward((1-pA).*IDSLSv+pA.*IDSLAv,Risk);
    
    for jj=1:length(RRTPCR)
        Unc_QEE(ii,jj)=min(q(R2Ag(jj,:)<=RRTPCR(jj)));
    end
end

[~,LBQEE,UBQEE]=Credible_Interval_High_Density(mode(Unc_QEE(:)),Unc_QEE(:),0.95,'discrete',[1 14]);
fprintf('Range of the quarantine duration for entry and exit for 18 tests: %3.0f to %3.0f (%3.0f - %3.0f) \n',[min(MLE_QEE) max(MLE_QEE) LBQEE UBQEE]);

for ii=1:Ntest
    [QDE,QDEE,STF,FPR] = Quarantine_Testing_Freq(testName{ii},pA,DurT,'Baseline');
    Quarantine_Duration_Exit{ii}=[num2str(QDE(1),'%d') ' (' num2str(QDE(2),'%d')  char(8211) num2str(QDE(3),'%d') ')'];
    Quarantine_Duration_Entry_Exit{ii}=[num2str(QDEE(1),'%d') ' (' num2str(QDEE(2),'%d')  char(8211) num2str(QDEE(3),'%d') ')'];
    Serial_Testing_Freq{ii}=[num2str(STF(1),'%d') ' (' num2str(STF(2),'%d')  char(8211) num2str(STF(3),'%d') ')'];
    False_Positive_Risk{ii}=[num2str(FPR(1),'%4.3g') ' (' num2str(FPR(2),'%4.3g')  char(8211) num2str(FPR(3),'%4.3g') ')'];
end



T=table(testName,Quarantine_Duration_Exit,Quarantine_Duration_Entry_Exit,Serial_Testing_Freq,False_Positive_Risk);

rmpath([pwd '\Delta_Variant']);
rmpath([pwd '\Delta_Variant\Results']);

writetable(T,'Table1.csv');



