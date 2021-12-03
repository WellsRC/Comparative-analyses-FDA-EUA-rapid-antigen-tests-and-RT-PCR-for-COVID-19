%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Produces the output for table 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

load('RAgTest_Name_Order.mat');

Ntest=length(testName);

Quarantine_Duration_Exit=zeros(Ntest,1);
Quarantine_Duration_Entry_Exit=zeros(Ntest,1);
Serial_Testing_Freq=zeros(Ntest,1);
False_Positive_Risk=zeros(Ntest,1);



beta0=zeros(Ntest,1);
beta1=zeros(Ntest,1);

pA=0.308; % Proportion of infections being asymptomatic
DurT=14; % Duration to evalaute false positive

for ii=1:Ntest
    [Quarantine_Duration_Exit(ii),Quarantine_Duration_Entry_Exit(ii),Serial_Testing_Freq(ii),False_Positive_Risk(ii)] = Quarantine_Testing_Freq_Hellewell(testName{ii},pA,DurT);
    
    load([testName{ii} '_LR_Parameters.mat'],'beta');
    beta0(ii)=beta(1);
    beta1(ii)=beta(2);
end
False_Positive_Risk=round(False_Positive_Risk,3); % Round for the presentation in the table;
T=table(testName,Quarantine_Duration_Exit,Quarantine_Duration_Entry_Exit,Serial_Testing_Freq,False_Positive_Risk);
T2=table(testName,beta0,beta1);
writetable(T,'Table1.csv');
writetable(T2,'Table_PPA_Coefficients.csv');
