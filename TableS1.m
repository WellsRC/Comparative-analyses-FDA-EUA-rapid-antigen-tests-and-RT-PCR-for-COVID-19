%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Produces the output for table 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

load('RAgTest_Name_Order.mat');

Ntest=length(testName);

Quarantine_Duration_Exit=zeros(Ntest,1);
Quarantine_Duration_Entry_Exit=zeros(Ntest,1);
Serial_Testing_Freq=zeros(Ntest,1);
False_Positive_Risk=zeros(Ntest,1);

pA=0.308; % Proportion of infections being asymptomatic
DurT=14; % Duration to evalaute false positive

for ii=1:Ntest
    [Quarantine_Duration_Exit(ii),Quarantine_Duration_Entry_Exit(ii),Serial_Testing_Freq(ii),False_Positive_Risk(ii)] = Quarantine_Testing_Freq_NatComm(testName{ii},pA,DurT);
end
False_Positive_Risk=round(False_Positive_Risk,3); % Round for the presentation in the table;
T=table(testName,Quarantine_Duration_Exit,Quarantine_Duration_Entry_Exit,Serial_Testing_Freq,False_Positive_Risk);

writetable(T,'TableS1.csv');