%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Produces the output for table 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
clear;
clc;
addpath([pwd '\Ag_Cut']);
addpath([pwd '\Ag_Cut\Results']);

load('RAgTest_Name_Order.mat');

Ntest=length(testName);

Quarantine_Duration_Exit=zeros(Ntest,1);
Quarantine_Duration_Entry_Exit=zeros(Ntest,1);
Serial_Testing_Freq=zeros(Ntest,1);
False_Positive_Risk=zeros(Ntest,1);


[pA,~,~,~,~] = BaselineParameters; % Proportion of infections being asymptomatic
DurT=14; % Duration to evalaute false positive

for ii=1:Ntest
    [Quarantine_Duration_Exit(ii),Quarantine_Duration_Entry_Exit(ii),Serial_Testing_Freq(ii),False_Positive_Risk(ii)] = Quarantine_Testing_Freq(testName{ii},pA,DurT,'Ag_Cut');
end
False_Positive_Risk=round(False_Positive_Risk,3); % Round for the presentation in the table;
T1=table(testName,Quarantine_Duration_Exit,Quarantine_Duration_Entry_Exit,Serial_Testing_Freq,False_Positive_Risk);


Quarantine_Duration_Exit=zeros(Ntest,1);
Quarantine_Duration_Entry_Exit=zeros(Ntest,1);
Serial_Testing_Freq=zeros(Ntest,1);
False_Positive_Risk=zeros(Ntest,1);


[pA,~,~,~,~] = BaselineParameters; % Proportion of infections being asymptomatic
DurT=14; % Duration to evalaute false positive

for ii=1:Ntest
    [Quarantine_Duration_Exit(ii),Quarantine_Duration_Entry_Exit(ii),Serial_Testing_Freq(ii),False_Positive_Risk(ii)] = Quarantine_Testing_Freq(testName{ii},pA,DurT,'Ag_Cut_Long');
end
False_Positive_Risk=round(False_Positive_Risk,3); % Round for the presentation in the table;
T2=table(testName,Quarantine_Duration_Exit,Quarantine_Duration_Entry_Exit,Serial_Testing_Freq,False_Positive_Risk);

rmpath([pwd '\Ag_Cut']);
rmpath([pwd '\Ag_Cut\Results']);

addpath([pwd '\Delta_Variant']);
addpath([pwd '\Delta_Variant\Results']);

load('RAgTest_Name_Order.mat');

Ntest=length(testName);

Quarantine_Duration_Exit=zeros(Ntest,1);
Quarantine_Duration_Entry_Exit=zeros(Ntest,1);
Serial_Testing_Freq=zeros(Ntest,1);
False_Positive_Risk=zeros(Ntest,1);



beta0=zeros(Ntest,1);
beta1=zeros(Ntest,1);

[pA,~,~,~,~] = BaselineParameters; % Proportion of infections being asymptomatic
DurT=14; % Duration to evalaute false positive

for ii=1:Ntest
    [Quarantine_Duration_Exit(ii),Quarantine_Duration_Entry_Exit(ii),Serial_Testing_Freq(ii),False_Positive_Risk(ii)] = Quarantine_Testing_Freq(testName{ii},pA,DurT,'Baseline_Ag_Cut');
end
False_Positive_Risk=round(False_Positive_Risk,3); % Round for the presentation in the table;
T3=table(testName,Quarantine_Duration_Exit,Quarantine_Duration_Entry_Exit,Serial_Testing_Freq,False_Positive_Risk);

rmpath([pwd '\Delta_Variant']);
rmpath([pwd '\Delta_Variant\Results']);



writetable(T1,'Table_Ag_Cut.xlsx','Sheet','Ag_Cut_6day');
writetable(T2,'Table_Ag_Cut.xlsx','Sheet','Ag_Cut_10day');
writetable(T3,'Table_Ag_Cut.xlsx','Sheet','Baseline');



