%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Produces the output for table 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

addpath([pwd '\Non_Delta']);
addpath([pwd '\Non_Delta\Results']);

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

for ii=1:Ntest
    [QDE,QDEE,STF,FPR] = Quarantine_Testing_Freq(testName{ii},pA,DurT,'Non_Delta');
    Quarantine_Duration_Exit{ii}=[num2str(QDE(1),'%d') ' (' num2str(QDE(2),'%d')  char(8211) num2str(QDE(3),'%d') ')'];
    Quarantine_Duration_Entry_Exit{ii}=[num2str(QDEE(1),'%d') ' (' num2str(QDEE(2),'%d')  char(8211) num2str(QDEE(3),'%d') ')'];
    Serial_Testing_Freq{ii}=[num2str(STF(1),'%d') ' (' num2str(STF(2),'%d')  char(8211) num2str(STF(3),'%d') ')'];
    False_Positive_Risk{ii}=[num2str(FPR(1),'%4.3g') ' (' num2str(FPR(2),'%4.3g')  char(8211) num2str(FPR(3),'%4.3g') ')'];
end



T=table(testName,Quarantine_Duration_Exit,Quarantine_Duration_Entry_Exit,Serial_Testing_Freq,False_Positive_Risk);

rmpath([pwd '\Non_Delta']);
rmpath([pwd '\Non_Delta\Results']);

writetable(T,'Table1_Non_Delta.csv');



