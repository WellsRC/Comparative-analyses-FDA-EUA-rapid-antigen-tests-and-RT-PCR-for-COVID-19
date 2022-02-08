clear;
clc;

load('RAgTest_Name_Order.mat');
Test=cell(length(testName)+1,1);
PNA=cell(length(testName)+1,1);
Sensitivity_Text=cell(length(testName)+1,1);
Test{1}='RT-PCR';
PNA{1}='N/A';

[S,AG_TP,RTPCR_RTP,RTPCR_TP,TP] = Test_Specificity ('RT-PCR',1);

Stemp=zeros(1000,1);

parfor ii=1:1000
    [Stemp(ii),~,~,~,~] = Test_Specificity ('RT-PCR',0);
end

Sensitivity=S;
[~,Sensitivity_LB,Sensitivity_UB]=Credible_Interval_High_Density(Sensitivity,Stemp,0.95,'continuous',[0 1]);
Sensitivity_Text{1}=[num2str(100.*Sensitivity,'%4.2f') '% (' num2str(100.*Sensitivity_LB,'%4.2f') '%' char(8211) num2str(100.*Sensitivity_UB,'%4.2f') '%)'];
for tn=1:length(testName)
    Test{tn+1}=testName{tn};

    [S,AG_TP,RTPCR_RTP,RTPCR_TP,TP] = Test_Specificity (testName{tn},1);
    PNA{tn+1}=[num2str(AG_TP) '/' num2str(RTPCR_RTP)];

    Stemp=zeros(1000,1);

    parfor ii=1:1000
        [Stemp(ii),~,~,~,~] = Test_Specificity (testName{tn},0);
    end
    
    Sensitivity=S;
    
    [~,Sensitivity_LB,Sensitivity_UB]=Credible_Interval_High_Density(Sensitivity,Stemp,0.95,'continuous',[0 1]);
    Sensitivity_Text{tn+1}=[num2str(100.*Sensitivity,'%4.2f') '% (' num2str(100.*Sensitivity_LB,'%4.2f') '%' char(8211) num2str(100.*Sensitivity_UB,'%4.2f') '%)'];
end

T=table(Test,PNA,Sensitivity_Text);

writetable(T,'TableS4_Specificity.csv');