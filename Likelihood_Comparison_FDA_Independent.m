clear;
clc;
addpath([pwd '\PPA_RA_Tests']);
LFDA=zeros(3,1);
LBoth=zeros(3,1);
LCommunity=zeros(3,1);

load('BinaxNOW (FDA)_LR_Parameters.mat','MLE')
LFDA(1)=MLE;
load('BinaxNOW_LR_Parameters.mat','MLE')
LBoth(1)=MLE;
load('BinaxNOW (Community)_LR_Parameters.mat','MLE')
LCommunity(1)=MLE;

load('CareStart (Anterior Nasal Swab - FDA)_LR_Parameters.mat','MLE')
LFDA(2)=MLE;
load('CareStart (Anterior Nasal Swab)_LR_Parameters.mat','MLE')
LBoth(2)=MLE;
load('CareStart (Anterior Nasal Swab - External)_LR_Parameters.mat','MLE')
LCommunity(2)=MLE;

load('Sofia (FDA)_LR_Parameters.mat','MLE')
LFDA(3)=MLE;
load('Sofia_LR_Parameters.mat','MLE')
LBoth(3)=MLE;
load('Sofia (CDC)_LR_Parameters.mat','MLE')
LCommunity(3)=MLE;


TestName={'BinaxNOW';'CareStart (AS)';'Sofia'};

LCombine=LFDA+LCommunity;

LDiff=LCombine-LBoth;

rmpath([pwd '\PPA_RA_Tests']);

SDiff=cell(3,1);

SDcut=chi2inv(0.95,2)./2; % Cut-off to determine if these curves/datasets are significantly different

SDiff(LDiff<=SDcut)={'No'};
SDiff(LDiff>SDcut)={'Yes'};

LFDA=round(LFDA,2);
LCommunity=round(LCommunity,2);
LBoth=round(LBoth,2);
LDiff=round(LDiff,2);
T=table(TestName,LFDA,LCommunity,LBoth,LDiff,SDiff);

writetable(T,'Compare_DataSets.csv');