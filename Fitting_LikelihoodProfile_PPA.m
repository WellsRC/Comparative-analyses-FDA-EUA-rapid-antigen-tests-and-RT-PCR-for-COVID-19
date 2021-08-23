%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, txt, ~] =xlsread('Rapid_Ag_PPA.xlsx');
TT={txt{1,:}};

testName=cell(length(unique(TT))-1,1);
cc=1;
for ii=1:length(TT)
   if(~isempty(TT{ii}))
       testName{cc}=TT{ii};
       cc=cc+1;
   end
end
save('RAgTest_Name.mat','testName');
Ntest=length(unique(TT))-1;
T=readtable('Rapid_Ag_PPA.xlsx','Range','A2:BU34');
DateOS=T.DaysFromSymptomOnset;
for ii=1:Ntest
   if(ii>1)
       eval(sprintf('w=T.Weight_%d;', ii-1));
       eval(sprintf('truepos=T.TestPositive_%d;', ii-1));
       eval(sprintf('totalpos=T.RT_PCRPositive_%d;', ii-1));
   else
       w=T.Weight;
       truepos=T.TestPositive;
       totalpos=T.RT_PCRPositive;
   end
   Dt=DateOS(~isnan(w))';
   opts= optimset('MaxIter',10^6,'MaxFunEvals',10^6,'TolFun',10^(-16),'TolX',10^(-16),'UseParallel',false,'Display','off');
   [beta,MLE]=fmincon(@(z)GenFit(z,Dt,truepos(~isnan(w))',totalpos(~isnan(w))',w(~isnan(w))'),[1 -0.1],[],[],[],[],[0 -100],[100 0 ],[],opts);
   MLE=-MLE;
   save([testName{ii} '_LR_Parameters.mat'],'MLE','beta','Dt','totalpos','truepos','w');
end