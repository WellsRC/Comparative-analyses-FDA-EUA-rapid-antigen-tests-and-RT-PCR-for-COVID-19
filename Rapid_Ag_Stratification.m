%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Stratify the type of test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
clear;
clc;
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

Ntest=length(unique(TT))-1;

for ii=1:Ntest
    
    load([testName{ii} '_LR_Parameters.mat'],'beta');
    PPA=LR([0 20 40],beta); % Evalaute at 20 days after symptom onset to 
    if(1-PPA(3)/PPA(1)<0.01)
       fprintf([testName{ii} ': Constant \n']); 
    elseif(PPA(2)<0.01)
       fprintf([testName{ii} ': Rapid decline \n']);        
    else
       fprintf([testName{ii} ': Gradual decline \n']);            
    end
end
