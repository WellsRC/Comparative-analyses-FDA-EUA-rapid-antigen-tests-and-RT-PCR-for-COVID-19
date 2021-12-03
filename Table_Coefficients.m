%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Produces the output for table 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

load('RAgTest_Name_Order.mat');

NSS=1000;
Ntest=length(testName);

MLE_beta0=zeros(Ntest,1);
MLE_beta1=zeros(Ntest,1);
beta0=zeros(Ntest,NSS);
beta1=zeros(Ntest,NSS);
b0=zeros(1,NSS);
b1=zeros(1,NSS);
for ii=1:Ntest
    [~,betaAg]=ParameterCOVIDTest(testName{ii},1);
    MLE_beta0(ii)=betaAg(1);
    MLE_beta1(ii)=betaAg(2);
    parfor jj=1:NSS        
        [~,betaAg]=ParameterCOVIDTest(testName{ii},0);
        b0(jj)=betaAg(1);
        b1(jj)=betaAg(2);        
    end
    beta0(ii,:)=b0;
    beta1(ii,:)=b1;
end
low95_beta0=prctile(beta0,2.5);
low95_beta1=prctile(beta1,2.5);

up95_beta0=prctile(beta0,97.5);
up95_beta1=prctile(beta1,97.5);

T2=table(testName,MLE_beta0,low95_beta0,up95_beta0,MLE_beta1,low95_beta1,up95_beta1);
writetable(T2,'Table_PPA_Coefficients.csv');
