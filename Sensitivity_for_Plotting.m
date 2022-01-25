function [MLE_RTPCR,MLE_Ag,U_RTPCR,U_Ag,MLE_PPA,U_PPA,Dt,totalpos,truepos,w,t,CCtestRTPCR,SymPRTPCR,CCtest,SymP] = Sensitivity_for_Plotting(testName,ts)

[betaRTPCR,betaAg]=ParameterCOVIDTest(testName,1);
NS=1000;
t=linspace(0,40,1001);

MLE_Ag = TestSensitivity(t,ts,betaAg,betaRTPCR);     
    
MLE_RTPCR = TestSensitivity(t,ts,[],betaRTPCR);

MLE_PPA=100.*LR(t,betaAg); % time zero is the time of symptom onset, which is what we want to be plotting from. 
    
U_Ag=zeros(NS,1001);
U_RTPCR=zeros(NS,1001);
U_PPA=zeros(NS,1001);
parfor ii=1:NS
    [betaRTPCR,betaAg]=ParameterCOVIDTest(testName,0);
    U_Ag(ii,:) = TestSensitivity(t,ts,betaAg,betaRTPCR);   
    U_RTPCR(ii,:) = TestSensitivity(t,ts,[],betaRTPCR);
    U_PPA(ii,:)=100.*LR(t,betaAg);
end

U_Ag(isnan(U_Ag))=0;
load([testName '_LR_Parameters.mat'],'Dt','totalpos','truepos','w')


[CCtestRTPCR,SymPRTPCR]=ColourTests('RTPCR');
[CCtest,SymP]= ColourTests(testName); 
end

