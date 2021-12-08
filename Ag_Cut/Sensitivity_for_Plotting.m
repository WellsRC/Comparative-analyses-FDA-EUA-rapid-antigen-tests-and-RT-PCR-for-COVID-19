function [MLE_RTPCR,MLE_Ag,U_RTPCR,U_Ag,MLE_PPA,U_PPA,Dt,totalpos,truepos,w,t,CCtestRTPCR,CCtest,SymP] = Sensitivity_for_Plotting(testName,ts)

[betaRTPCR,betaAg]=ParameterCOVIDTest(testName,1);
NS=1000;
t=linspace(0,40,1001);

MLE_Ag = TestSensitivity(t,ts,inf,betaAg,betaRTPCR);     
    
MLE_RTPCR = TestSensitivity(t,ts,inf,[],betaRTPCR);

MLE_PPA=LR(t-ts,betaAg);
    
U_Ag=zeros(NS,1001);
U_RTPCR=zeros(NS,1001);
U_PPA=zeros(NS,1001);
parfor ii=1:NS
    [betaRTPCR,betaAg]=ParameterCOVIDTest(testName,0);
    U_Ag(ii,:) = TestSensitivity(t,ts,inf,betaAg,betaRTPCR);   
    U_RTPCR(ii,:) = TestSensitivity(t,ts,inf,[],betaRTPCR);
    U_PPA(ii,:)=LR(t,betaAg);
end


load([testName '_LR_Parameters.mat'],'Dt','totalpos','truepos','w')


[CCtestRTPCR,~]=ColourTests('RTPCR');
[CCtest,SymP]= ColourTests(testName); %
end

