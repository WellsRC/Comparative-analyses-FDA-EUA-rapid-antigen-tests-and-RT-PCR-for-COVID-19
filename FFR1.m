%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency testing: BD Veritor 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parpool(16); % Parallel pool
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEW Sensitivity curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delayRTPCR=1;
load('RAgTest_Name.mat','testName');
NumTests=length(testName);
for TestN=1:NumTests
    tL=2.9;
    [pA,~,R0,ts,td] = BaselineParameters(tL);

    SelfIsolate=1; %If sympmatics self-isolate
    R0S=R0;
    R0A=R0;
    RTotA=zeros(14,1);
    RTotS=zeros(14,1);
    RTotARF=zeros(14,1);
    RTotSRF=zeros(14,1);

    load('MLE-Estimate-RTPCR-Hill_Incubation_8_29_days.mat','beta')
    betaRTPCR=beta;
    load([testName{TestN} '_LR_Parameters.mat'],'beta');

    testtype=cell(14,1);
    timetoff=cell(14,1);
    NT=zeros(14,1);
    for dT=1:14
        timetoff{dT}=[0:dT:(floor(td))];
        NT(dT)=length([timetoff{dT}]);
        temp=cell(NT(dT),1);

        for ii=1:length([timetoff{dT}])
            temp{ii}=beta;
        end
        testtype{dT}=temp;
    end

    parfor dT=1:14  

        [RS,RA,RSRF,RARF] = SerialTestingRTPCRFollowup(testtype{dT},[timetoff{dT}],R0S,R0A,ts,tL,td,SelfIsolate,NT(dT),dT,betaRTPCR,delayRTPCR);

        RTotS(dT)=sum(RS);
        RTotA(dT)=sum(RA);
        
        
        RTotSRF(dT)=sum(RSRF);
        RTotARF(dT)=sum(RARF);
    end
    save(['Testing_Frequency_RTPCR_Follow_up_' testName{TestN} '_Hellewell.mat']);
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OLD Sensitivity curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delayRTPCR=1;
for TestN=1:NumTests
    tL=2.9;
    [pA,IncubationI,R0,ts,td] = BaselineParameters(tL);

    SelfIsolate=1; %If sympmatics self-isolate
    R0S=R0;
    R0A=R0;
    
    RTotA=zeros(14,1);
    RTotS=zeros(14,1);
    RTotARF=zeros(14,1);
    RTotSRF=zeros(14,1);

    
    load([testName{TestN} '_LR_Parameters.mat'],'beta');

    testtype=cell(14,1);
    timetoff=cell(14,1);
    NT=zeros(14,1);
    for dT=1:14
        timetoff{dT}=[0:dT:(floor(td))];
        NT(dT)=length([timetoff{dT}]);
        temp=cell(NT(dT),1);

        for ii=1:length([timetoff{dT}])
            temp{ii}=beta;
        end
        testtype{dT}=temp;
    end

    parfor dT=1:14  

        [RS,RA,RSRF,RARF] = SerialTestingRTPCRFollowupOLD(testtype{dT},[timetoff{dT}],R0S,R0A,ts,tL,td,SelfIsolate,NT(dT),dT,delayRTPCR);

        
        RTotS(dT)=sum(RS);
        RTotA(dT)=sum(RA);
        
        
        RTotSRF(dT)=sum(RSRF);
        RTotARF(dT)=sum(RARF);
    end
    save(['Testing_Frequency_RTPCR_Follow_up_' testName{TestN} '_NatComm.mat']);
end

    