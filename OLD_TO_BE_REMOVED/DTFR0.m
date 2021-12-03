%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency testing: RTPCR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parpool(6); % Parallel pool
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEW Sensitivity curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tL=2.9;
[pA,~,R0,ts,td] = BaselineParameters(tL);

SelfIsolate=1; %If sympmatics self-isolate
R0S=R0;
R0A=R0;
RTotA=zeros(14,1);
RTotS=zeros(14,1);

load('MLE-Estimate-RTPCR-Hill_Incubation_8_29_days.mat','beta')
betaRTPCR=beta;

testtype=cell(14,1);
timetoff=cell(14,1);
NT=zeros(14,1);
for dT=1:14
    timetoff{dT}=[0:dT:(floor(td))];
    NT(dT)=length([timetoff{dT}]);
    temp=cell(NT(dT),1);
    testtype{dT}=temp;
end

for delayTR=0:5
    parfor dT=1:14  

        [RS,RA] = SerialTestingDelay(testtype{dT},[timetoff{dT}],R0S,R0A,ts,tL,td,SelfIsolate,NT(dT),dT,betaRTPCR,delayTR);

        RTotS(dT)=sum(RS);
        RTotA(dT)=sum(RA);
    end
    save([num2str(delayTR) '-day_Delay_Testing_Frequency_RTPCR_Hellewell.mat']);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OLD Sensitivity curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tL=2.9;
[pA,IncubationI,R0,ts,td] = BaselineParameters(tL);

SelfIsolate=1; %If sympmatics self-isolate
R0S=R0;
R0A=R0;
RTotA=zeros(14,1);
RTotS=zeros(14,1);

testtype=cell(14,1);
timetoff=cell(14,1);
NT=zeros(14,1);
for dT=1:14
    timetoff{dT}=[0:dT:(floor(td))];
    NT(dT)=length([timetoff{dT}]);
    temp=cell(NT(dT),1);
    testtype{dT}=temp;
end
for delayTR=0:5
    parfor dT=1:14  
        [RS,RA] = SerialTestingDelayOLD(testtype{dT},[timetoff{dT}],R0S,R0A,ts,tL,td,SelfIsolate,NT(dT),dT,delayTR);

        RTotS(dT)=sum(RS);
        RTotA(dT)=sum(RA);
    end
    save([num2str(delayTR) '-day_Delay_Testing_Frequency_RTPCR_NatComm.mat']);
end