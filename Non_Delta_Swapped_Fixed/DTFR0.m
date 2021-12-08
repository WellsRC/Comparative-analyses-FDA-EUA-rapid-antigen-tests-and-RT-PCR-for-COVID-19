%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency testing: RTPCR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parpool(16); % Parallel pool
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEW Sensitivity curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,~,R0,ts,td] = BaselineParameters;

SelfIsolate=1; %If sympmatics self-isolate
R0S=R0;
R0A=R0;
RTotA=zeros(14,1);
RTotS=zeros(14,1);


[betaRTPCR,~]=ParameterCOVIDTest([],1);

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

        [RS,RA] = SerialTestingDelay(testtype{dT},[timetoff{dT}],R0S,R0A,ts,td,SelfIsolate,NT(dT),dT,betaRTPCR,delayTR);

        RTotS(dT)=sum(RS);
        RTotA(dT)=sum(RA);
    end
    save([num2str(delayTR) '-day_Delay_Testing_Frequency_RTPCR_General_Swapped.mat']);
end
