% Random entry in qiaratine with testing on exit
clear;

pobj=parpool(16); % Parallel pool
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RT-PCR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('RAgTest_Name.mat','testName');
NumTests=length(testName);
for TestN=1:NumTests
    q=[1:14]; % Quarantine durations consideredd
    
    SelfIsolate=1; % Self-isolation
    tL=[2.9]; % vecotor for the incbation periods to be integrated over

    % Allcoate memory for output
    IDSLS=zeros(1,length(q)); 
    IDSLA=zeros(1,length(q)); 

    % Get Basline parameters
    [pA,~,R0,ts] = BaselineParameters(tL); % Does not matter here what ts is fed in 

    IDSLS=IDSLS(:); % Vectorize the matrix
    IDSLA=IDSLA(:); % Vectorize the matrix

    td=ts+20; % Asymptomatic increase 20 days from symptom onset

    R0S=R0; % Set R0 for symptomatic
    R0A=R0; % Set R0 for asymptomatic


    load('MLE-Estimate-RTPCR-Hill_Incubation_8_29_days.mat','beta')
    betaRTPCR=beta;
    load([testName{TestN} '_LR_Parameters.mat'],'beta');
    testtype=cell(2,1);
    testtype{1}=beta;
    testtype{2}=beta;

    parfor jj=1:14 
        IDSLS(jj)=((1./ts).*integral2(@(u,t)InfectiousnessfromInfectionTesting(t+u,u,[0 q(jj)],testtype,R0S,R0A,0,ts,tL,td,SelfIsolate,betaRTPCR),0,ts,q(jj),inf));
        IDSLA(jj)=((1./td).*integral2(@(u,t)InfectiousnessfromInfectionTesting(t+u,u,[0 q(jj)],testtype,R0S,R0A,1,ts,tL,td,0,betaRTPCR),0,td,q(jj),inf));  
    end

    save(['TestingonEntryExit_' testName{TestN} '_NoDelay_Hellewell.mat']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RT-PCR (OLD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for TestN=1:NumTests    
    q=[1:14]; % Quarantine durations consideredd
    
    SelfIsolate=1; % Self-isolation
    tL=[2.9]; % vecotor for the incbation periods to be integrated over

    % Allcoate memory for output
    IDSLS=zeros(1,length(q)); 
    IDSLA=zeros(1,length(q)); 

    % Get Basline parameters
    [pA,~,R0,ts] = BaselineParameters(tL); % Does not matter here what ts is fed in 

    IDSLS=IDSLS(:); % Vectorize the matrix
    IDSLA=IDSLA(:); % Vectorize the matrix

    td=ts+20; % Asymptomatic increase 20 days from symptom onset

    R0S=R0; % Set R0 for symptomatic
    R0A=R0; % Set R0 for asymptomatic

    load([testName{TestN} '_LR_Parameters.mat'],'beta');
    testtype=cell(2,1);
    testtype{1}=beta;
    testtype{2}=beta;

    parfor jj=1:14 
        IDSLS(jj)=((1./ts).*integral2(@(u,t)InfectiousnessfromInfectionTestingOLD(t+u,u,[0 q(jj)],testtype,R0S,R0A,0,ts,tL,td,SelfIsolate),0,ts,q(jj),inf));
        IDSLA(jj)=((1./td).*integral2(@(u,t)InfectiousnessfromInfectionTestingOLD(t+u,u,[0 q(jj)],testtype,R0S,R0A,1,ts,tL,td,0),0,td,q(jj),inf));  
    end

    save(['TestingonEntryExit_' testName{TestN} '_NoDelay_NatComm.mat']);
end
clear;

