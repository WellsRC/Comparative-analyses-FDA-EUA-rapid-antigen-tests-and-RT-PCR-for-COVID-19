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
    

    % Allcoate memory for output
    IDSLS=zeros(1,length(q)); 
    IDSLA=zeros(1,length(q)); 

    % Get Basline parameters
    [~,~,R0,ts,td] = BaselineParameters;
    AgCutoff=ts+10;
    IDSLS=IDSLS(:); % Vectorize the matrix
    IDSLA=IDSLA(:); % Vectorize the matrix

    R0S=R0; % Set R0 for symptomatic
    R0A=R0; % Set R0 for asymptomatic


    [betaRTPCR,betaAg]=ParameterCOVIDTest(testName{TestN},1);
    testtype=cell(1,1);
    testtype{1}=betaAg;

    parfor jj=1:14 
        IDSLS(jj)=((1./ts).*integral2(@(u,t)InfectiousnessfromInfectionTesting(t+u,u,[q(jj)],testtype,R0S,R0A,0,ts,td,AgCutoff,SelfIsolate,betaRTPCR),0,ts,q(jj),inf));
        IDSLA(jj)=((1./td).*integral2(@(u,t)InfectiousnessfromInfectionTesting(t+u,u,[q(jj)],testtype,R0S,R0A,1,ts,td,AgCutoff,0,betaRTPCR),0,td,q(jj),inf));  
    end

    save(['TestingonExit_' testName{TestN} '_NoDelay_DeltaVOC_AgCutoff.mat']);
end



load('RTPCR_Parameter_Uncertainty.mat','betaRTPCRv');

NSS=length(betaRTPCRv(:,1));

for TestN=1:NumTests
    q=[1:14]; % Quarantine durations consideredd
    
    SelfIsolate=1; % Self-isolation
    

    % Allcoate memory for output
    IDSLS=zeros(NSS,1); 
    IDSLA=zeros(NSS,1); 

    % Get Basline parameters
    [~,~,R0,ts,td] = BaselineParameters;
    AgCutoff=ts+10;

    IDSLS=IDSLS(:); % Vectorize the matrix
    IDSLA=IDSLA(:); % Vectorize the matrix
    
    IDSLSv=zeros(NSS,length(q)); % Vectorize the matrix
    IDSLAv=zeros(NSS,length(q)); % Vectorize the matrix
    
    R0S=R0; % Set R0 for symptomatic
    R0A=R0; % Set R0 for asymptomatic

    for jj=1:14 
        parfor ns=1:NSS
            [~,betaAg]=ParameterCOVIDTest(testName{TestN},0);
            testtype=cell(2,1);
            testtype{1}=betaAg;
            testtype{2}=betaAg;

            IDSLS(ns)=((1./ts).*integral2(@(u,t)InfectiousnessfromInfectionTesting(t+u,u,[q(jj)],testtype,R0S,R0A,0,ts,td,AgCutoff,SelfIsolate,betaRTPCRv(ns,:)),0,ts,q(jj),inf));
            IDSLA(ns)=((1./td).*integral2(@(u,t)InfectiousnessfromInfectionTesting(t+u,u,[q(jj)],testtype,R0S,R0A,1,ts,td,AgCutoff,0,betaRTPCRv(ns,:)),0,td,q(jj),inf)); 
        end
        IDSLSv(:,jj)=IDSLS;
        IDSLAv(:,jj)=IDSLA;
    end
    save(['TestingonExit_' testName{TestN} '_NoDelay_DeltaVOC_AgCutoff_Uncertainty.mat']);
end
