% Random entry in qiaratine with testing on exit
clear;

pobj=parpool(16); % Parallel pool
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RT-PCR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q=[1:14]; % Quarantine durations consideredd
SelfIsolate=1; % Self-isolation


% Allcoate memory for output
IDSLS=zeros(1,length(q)); 
IDSLA=zeros(1,length(q)); 

% Get Basline parameters
[~,~,R0,ts,td] = BaselineParameters;

IDSLS=IDSLS(:); % Vectorize the matrix
IDSLA=IDSLA(:); % Vectorize the matrix

R0S=R0; % Set R0 for symptomatic
R0A=R0; % Set R0 for asymptomatic


[betaRTPCR,~]=ParameterCOVIDTest([],1);
AgCutoffPSO=10;
testtype=cell(1,1);
testtype{1}=[];


parfor jj=1:14 
    IDSLS(jj)=((1./ts).*integral2(@(u,t)InfectiousnessfromInfectionTesting(t+u,u,[max(q(jj)-1,0)],testtype,AgCutoffPSO,R0S,R0A,0,ts,td,SelfIsolate,betaRTPCR),0,ts,q(jj),inf));
    IDSLA(jj)=((1./td).*integral2(@(u,t)InfectiousnessfromInfectionTesting(t+u,u,[max(q(jj)-1,0)],testtype,AgCutoffPSO,R0S,R0A,1,ts,td,0,betaRTPCR),0,td,q(jj),inf));  
end

save(['TestingonExit_RTPCR_24hrDelay_AgCutoff=' num2str(AgCutoffPSO) '.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RT-PCR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q=[1:14]; % Quarantine durations consideredd
SelfIsolate=1; % Self-isolation


% Allcoate memory for output
IDSLS=zeros(1,length(q)); 
IDSLA=zeros(1,length(q)); 

% Get Basline parameters
[~,~,R0,ts,td] = BaselineParameters;

IDSLS=IDSLS(:); % Vectorize the matrix
IDSLA=IDSLA(:); % Vectorize the matrix

R0S=R0; % Set R0 for symptomatic
R0A=R0; % Set R0 for asymptomatic


[betaRTPCR,~]=ParameterCOVIDTest([],1);
AgCutoffPSO=fminbnd(@(x)(integral(@(t)ViralShedding_Symptomatic(t,td),0,ts+x)-0.99).^2,0,10);
testtype=cell(1,1);
testtype{1}=[];


parfor jj=1:14 
    IDSLS(jj)=((1./ts).*integral2(@(u,t)InfectiousnessfromInfectionTesting(t+u,u,[max(q(jj)-1,0)],testtype,AgCutoffPSO,R0S,R0A,0,ts,td,SelfIsolate,betaRTPCR),0,ts,q(jj),inf));
    IDSLA(jj)=((1./td).*integral2(@(u,t)InfectiousnessfromInfectionTesting(t+u,u,[max(q(jj)-1,0)],testtype,AgCutoffPSO,R0S,R0A,1,ts,td,0,betaRTPCR),0,td,q(jj),inf));  
end

save(['TestingonExit_RTPCR_24hrDelay_AgCutoff=' num2str(AgCutoffPSO) '.mat']);