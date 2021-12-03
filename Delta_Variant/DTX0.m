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

testtype=cell(1,1);
testtype{1}=[];


parfor jj=1:14 
    IDSLS(jj)=((1./ts).*integral2(@(u,t)InfectiousnessfromInfectionTesting(t+u,u,[max(q(jj)-1,0)],testtype,R0S,R0A,0,ts,td,SelfIsolate,betaRTPCR),0,ts,q(jj),inf));
    IDSLA(jj)=((1./td).*integral2(@(u,t)InfectiousnessfromInfectionTesting(t+u,u,[max(q(jj)-1,0)],testtype,R0S,R0A,1,ts,td,0,betaRTPCR),0,td,q(jj),inf));  
end

save('TestingonExit_RTPCR_24hrDelay_DeltaVOC.mat');



betaRTPCRv=zeros(1000,3);
for jj=1:1000
    [betaRTPCRv(jj,:),~]=ParameterCOVIDTest([],0);
end
save('RTPCR_Parameter_Uncertainty.mat','betaRTPCRv');
testtype=cell(1,1);
testtype{1}=[];

IDSLSv=zeros(1000,length(IDSLS));
IDSLAv=zeros(1000,length(IDSLA));

for jj=1:14 
    parfor zz=1:1000
        IDSLSt(zz)=((1./ts).*integral2(@(u,t)InfectiousnessfromInfectionTesting(t+u,u,[max(q(jj)-1,0)],testtype,R0S,R0A,0,ts,td,SelfIsolate,betaRTPCRv(zz,:)),0,ts,q(jj),inf));
        IDSLAt(zz)=((1./td).*integral2(@(u,t)InfectiousnessfromInfectionTesting(t+u,u,[max(q(jj)-1,0)],testtype,R0S,R0A,1,ts,td,0,betaRTPCRv(zz,:)),0,td,q(jj),inf));  
    end
    IDSLSv(:,jj)=IDSLSt;
    IDSLAv(:,jj)=IDSLAt;
end

save('TestingonExit_RTPCR_24hrDelay_DeltaVOC_Uncertainty.mat','IDSLSv','IDSLAv','betaRTPCRv','R0','ts','td','q');
