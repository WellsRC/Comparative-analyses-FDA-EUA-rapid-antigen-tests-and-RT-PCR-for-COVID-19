% Random entry in qiaratine with testing on exit
clear;

pobj=parpool(20); % Parallel pool
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RT-PCR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qt=[0:14]; % Quarantine durations consideredd
durTt=[1:30];
SelfIsolate=1; % Self-isolation
tL=[2.9]; % vecotor for the incbation periods to be integrated over

% Allcoate memory for output
IDSL=zeros(length(durTt),length(qt)); % Assumes asymptomatic enter over infections period [0,ts+21]

% Get Basline parameters
[pA,~,R0,ts] = BaselineParameters(tL); % Does not matter here what ts is fed in 

[q,durT]=meshgrid(qt,durTt); % Create a mesh grid of the paramters being changes
q=q(:); % Vectorize the matrix
durT=durT(:); % Vectorize the matrix
IDSL=IDSL(:); % Vectorize the matrix

td=ts+20; % Asymptomatic increase 20 days from symptom onset

R0S=R0; % Set R0 for symptomatic
R0A=R0; % Set R0 for asymptomatic

load('MLE-Estimate-RTPCR-Hill.mat','beta');
testtype=cell(1,1);
testtype{1}=[];

parfor jj=1:450 
    IDSL(jj)=(1-pA).*((1./ts).*integral2(@(u,t)InfectiousnessfromInfectionTesting(t+u,u,[],testtype,R0S,R0A,0,ts,tL,td,SelfIsolate,beta),0,ts,q(jj),q(jj)+durT(jj)))+pA.*((1./td).*integral2(@(u,t)InfectiousnessfromInfectionTesting(t+u,u,[],testtype,R0S,R0A,1,ts,tL,td,0,beta),0,td,q(jj),q(jj)+durT(jj)));  
end

save('NoTest_Quaratnine.mat');

clear;