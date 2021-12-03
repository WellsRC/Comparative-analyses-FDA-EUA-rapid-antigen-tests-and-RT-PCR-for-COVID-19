function R = InfectiousnessfromInfectionTestingFreqRTPCRFollowup(ttemp,utemp,timet,testtype,R0S,R0A,pA,ts,td,SelfIsolate,betaPCR)
%InfectiousnessfromInfection returns the infectiousness at time t given total virus
%shed and R0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% t - time post-infection
% u - time of testing
% R0S - Reproductive number symptomatic
% R0A - Reproductive number asymptomatic
% pA - proportion asymptomatic
% ts - time from infection to symptom onset (i.e. incubation period)
% tL - Duration of the latent period
% SelfIsolate - 0 no self-isolation of symptomatic otherwise self-isolation
% upon symptom onset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S - Probability of that it is a true positive based on the viral load

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u=utemp(:)';
t=ttemp(:)';
% Computation for asymptomatic
RA=R0A.*ViralShedding_Asymptomatic(t,td);

Senstd=inf; % Need to use infity as the positve will not depend on the truncated values

SA=zeros(length(timet),length(u));
SARTPCR=zeros(length(timet),length(u));
for jj=1:length(timet)
    SA(jj,:) = TestSensitivity(u+timet(jj),ts,Senstd,[testtype{jj}],betaPCR);
    SARTPCR(jj,:) = TestSensitivity(u+timet(jj),ts,Senstd,[],betaPCR);
end
% Computation for symptomatic
if(SelfIsolate==0)
    RS=R0S.*ViralShedding_Symptomatic(t,td);
else
    RS=zeros(size(t));
    RS(t<=ts)=R0S.*ViralShedding_Symptomatic(t(t<=ts),td);
end

SS=zeros(length(timet),length(u));
SSRTPCR=zeros(length(timet),length(u));
for jj=1:length(timet)
    SS(jj,:) = TestSensitivity(u+timet(jj),ts,Senstd,[testtype{jj}],betaPCR);
    SSRTPCR(jj,:) = TestSensitivity(u+timet(jj),ts,Senstd,[],betaPCR);
end

%Combine asymptmatic and symptomaic 
R=(1-pA).*RS.*prod(SS.*(1-SSRTPCR),1)+pA.*RA.*prod(SA.*(1-SARTPCR),1); % Used to only evalaute the RT-PCR follow-up and not the false-negative Rapid Ag test
R=R';
R=reshape(R,size(utemp));
end

