function R = MissedTesting(u,timet,pA,ts,tL,td,testtype)
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



SA=zeros(length(timet),length(u));
for jj=1:length(timet)
    SA(jj,:) = SensitivityvsViralLoad(ViralShedding_Asymptomatic(u+timet(jj),tL),ts,tL,td,[testtype{jj}]);
    %SA(t-u<0)=0; % Should not occur as we feed in t+u but have condition just in case
end
% Computation for symptomatic


SS=zeros(length(timet),length(u));
for jj=1:length(timet)
    SS(jj,:) = SensitivityvsViralLoad(ViralShedding_Symptomatic(u+timet(jj),tL),ts,tL,td,[testtype{jj}]);
    
    SS(jj,u+timet(jj)>=ts)=0;
end


%Combine asymptmatic and symptomaic 
R=(1-pA).*prod((1-SS),1)+pA.*prod((1-SA),1);

test=0;
end

