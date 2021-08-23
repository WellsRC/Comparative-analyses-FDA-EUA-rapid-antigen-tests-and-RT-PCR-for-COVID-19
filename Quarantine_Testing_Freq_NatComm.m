function [QDurE,QDurEE,FreqT,FPProb] = Quarantine_Testing_Freq_NatComm(testName,pA,DurT)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Determine the equivilant quarantine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
Risk=1;
load('TestingonExit_RTPCR_24hrDelay_NatComm.mat','q','IDSLS','IDSLA');
RRTPCR=Probability_Onward((1-pA).*IDSLS(q==7)+pA.*IDSLA(q==7),Risk);

load(['TestingonExit_' testName '_NoDelay_NatComm.mat'],'q','IDSLS','IDSLA')
RAg=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

QDurE=min(q(RAg<RRTPCR));

load(['TestingonEntryExit_' testName '_NoDelay_NatComm.mat'],'q','IDSLS','IDSLA')
RAg=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

QDurEE=min(q(RAg<RRTPCR));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Determine the minimal testing frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
ft=[1:14];

load(['Testing_Frequency_' testName '_NatComm.mat'],'RTotA','RTotS')
R=(1-pA).*RTotS+pA.*RTotA;
FreqT=max(ft(R<1)); % take the maximum time between tests


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Risk of false positve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
FPProb=CalcFalsePositive(testName,DurT,FreqT); % Calacute the probability of a 
end

