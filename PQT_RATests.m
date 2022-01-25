function [q,RRTPCR,RAgX,RAgEX] = PQT_RATests(testName,Scenario,pA)

Risk=1;

load(['TestingonExit_RTPCR_24hrDelay_' Scenario '.mat'],'IDSLA','IDSLS','q')
RRTPCR=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);
load(['TestingonExit_' testName '_NoDelay_' Scenario '.mat'],'IDSLA','IDSLS');
RAgX=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

load(['TestingonEntryExit_' testName '_NoDelay_' Scenario '.mat'],'IDSLA','IDSLS');
RAgEX=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

end

