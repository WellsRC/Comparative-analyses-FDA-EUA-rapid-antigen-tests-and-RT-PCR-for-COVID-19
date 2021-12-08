function [RS,RA,RSRF,RARF] = SerialTestingRTPCRFollowup(testtype,timetoff,R0S,R0A,ts,td,SelfIsolate,NT,dT,betaPCR,delayRTPCR)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


RA=zeros(length(timetoff)+1,1);
RS=zeros(length(timetoff)+1,1);

RARF=zeros(length(timetoff)+1,1);
RSRF=zeros(length(timetoff)+1,1);

RS(1)=integral2(@(u,tx)InfectiousnessfromInfectionTestingFreq(tx,u,[],testtype,R0S,R0A,0,ts,td,SelfIsolate,betaPCR),0,dT,0,@(u)(u))./dT;
RA(1)=integral2(@(u,tx)InfectiousnessfromInfectionTestingFreq(tx,u,[],testtype,R0S,R0A,1,ts,td,0,betaPCR),0,dT,0,@(u)(u))./dT;

parfor ii=0:(NT-2)
    dayT=[timetoff(1:(ii+1))]; % have to add +1 to ii as ii starts at zero
    % have to add 2 to ii for RS and RA as the ii starts at zero and
    % need index to start at 2
    RS(ii+2)=integral2(@(u,tx)InfectiousnessfromInfectionTestingFreq(tx,u,dayT,testtype,R0S,R0A,0,ts,td,SelfIsolate,betaPCR),0,dT,@(u)(u+dT.*(ii)),@(u)(u+dT.*(ii+1)))./dT;
    RA(ii+2)=integral2(@(u,tx)InfectiousnessfromInfectionTestingFreq(tx,u,dayT,testtype,R0S,R0A,1,ts,td,0,betaPCR),0,dT,@(u)(u+dT.*(ii)),@(u)(u+dT.*(ii+1)))./dT;
    
    
    RSRF(ii+2)=integral2(@(u,tx)InfectiousnessfromInfectionTestingFreqRTPCRFollowup(tx,u,dayT,testtype,R0S,R0A,0,ts,td,SelfIsolate,betaPCR),0,dT,@(u)min((u+dT.*(ii)+delayRTPCR),(u+dT.*(ii+1))),@(u)(u+dT.*(ii+1)))./dT;
    RARF(ii+2)=integral2(@(u,tx)InfectiousnessfromInfectionTestingFreqRTPCRFollowup(tx,u,dayT,testtype,R0S,R0A,1,ts,td,0,betaPCR),0,dT,@(u)min((u+dT.*(ii)+delayRTPCR),(u+dT.*(ii+1))),@(u)(u+dT.*(ii+1)))./dT;
end

dayT=[timetoff];
RS(end)=integral2(@(u,tx)InfectiousnessfromInfectionTestingFreq(tx,u,dayT,testtype,R0S,R0A,0,ts,td,SelfIsolate,betaPCR),0,dT,@(u)(u+dT.*(NT-1)),td)./dT;
RA(end)=integral2(@(u,tx)InfectiousnessfromInfectionTestingFreq(tx,u,dayT,testtype,R0S,R0A,1,ts,td,0,betaPCR),0,dT,@(u)(u+dT.*(NT-1)),td)./dT;

RSRF(end)=integral2(@(u,tx)InfectiousnessfromInfectionTestingFreqRTPCRFollowup(tx,u,dayT,testtype,R0S,R0A,0,ts,td,SelfIsolate,betaPCR),0,dT,@(u)min((u+dT.*(NT-1)+delayRTPCR),td),td)./dT;
RARF(end)=integral2(@(u,tx)InfectiousnessfromInfectionTestingFreqRTPCRFollowup(tx,u,dayT,testtype,R0S,R0A,1,ts,td,0,betaPCR),0,dT,@(u)min((u+dT.*(NT-1)+delayRTPCR),td),td)./dT;
        
end

