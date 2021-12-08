function [RS,RA] = SerialTesting(testtype,timetoff,R0S,R0A,ts,td,AgCutoff,SelfIsolate,NT,dT,betaPCR)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


RA=zeros(length(timetoff)+1,1);
RS=zeros(length(timetoff)+1,1);

RS(1)=integral2(@(u,tx)InfectiousnessfromInfectionTestingFreq(tx,u,[],testtype,R0S,R0A,0,ts,td,AgCutoff,SelfIsolate,betaPCR),0,dT,0,@(u)(u))./dT;
RA(1)=integral2(@(u,tx)InfectiousnessfromInfectionTestingFreq(tx,u,[],testtype,R0S,R0A,1,ts,td,AgCutoff,0,betaPCR),0,dT,0,@(u)(u))./dT;

parfor ii=0:(NT-2)
    dayT=[timetoff(1:(ii+1))]; % have to add +1 to ii as ii starts at zero
    % have to add 2 to ii for RS and RA as the ii starts at zero and
    % need index to start at 2
    RS(ii+2)=integral2(@(u,tx)InfectiousnessfromInfectionTestingFreq(tx,u,dayT,testtype,R0S,R0A,0,ts,td,AgCutoff,SelfIsolate,betaPCR),0,dT,@(u)(u+dT.*(ii)),@(u)(u+dT.*(ii+1)))./dT;
    RA(ii+2)=integral2(@(u,tx)InfectiousnessfromInfectionTestingFreq(tx,u,dayT,testtype,R0S,R0A,1,ts,td,AgCutoff,0,betaPCR),0,dT,@(u)(u+dT.*(ii)),@(u)(u+dT.*(ii+1)))./dT;
end

dayT=[timetoff];
RS(end)=integral2(@(u,tx)InfectiousnessfromInfectionTestingFreq(tx,u,dayT,testtype,R0S,R0A,0,ts,td,AgCutoff,SelfIsolate,betaPCR),0,dT,@(u)(u+dT.*(NT-1)),td)./dT;
RA(end)=integral2(@(u,tx)InfectiousnessfromInfectionTestingFreq(tx,u,dayT,testtype,R0S,R0A,1,ts,td,AgCutoff,0,betaPCR),0,dT,@(u)(u+dT.*(NT-1)),td)./dT;
        
end

