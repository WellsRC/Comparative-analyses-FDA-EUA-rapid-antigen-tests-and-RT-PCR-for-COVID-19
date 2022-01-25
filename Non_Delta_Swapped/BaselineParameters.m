function [pA,IncubationI,R0,ts,td] = BaselineParameters
%BASELINEPARAMETERS 

pA=0.351;

R0=3.2;

ts=5.723;

td=ts+20;

IncubationI=integral(@(t)ViralShedding_Symptomatic(t,td),0,ts);
end

