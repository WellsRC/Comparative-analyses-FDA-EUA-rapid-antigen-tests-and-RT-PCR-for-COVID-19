function [pA,IncubationI,R0,ts,td] = BaselineParameters
%BASELINEPARAMETERS 

pA=0.351;

R0=2.5;

ts=5.72;

td=ts+20;

IncubationI=integral(@(t)ViralShedding_Symptomatic(t,td),0,ts);
end

