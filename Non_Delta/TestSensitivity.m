function S = TestSensitivity(t,ts,testtype,beta)
%SensitivityvsViralLoad(V,asym) - Returns the sensitivity for a given viral
%load

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ts- incubation period
% testtype - test type

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S - Probability of that it is a true positive based on the viral load

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Gx =[6.127 3.277];
mm=Gx(1).*((Gx(2)-1)/Gx(2))^(1/Gx(2));

 [S]=PCRSens(t,beta);
 V=ViralShedding_Symptomatic(t,inf); % using inf as we need to construc the mapping
if(~isempty(testtype))
    tt=[mm:0.1:46.2]; % made coarsered to improve the search tt=[ts:0.001:90]; 
    Vx=ViralShedding_Symptomatic(tt,inf); % using inf as we need to construc the mapping
    PPA=LR(tt-ts,testtype);
    if(~isempty(t(t<mm)))
        S(t<mm)=pchip(Vx,PPA,V(t<mm)).*PCRSens(t(t<mm),beta);
    end
    S(t>=mm)=PCRSens(t(t>=mm),beta).*LR(t(t>=mm)-ts,testtype);
end

end

