function [ProbFalsePosSerial] = CalcFalsePositive(testName,DurT,FreqT)

PFP=Test_Specificity(testName); % Obtain the specificity

ProbFalsePosSerial=0;
temp=[1:FreqT:(FreqT*(DurT-1)+1)];
for ii=1:FreqT 
   tempNT=length(temp((temp>=1+14*(ii-1)) & (temp<=14*ii)));
   ProbFalsePosSerial=ProbFalsePosSerial+1-PFP^(tempNT);
end
ProbFalsePosSerial=ProbFalsePosSerial./FreqT;

end

