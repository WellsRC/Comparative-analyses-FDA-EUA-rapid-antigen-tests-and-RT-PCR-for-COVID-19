function L = LikelihoodLAMPCurve(x,Time,Positive,TotalTest,tL)
beta=10.^x(end-1:end);

p=PCRSens(Time,beta,tL);

Avg=integral(@(x)PCRSens(x,beta,tL),0,40)./40;
L2=log((p.^Positive).*((1-p).^(TotalTest-Positive)));
L1=log((Avg.^5).*(1-Avg).^(12-5));

L=-sum(L2)-L1;
end

