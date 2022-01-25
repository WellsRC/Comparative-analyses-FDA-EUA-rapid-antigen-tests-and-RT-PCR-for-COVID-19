function L = LikelihoodPCRCurve_Gamma(x,PtID,dstart,dlast,TPtID,TDate,TResult,Gx)
beta=10.^x(end-1:end);

betat=log((Gx(1)-1).*Gx(2))+beta(1)^2;
beta=[betat beta];
TI=x(1:end-2);
[~,Lstart] = DistIncubationG(dstart-TI);
[~,Llast] = DistIncubationG(dlast-TI);
L1=log(Lstart-Llast);

[~,b]=ismember(TPtID,PtID);

p=PCRSens(TDate'-TI(b(b>0)),beta);
L2=log((p.^TResult).*((1-p).^(1-TResult)));

L= -sum(L1)-sum(L2);
end

