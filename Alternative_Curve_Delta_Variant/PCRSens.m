function p = PCRSens(t,beta)

mmm=1./beta(1);
t(t<0)=0;
p=beta(3).*tpdf(log(t.*beta(1)),beta(2))./tpdf(log(mmm.*beta(1)),beta(2));
end

