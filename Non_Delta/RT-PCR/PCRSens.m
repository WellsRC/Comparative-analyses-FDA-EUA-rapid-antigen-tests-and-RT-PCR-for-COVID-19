function p = PCRSens(t,beta)

mmm=exp(beta(1)-beta(2)^2);
p=beta(3).*lognpdf(t,beta(1),beta(2))./lognpdf(mmm,beta(1),beta(2));
end

