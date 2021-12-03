function F = GenFit(x,T,P,N,wp)
L=zeros(size(wp));

L(wp>=0)=((LR(T(wp>=0),x).^P(wp>=0)).*((1-LR(T(wp>=0),x)).^(N(wp>=0)-P(wp>=0)))).^(wp(wp>=0));


funkown=find(wp<0);
for ii=1:length(funkown)
    Pvec=[P(funkown(ii)):N(funkown(ii))];
    L(funkown(ii))=sum((LR(T(funkown(ii)),x).^Pvec).*((1-LR(T(funkown(ii)),x)).^(N(funkown(ii))-Pvec)));
end


L(L==0)=10^(-300);


F=-sum(log(L));
end

