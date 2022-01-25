function [mle,lb,ub] = Credible_Interval_High_Density(mle,data,alphaZ,datatype,databnd)

if(strcmp(datatype,'continuous'))
    pd = fitdist(data(:),'Kernel','Kernel','epanechnikov','support',[databnd(1) databnd(2)]);
    
    if((mle>databnd(1)) && (databnd(2)>mle))
        options = optimset('TolX',10^(-16));
        [bnd]=fminbnd(@(bp)(pdf(pd,bp)-pdf(pd,icdf(pd,min(cdf(pd,bp)+alphaZ,1)))).^2+((cdf(pd,icdf(pd,min(cdf(pd,bp)+alphaZ,1)))-cdf(pd,bp))-alphaZ).^2,databnd(1),mle,options);

        lb=bnd;
        ub=icdf(pd,min(cdf(pd,bnd)+alphaZ,1));
    elseif (mle==databnd(1))
        lb=mle;
        ub=prctile(data,95);
    else        
        ub=mle;
        lb=prctile(data,5);
    end
elseif(strcmp(datatype,'discrete'))
    xc=[databnd(1):databnd(2)];
    y=histcounts(data,[(databnd(1)-0.5):(databnd(2)+0.5)]);
    y=y./sum(y);
    z=flip(unique(y));
    z=z(z>0);
    AUC=zeros(size(z));
    for kk=1:length(z)
       AUC(kk)=sum(y(y>=z(kk)))-alphaZ;
    end
    zmin=z(AUC==min(AUC(AUC>=0)));
    lb=min(xc(y>=zmin));
    ub=max(xc(y>=zmin));
end
end

