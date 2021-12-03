function [pA,R0v,QD] = Duration_Quarantine_Test(testName)

load(['TestingonExit_' testName '_NoDelay_Hellewell.mat'],'IDSLA','IDSLS','R0','q')
[pAb,~,~,~,~] = BaselineParameters;
R0v=linspace(1.2,5,101);
pA=[0.1:0.001:0.95];

QD=zeros(length(pA),length(R0v));
for jj=1:length(pA)
    for ii=1:length(R0v)
        findx=find((pA(jj).*IDSLA+(1-pA(jj)).*IDSLS).*R0v(ii)/R0<=(pAb.*IDSLA(q==7)+(1-pAb).*IDSLS(q==7)));
        if(~isempty(findx))
           QD(jj,ii)=min(findx); 
        end
    end
end
end

