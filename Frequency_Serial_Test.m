function [pA,R0v,FS] = Frequency_Serial_Test(testName)

load(['Testing_Frequency_' testName '_DeltaVOC.mat'],'RTotA','RTotS','R0')

R0v=linspace(1.2,5,101);
pA=[0.1:0.001:0.95];

FS=zeros(length(pA),length(R0v));
for jj=1:length(pA)
    for ii=1:length(R0v)
        findx=find((pA(jj).*RTotA+(1-pA(jj)).*RTotS).*R0v(ii)/R0<1);
        if(~isempty(findx))
           FS(jj,ii)=findx(end); 
        end
    end
end

end

