 clear;
rng default
%https://cmmid.github.io/topics/covid19/pcr-positivity-over-time.html
load('HCW-Test-Positive.mat');

PtID=unique([Symptom{:,1}]);
dlast=zeros(size(PtID));
dstart=zeros(size(PtID));
UB=zeros(size(PtID));


TPtID=[Test{:,1}];
TDate=[datenum({Test{:,2}})];
temp={Test{:,3}};
TResult=double(strcmp(temp,'TRUE'));
temp={Test{:,5}};
for ii=1:length(temp)
   if(strcmp(temp{ii},'NA'))
       temp{ii}=-1;
   end
end
Ct=cell2mat(temp);


for ii=1:length(PtID)
    temp=[Symptom{:,1}];
    f=find(temp==PtID(ii));
    temp={Symptom{f,3}};
    g=strcmp(temp,'TRUE');
    h=find(g==1,1);
    dstart(ii)=datenum({Symptom{f(h),2}});
    dlast(ii)=datenum({Symptom{f(h)-1,2}});
    f=find(PtID(ii)==TPtID);
    temp=TResult(f);
    g=find(temp==1,1);
    temp=TDate(f(g));
    
    temp2=Ct(f);
    g=find(temp2>0,1);
    temp2=TDate(f(g));
    if(~isempty(temp)&&~isempty(temp2))
        UB(ii)=min(min(temp,temp2),dstart(ii));    
    elseif(~isempty(temp))
        UB(ii)=min(temp,dstart(ii));    
    elseif(~isempty(temp2))
        UB(ii)=min(temp2,dstart(ii));
    else
        UB(ii)=dstart(ii);
    end
end

LB=[UB-30 log10(0.01) -6];
UB2=[UB log10(10) log10(1)];

Gx =[4.230765931945100   0.959892250476516];
% ENSURE DISTRIBUTION IS SET CORRECTLY
options = optimoptions('ga','PlotFcn', @gaplotbestf,'MaxGenerations',5*10^3,'FunctionTolerance',10^(-16),'UseParallel',true);

[part,fvalt]=ga(@(x)LikelihoodPCRCurve_Gamma([x],PtID,dstart,dlast,TPtID,TDate,TResult,Gx),length(LB),[],[],[],[],LB,UB2,[],options);

opts= optimset('MaxIter',10^4,'MaxFunEvals',10^4,'TolFun',10^(-16),'TolX',10^(-16),'UseParallel',false,'Display','off');
[par,MLE]=fmincon(@(x)LikelihoodPCRCurve_Gamma(x,PtID,dstart,dlast,TPtID,TDate,TResult,Gx),part,[],[],[],[],LB,UB2,[],opts);

beta=10.^par(end-1:end);
if(fvalt<MLE)
    beta=10.^part(end-1:end);
    MLE=fvalt;
    par=part;
end


betat=log((Gx(1)-1).*Gx(2))+beta(1)^2;
beta=[betat beta];
save('MLE-Estimate-RTPCR-Gamma.mat');
