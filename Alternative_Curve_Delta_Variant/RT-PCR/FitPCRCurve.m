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

LB=[UB-30 6 -6];
UB2=[UB 8 log10(1)];

Gx =[4.230765931945100   0.959892250476516];

options = optimoptions('ga','PlotFcn', @gaplotbestf,'MaxGenerations',5*10^3,'FunctionTolerance',10^(-16),'UseParallel',true);

[part,fvalt]=ga(@(x)LikelihoodPCRCurve([x],PtID,dstart,dlast,TPtID,TDate,TResult,Gx),length(LB),[],[],[],[],LB,UB2,[],options);

opts= optimset('MaxIter',10^4,'MaxFunEvals',10^4,'TolFun',10^(-16),'TolX',10^(-16),'UseParallel',false,'Display','off');
[par,MLE]=fmincon(@(x)LikelihoodPCRCurve(x,PtID,dstart,dlast,TPtID,TDate,TResult,Gx),part,[],[],[],[],LB,UB2,[],opts);

beta=10.^par(end-1:end);
if(fvalt<MLE)
    beta=10.^part(end-1:end);
    MLE=fvalt;
    par=part;
end


betat=1./((Gx(1)-1).*Gx(2));
beta=[betat beta];
save('MLE-Estimate-RTPCR_Alternative.mat');

parv=par(1:end-2);
Pk=log10(linspace(10^(-1),1,1001));
b1=linspace(6,9,1001);

[Pk,b1]=meshgrid(Pk,b1);
Pk=Pk(:);
b1=b1(:);
L=zeros(length(Pk),1);
parfor ii=1:1002001
    L(ii)=-LikelihoodPCRCurve([parv b1(ii) Pk(ii)],PtID,dstart,dlast,TPtID,TDate,TResult,Gx);
end

b1=10.^b1;
Pk=10.^Pk;
betaRTU=[1./((Gx(1)-1).*Gx(2)).*ones(size(b1)) b1 Pk];
betaRTU=betaRTU(~isinf(L),:);
L=L(~isinf(L));
save('Uncertainty-BetaEstimate-RTPCR_Alternative.mat');