clear;
close all;
figure('units','normalized','outerposition',[0 0 1 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Baseline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
subplot('Position',[0.061491596638655,0.592256332320174,0.424600840336134,0.390000000000001]);
addpath([pwd '\Delta_Variant']);
addpath([pwd '\Delta_Variant\Results']);
[~,~,~,ts,~] = BaselineParameters;
ts1=ts;
load('MLE-Estimate-RTPCR.mat','TPtID','TDate','TResult','PtID','par','beta');
TI=par(1:end-2);
[~,b]=ismember(TPtID,PtID);
dt=round(TDate'-TI(b(b>0)));
t=linspace(0,50,1001);

S1 = TestSensitivity(t,ts,[],beta);

load('RTPCR_Parameter_Uncertainty.mat','betaRTPCRv');
SU=zeros(length(betaRTPCRv(:,1)),length(t));
for ii=1:length(betaRTPCRv(:,1))
   SU(ii,:)=  TestSensitivity(t,ts,[],betaRTPCRv(ii,:));
end

patch([t flip(t)],[prctile(SU,2.5) flip(prctile(SU,97.5))],'k','LineStyle','none','Facealpha',0.2); hold on
LL1=plot(t,S1,'k','LineWidth',2);
plot(ts.*ones(101,1),linspace(0,1,101),':','color',[0.75 0.75 0.75],'LineWidth',2);

L=[1:40];
S=zeros(1,40);

for ii=1:40
   S(ii)=sum(TResult(dt<=(ii+0.5)&dt>(ii-0.5)))./length(dt(dt<=(ii+0.5)&dt>(ii-0.5))); 
end
LL2=scatter(L,S,40,'k','filled');
legend([LL2,LL1],{'Avg Test Result:Binned','Estimated RT-PCR'},'Fontsize',20);

box off;
set(gca,'LineWidth',2,'tickdir','out','XTick',[0:5:50],'Xminortick','on','YTick',[0:0.1:1],'YminorTick','on','Fontsize',20);
ylabel('Diagnostic sensitivity','Fontsize',22);
xlabel('Days post-infection','Fontsize',22);

text(-6.30407911001236,0.995,'A','Fontsize',32,'FontWeight','bold');

rmpath([pwd '\Delta_Variant']);
rmpath([pwd '\Delta_Variant\Results']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Non-Delta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
subplot('Position',[0.561491596638655,0.592256332320174,0.424600840336134,0.390000000000001]);
addpath([pwd '\Non_Delta']);
addpath([pwd '\Non_Delta\Results']);
[~,~,~,ts,~] = BaselineParameters;
load('MLE-Estimate-RTPCR_Non_Delta.mat','TPtID','TDate','TResult','PtID','par','beta');
TI=par(1:end-2);
[~,b]=ismember(TPtID,PtID);
dt=round(TDate'-TI(b(b>0)));
t=linspace(0,50,1001);
ts2=ts;
S2 = TestSensitivity(t,ts,[],beta);

load('RTPCR_Parameter_Uncertainty.mat','betaRTPCRv');
SU=zeros(length(betaRTPCRv(:,1)),length(t));
for ii=1:length(betaRTPCRv(:,1))
   SU(ii,:)=  TestSensitivity(t,ts,[],betaRTPCRv(ii,:));
end

patch([t flip(t)],[prctile(SU,2.5) flip(prctile(SU,97.5))],'k','LineStyle','none','Facealpha',0.2); hold on
LL1=plot(t,S2,'k','LineWidth',2);
plot(ts.*ones(101,1),linspace(0,1,101),':','color',[0.75 0.75 0.75],'LineWidth',2);

L=[1:40];
S=zeros(1,40);

for ii=1:40
   S(ii)=sum(TResult(dt<=(ii+0.5)&dt>(ii-0.5)))./length(dt(dt<=(ii+0.5)&dt>(ii-0.5))); 
end
LL2=scatter(L,S,40,'k','filled');
legend([LL2,LL1],{'Avg Test Result:Binned','Estimated RT-PCR'},'Fontsize',20);

box off;
set(gca,'LineWidth',2,'tickdir','out','XTick',[0:5:50],'Xminortick','on','YTick',[0:0.1:1],'YminorTick','on','Fontsize',20);
ylabel('Diagnostic sensitivity','Fontsize',22);
xlabel('Days post-infection','Fontsize',22);

text(-6.30407911001236,0.995,'B','Fontsize',32,'FontWeight','bold');

rmpath([pwd '\Non_Delta']);
rmpath([pwd '\Non_Delta\Results']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Alternative curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
subplot('Position',[0.061491596638655,0.114,0.424600840336134,0.390000000000001]);
addpath([pwd '\Alternative_Curve_Delta_Variant']);
addpath([pwd '\Alternative_Curve_Delta_Variant\Results']);
[~,~,~,ts,~] = BaselineParameters;
load('MLE-Estimate-RTPCR_Alternative.mat','TPtID','TDate','TResult','PtID','par','beta');
TI=par(1:end-2);
[~,b]=ismember(TPtID,PtID);
dt=round(TDate'-TI(b(b>0)));
t=linspace(0,50,1001);
ts3=ts;
S3 = TestSensitivity(t,ts,[],beta);

load('RTPCR_Parameter_Uncertainty.mat','betaRTPCRv');
SU=zeros(length(betaRTPCRv(:,1)),length(t));
for ii=1:length(betaRTPCRv(:,1))
   SU(ii,:)=  TestSensitivity(t,ts,[],betaRTPCRv(ii,:));
end

patch([t flip(t)],[prctile(SU,2.5) flip(prctile(SU,97.5))],'k','LineStyle','none','Facealpha',0.2); hold on
LL1=plot(t,S3,'k','LineWidth',2);
plot(ts.*ones(101,1),linspace(0,1,101),':','color',[0.75 0.75 0.75],'LineWidth',2);

L=[1:40];
S=zeros(1,40);

for ii=1:40
   S(ii)=sum(TResult(dt<=(ii+0.5)&dt>(ii-0.5)))./length(dt(dt<=(ii+0.5)&dt>(ii-0.5))); 
end
LL2=scatter(L,S,40,'k','filled');
legend([LL2,LL1],{'Avg Test Result:Binned','Estimated RT-PCR'},'Fontsize',20);

box off;
set(gca,'LineWidth',2,'tickdir','out','XTick',[0:5:50],'Xminortick','on','YTick',[0:0.1:1],'YminorTick','on','Fontsize',20);
ylabel('Diagnostic sensitivity','Fontsize',22);
xlabel('Days post-infection','Fontsize',22);

text(-6.30407911001236,0.995,'C','Fontsize',32,'FontWeight','bold');

rmpath([pwd '\Alternative_Curve_Delta_Variant']);
rmpath([pwd '\Alternative_Curve_Delta_Variant\Results']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% All together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
subplot('Position',[0.561491596638655,0.114,0.424600840336134,0.390000000000001]);

plot(t-ts1,S1,'k',t-ts2,S2,'r',t-ts3,S3,'b','LineWidth',2); hold on
plot(zeros(101,1),linspace(0,1,101),':','color',[0.75 0.75 0.75],'LineWidth',2);
legend({'RT-PCR: log-Normal (4.4 day incubation period)','RT-PCR: log-Normal (5.72 day incubation period)','RT-PCR: log-Students t (4.4 day incubation period)'},'Fontsize',19);
legend boxoff;
box off;
xlim([-5.5 40]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[-5:5:50],'Xminortick','on','YTick',[0:0.1:1],'YminorTick','on','Fontsize',20);
ylabel('Diagnostic sensitivity','Fontsize',22);
xlabel('Days since symptom onset','Fontsize',22);

text(-11.2533992,0.995,'D','Fontsize',32,'FontWeight','bold');
print(gcf,['RT-PCR_Curves.png'],'-dpng','-r600');