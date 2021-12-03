close all;
clear;
clc;
load('Nasal_Swab_RTPCR-MinDayOnset-MaxDayOnset-SampleSiz-Positive.mat');

Tmin=RTPCRA(:,1)';
Tmax=RTPCRA(:,2)';
wtemp=1./(Tmax-Tmin+1);
w=[];
T=[];
P=[];
N=[];
Ptemp=RTPCRA(:,4)';
Ntemp=RTPCRA(:,3)';
for ii=1:length(Tmin)
   for jj=Tmin(ii):Tmax(ii)
      w=[w;wtemp(ii)];
      T=[T;jj];
      P=[P;Ptemp(ii)];
      N=[N;Ntemp(ii)];
   end
end

TU=unique(T);
PRD=TU;
for ii=1:length(TU)
    ff=find(T==TU(ii));
    PRD(ii)=sum(w(ff).*P(ff))./sum(w(ff).*N(ff));
end

scatter(TU,PRD,20,'k','filled');
hold on;


RTPCRmodelv=[1 2 3 5 6 7];
t=linspace(-20,80,100001);

CM=[hex2rgb('#e41a1c')
hex2rgb('#377eb8')
hex2rgb('#4daf4a')
hex2rgb('#984ea3')
hex2rgb('#ff7f00')
hex2rgb('#a65628')];

for rr=1:6
    [S,mm] = RTPCRSensitivity(t,RTPCRmodelv(rr));
    plot(t,S,'-','color',CM(rr,:),'LineWidth',2);
end
xlim([-20 80])
ylim([0 1])
ll=legend({'Pooled data','Weibull','Gamma','Skew-normal','log-Normal','Transformed Gamma','Transformed Weibull'});
xlabel('Days post-symptom onset')
ylabel('Diagnostic sensitivity');
set(gca,'LineWidth',2,'tickdir','out','XTick',[-20:10:80],'Ytick',[0:0.1:1],'Xminortick','on','Yminortick','on','Fontsize',14)
% print(gcf,'Nasal_Sens_RT_PCR','-dpng','-r300');
% xlim([-5 0])
% ll.Location='NorthWest';
% set(gca,'LineWidth',2,'tickdir','out','XTick',[-5:1:0],'Ytick',[0:0.1:1],'Xminortick','on','Yminortick','on','Fontsize',14)
print(gcf,'Nasal_Sens_RT_PCR_Presymptomatic','-dpng','-r300');