clear;
close all;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% % Hellewell et al (RTPCR)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
ts=4.4;

figure('units','normalized','outerposition',[0.1 0.1 0.6 0.6]);
subplot('Position',[0.095070422535211,0.177867863572038,0.881161971830986,0.746456460752287]);
load('MLE-Estimate-RTPCR.mat','beta','TPtID','TDate','TResult','PtID','par');

betaRTPCR=beta;
TI=par(1:end-2);
[~,b]=ismember(TPtID,PtID);
dt=round(TDate'-TI(b(b>0)));
t=linspace(0,40,1001);

[CCtestRTPCR,~]= ColourTests('RTPCR');
p = PCRSens(t,betaRTPCR);
plot(t,p,'color',CCtestRTPCR,'LineWidth',2)
hold on;

L=[1:40];
S=zeros(1,40);

for ii=1:40
   S(ii)=sum(TResult(dt<=(ii+0.5)&dt>(ii-0.5)))./length(dt(dt<=(ii+0.5)&dt>(ii-0.5))); 
end
scatter(L,S,40,CCtestRTPCR,'filled')

legend({'Estimated','Avg Test Result:Binned'},'Fontsize',20,'location','NorthEast');
box off;
set(gca,'LineWidth',2,'tickdir','out','Fontsize',22,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:0.1:1],'Ylim',[0 1]);
xlabel('Days post-infection','Fontsize',26);
ylabel({'Diagnostic sensitivty'},'Fontsize',26);
title('RT-PCR','Fontsize',26);
print(gcf,['RTPCR_Sensitivity.png'],'-dpng','-r600');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Rapid antigen tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

load('RAgTest_Name_Order.mat')
Ntest=length(testName);
for ii=1:Ntest
    SIFigure_Test(testName{ii},ts,testName)
end