%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Plots the probability of PQT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
close all;
% Colur blind palete
% CTest=[hex2rgb('#C6D4E1');hex2rgb('#0F2080');hex2rgb('#F5793A');hex2rgb('#A95AA1');hex2rgb('#85C0F9')];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%% Hellewell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55

CTest=[hex2rgb('#231B12');hex2rgb('#375E97');hex2rgb('#486824');hex2rgb('#F9A603');hex2rgb('#CF3721');hex2rgb('#810f7c')];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Exit 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Risk=1;
R=zeros(14,6);
load('TestingonExit_RTPCR_24hrDelay_Hellewell.mat')
R(:,1)=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

load('TestingonExit_LumiraDX (Anterior Nasal Swab)_NoDelay_Hellewell.mat')
R(:,2)=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

load('TestingonExit_Sofia (FDA)_NoDelay_Hellewell.mat')
R(:,3)=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

load('TestingonExit_BinaxNOW (FDA)_NoDelay_Hellewell.mat')
R(:,4)=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

load('TestingonExit_BD Veritor_NoDelay_Hellewell.mat')
R(:,5)=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

load('TestingonExit_CareStart (Anterior Nasal Swab - FDA)_NoDelay_Hellewell.mat')
R(:,6)=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);


load('RAgTest_PlotOrder.mat');
% testName=OrderPlotTest;
NumTest=length(testName);
RAllAgTest=zeros(14,NumTest);

for kk=1:NumTest
    load(['TestingonExit_' testName{kk} '_NoDelay_Hellewell.mat'],'IDSLS','IDSLA');
    RAllAgTest(:,kk)=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);
end

AgBetter=zeros(14,1);

for qq=1:14
    fAgBetter=find(RAllAgTest(qq,:)<=R(qq,1));
    AgBetter(qq)=length(fAgBetter);
end


figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.18114406779661./2,0.555,0.79885593220339./2,0.39]);

MF={'-p','-.s','-.d','-.^','-.o','-.h'};
for ii=1:6
   plot(q,sqrt(R(:,ii)),MF{ii},'LineWidth',2,'MarkerFacecolor',CTest(ii,:),'color',CTest(ii,:),'Markersize',10); hold on;
end
box off;
grid on;
% legend({'RT-PCR','LumiraDx','Sofia','BinaxNOW','BD Veritor'},'Fontsize',18,'NumColumns',3,'Location','NorthEast');
% legend boxoff;
yytick=[0 0.001 0.01 0.025 0.05 0.075 0.1:0.05:0.3];
set(gca,'LineWidth',2,'Tickdir','out','Fontsize',20,'XTick',[1:14],'Xlim',[1 14],'Xminortick','off','YTick',sqrt(yytick),'ylim',[0 sqrt(yytick(end))],'YTickLabel',num2str(yytick'),'Yminortick','on');
% xlabel('Duration of quarantine (days)','Fontsize',26);
ylabel({'Probability of','post-quarantine transmission'},'Fontsize',26);
title('Single rapid antigen test on exit','Fontsize',26);
text(-1.89,0.585,'A','Fontsize',34,'FontWeight','bold');

subplot('Position',[0.18114406779661./2,0.105,0.79885593220339./2,0.39]);


   plot(q,AgBetter,'k-o','LineWidth',2,'MarkerFacecolor','k','color','k','Markersize',10); hold on;

box off;
grid on;
%%%legend({'LumiraDx','Sofia','BinaxNOW','BD Veritor'},'Fontsize',18,'NumColumns',3,'Location','NorthEast');
% % % legend boxoff;

YTickLabelText=cell(NumTest+1,1);
for kk=0:NumTest
    if(rem(kk,2)==0)
        YTickLabelText{kk+1}=[ num2str(kk) '/' num2str(NumTest)];
    else
        YTickLabelText{kk+1}=' ';
    end
end

set(gca,'LineWidth',2,'Tickdir','out','Fontsize',20,'XTick',[1:14],'Xlim',[1 14],'YTick',[0:NumTest],'Yticklabel',{YTickLabelText{:}},'Yminortick','off','Xminortick','off');
xlabel('Duration of quarantine (days)','Fontsize',26);
ylim([0 NumTest]);
ylabel({'Fraction of antigen tests','outperforming RT-PCR'},'Fontsize',24);
text(-1.89,18.275,'C','Fontsize',34,'FontWeight','bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Entry and exit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Risk=1;
R=zeros(14,6);
load('TestingonExit_RTPCR_24hrDelay_Hellewell.mat')
R(:,1)=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

load('TestingonEntryExit_LumiraDX (Anterior Nasal Swab)_NoDelay_Hellewell.mat')
R(:,2)=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

load('TestingonEntryExit_Sofia (FDA)_NoDelay_Hellewell.mat')
R(:,3)=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

load('TestingonEntryExit_BinaxNOW (FDA)_NoDelay_Hellewell.mat')
R(:,4)=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

load('TestingonEntryExit_BD Veritor_NoDelay_Hellewell.mat')
R(:,5)=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

load('TestingonEntryExit_CareStart (Anterior Nasal Swab - FDA)_NoDelay_Hellewell.mat')
R(:,6)=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);


load('RAgTest_PlotOrder.mat');
% testName=OrderPlotTest;
NumTest=length(testName);
RAllAgTest=zeros(14,NumTest);

for kk=1:NumTest
    load(['TestingonEntryExit_' testName{kk} '_NoDelay_Hellewell.mat'],'IDSLS','IDSLA');
    RAllAgTest(:,kk)=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);
end

AgBetter=zeros(14,1);

for qq=1:14
    fAgBetter=find(RAllAgTest(qq,:)<=R(qq,1));
    AgBetter(qq)=length(fAgBetter);
end

subplot('Position',[0.59,0.555,0.79885593220339./2,0.39]);

MF={'-p','-.s','-.d','-.^','-.o','-.h'};
for ii=1:6
   plot(q,sqrt(R(:,ii)),MF{ii},'LineWidth',2,'MarkerFacecolor',CTest(ii,:),'color',CTest(ii,:),'Markersize',10); hold on;
end
box off;
grid on;
% legend({'RT-PCR','LumiraDx','Sofia','BinaxNOW','BD Veritor'},'Fontsize',18,'NumColumns',3,'Location','NorthEast');
% legend boxoff;
yytick=[0 0.001 0.01 0.025 0.05 0.075 0.1:0.05:0.3];
set(gca,'LineWidth',2,'Tickdir','out','Fontsize',20,'XTick',[1:14],'Xlim',[1 14],'Xminortick','off','YTick',sqrt(yytick),'ylim',[0 sqrt(yytick(end))],'YTickLabel',num2str(yytick'),'Yminortick','on');
% xlabel('Duration of quarantine (days)','Fontsize',26);
ylabel({'Probability of','post-quarantine transmission'},'Fontsize',26);
title('Single rapid antigen test on both entry and exit','Fontsize',26);
text(-1.89,0.585,'B','Fontsize',34,'FontWeight','bold');

legend({'RT-PCR','LumiraDx','Sofia','BinaxNOW','BD Veritor','CareStart'},'Fontsize',18,'NumColumns',3,'Location','NorthEast');
legend boxoff;


subplot('Position',[0.59,0.105,0.79885593220339./2,0.39]);
plot(q,AgBetter,'k-o','LineWidth',2,'MarkerFacecolor','k','color','k','Markersize',10); hold on;

box off;
grid on;
%%%legend({'LumiraDx','Sofia','BinaxNOW','BD Veritor'},'Fontsize',18,'NumColumns',3,'Location','NorthEast');
% % % legend boxoff;

YTickLabelText=cell(NumTest+1,1);
for kk=0:NumTest
    if(rem(kk,2)==0)
        YTickLabelText{kk+1}=[ num2str(kk) '/' num2str(NumTest)];
    else
        YTickLabelText{kk+1}=' ';
    end
end
set(gca,'LineWidth',2,'Tickdir','out','Fontsize',20,'XTick',[1:14],'Xlim',[1 14],'YTick',[0:NumTest],'Yticklabel',{YTickLabelText{:}},'Yminortick','off','Xminortick','off');
xlabel('Duration of quarantine (days)','Fontsize',26);
ylim([0 NumTest]);
ylabel({'Fraction of antigen tests','outperforming RT-PCR'},'Fontsize',24);
text(-1.89,18.275,'D','Fontsize',34,'FontWeight','bold');

print(gcf,'Figure2','-dpng','-r300');
