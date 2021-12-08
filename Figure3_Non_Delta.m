%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Plots the R_Eff of the various testing frequencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
close all;

R0Sim=3.2;
addpath([pwd '\Non_Delta']);
addpath([pwd '\Non_Delta\Results']);

[pA,~,R0,ts,td] = BaselineParameters;
CTest=[hex2rgb('#231B12');hex2rgb('#375E97');hex2rgb('#486824');hex2rgb('#F9A603');hex2rgb('#CF3721');hex2rgb('#810f7c')];
Freq=[1:14];

RNoTest=integral(@(t)InfectiousnessfromInfection(t,R0Sim,R0Sim,pA,ts,td,1),0,inf);

R=zeros(length(Freq),6);
R_UN=zeros(length(Freq),1000,6);
% RT-PCR with one day delay
load([num2str(1) '-day_Delay_Testing_Frequency_RTPCR_General.mat']);
R(:,1)=(1-pA).*RTotS+pA.*RTotA;
RRTPCRDelay=R(:,1);

% load([num2str(1) '-day_Delay_Testing_Frequency_RTPCR_General_Uncertainty.mat']);
% R_UN(:,:,1)=(1-pA).*RTotSv+pA.*RTotAv;


load('Testing_Frequency_LumiraDX (Anterior Nasal Swab)_General.mat')
R(:,2)=(1-pA).*RTotS+pA.*RTotA;
% load('Testing_Frequency_LumiraDX (Anterior Nasal Swab)_General_Uncertainty.mat')
% R_UN(:,:,2)=(1-pA).*RTotSv+pA.*RTotAv;

load('Testing_Frequency_Sofia (FDA)_General.mat')
R(:,3)=(1-pA).*RTotS+pA.*RTotA;
% load('Testing_Frequency_Sofia (FDA)_General_Uncertainty.mat')
% R_UN(:,:,3)=(1-pA).*RTotSv+pA.*RTotAv;


load('Testing_Frequency_BinaxNOW (FDA)_General.mat')
R(:,4)=(1-pA).*RTotS+pA.*RTotA;
% load('Testing_Frequency_BinaxNOW (FDA)_General_Uncertainty.mat')
% R_UN(:,:,4)=(1-pA).*RTotSv+pA.*RTotAv;

load('Testing_Frequency_BD Veritor_General.mat')
R(:,5)=(1-pA).*RTotS+pA.*RTotA;
% load('Testing_Frequency_BD Veritor_General_Uncertainty.mat')
% R_UN(:,:,5)=(1-pA).*RTotSv+pA.*RTotAv;

load('Testing_Frequency_CareStart (Anterior Nasal Swab - FDA)_General.mat')
R(:,6)=(1-pA).*RTotS+pA.*RTotA;
% load('Testing_Frequency_CareStart (Anterior Nasal Swab - FDA)_General_Uncertainty.mat')
% R_UN(:,:,6)=(1-pA).*RTotSv+pA.*RTotAv;

R=(R0Sim./R0).*R;
RRTPCRDelay=(R0Sim./R0).*RRTPCRDelay;
figure('units','normalized','outerposition',[0.05 0.05 1 1]);
subplot('Position',[0.097689075630252,0.586605876393113,0.37,0.4]);

b=bar(Freq,R,'LineStyle','none');
hold on;
plot(linspace(0.5,14.5,1001),ones(1001,1),'-.','color',[0.75 0.75 0.75],'LineWidth',2)
plot(linspace(0.5,14.5,1001),RNoTest.*ones(1001,1),'-','color',[0.75 0.75 0.75],'LineWidth',2)
box off;
% xlabel('Days between tests','Fontsize',24);
text(18.07,1.282./1.3*3.2,'A','Fontsize',34,'FontWeight','bold');
xlim([0.5 14.5]);
for ii=1:6
   b(ii).FaceColor=CTest(ii,:); 
%    errorbar(b(ii).XEndPoints,b(ii).YEndPoints,b(ii).YEndPoints-prctile(R_UN(:,:,ii),2.5),prctile(R_UN(:,:,ii),97.5)-b(ii).YEndPoints,'.','Markersize',10^(-16),'color',[0.75 0.75 0.75]);
end
legend({'RT-PCR (1-day delay)','LumiraDx','Sofia','BinaxNOW','BD Veritor','CareStart'},'Fontsize',18,'NumColumns',3,'Position',[0.100665266106443,0.921816961567413,0.326155454859513,0.067375884652742]);
legend boxoff;
set(gca,'LineWidth',2,'Tickdir','out','Fontsize',20,'XTick',[1:14],'Xminortick','off','YTick',[0:0.25:2.5],'ylim',[0 3.2],'Yminortick','on','xdir', 'reverse');

ylabel({'Effective','reproduction number'},'Fontsize',24);

CTest=[hex2rgb('#231B12');hex2rgb('525252');hex2rgb('737373');hex2rgb('969696');hex2rgb('bdbdbd');hex2rgb('d9d9d9')];
Freq=[1:14];
R=zeros(length(Freq),6);
RU=zeros(length(Freq),1000,6);

for zz=0:5
    load([num2str(zz) '-day_Delay_Testing_Frequency_RTPCR_General.mat']);
    R(:,zz+1)=(1-pA).*RTotS+pA.*RTotA;
    
    
    load([num2str(zz) '-day_Delay_Testing_Frequency_RTPCR_General_Uncertainty.mat']);
    RU(:,:,zz+1)=(1-pA).*RTotSv+pA.*RTotAv;
end

R=(R0Sim./R0).*R;

RU=(R0Sim./R0).*RU;

subplot('Position',[0.572058823529412,0.586605876393113,0.37,0.4]);


b=bar(Freq,R,'LineStyle','none');
hold on;
plot(linspace(0.5,14.5,1001),ones(1001,1),'-.','color',[0.75 0.75 0.75],'LineWidth',2)

plot(linspace(0.5,14.5,1001),RNoTest.*ones(1001,1),'-','color',[0.75 0.75 0.75],'LineWidth',2)
box off;

text(18.07,1.282/1.3*3.2,'B','Fontsize',34,'FontWeight','bold');
xlim([0.5 14.5]);
for ii=1:6
   b(ii).FaceColor=CTest(ii,:); 
   errorbar(b(ii).XEndPoints,b(ii).YEndPoints,b(ii).YEndPoints'-prctile(RU(:,:,ii),2.5,2),prctile(RU(:,:,ii),97.5,2)-b(ii).YEndPoints','.','Markersize',10^(-16),'LineWidth',2,'color',[0.75 0.75 0.75]);
   if(ii==3)
       testFreq=[b(ii).YEndPoints' prctile(RU(:,:,ii),2.5,2) prctile(RU(:,:,ii),97.5,2)];
       fprintf('Bounds for a delay of 2 days: %3.2f (%3.2f - %3.2f) \n',testFreq(1,:))
   elseif(ii==2)
       testFreq=[b(ii).YEndPoints' prctile(RU(:,:,ii),2.5,2) prctile(RU(:,:,ii),97.5,2)];
       fprintf('Bounds for a delay of 24 hrs: %3.2f (%3.2f - %3.2f) \n',testFreq(1,:))       
   elseif(ii==1)
       D=zeros(1000,1);
       for jj=1:1000
           D(jj)=find(RU(:,jj,ii)<1,1,'last');
       end
       MLED=find(R(:,ii)<1,1,'last');
       fprintf('Bounds for time between tests no-delay: %d (%d - %d) \n',[MLED prctile(D,[2.5 97.5])]);
       
       testFreq=[b(ii).YEndPoints' prctile(RU(:,:,ii),2.5,2) prctile(RU(:,:,ii),97.5,2)];
       fprintf('Bounds for NO delay and testing every three days: %3.2f (%3.2f - %3.2f) \n',testFreq(3,:))  
   end
end
legend({'No delay','1-day delay','2-day delay','3-day delay','4-day delay','5-day delay'},'Fontsize',18,'NumColumns',3,'Location','NorthWest');
legend boxoff;
set(gca,'LineWidth',2,'Tickdir','out','Fontsize',20,'XTick',[1:14],'Xminortick','off','YTick',[0:0.25:2.5],'ylim',[0 3.2],'Yminortick','on','xdir', 'reverse' );

ylabel({'Effective','reproduction number'},'Fontsize',24);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Comapre test to RT-PCR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5



%%%%%%%%%%%%%%%%%%%%%%%
% All tests
%%%%%%%%%%%%%%%%%%%%%%
load('RAgTest_PlotOrder.mat');
% testName=OrderPlotTest;
NumTest=length(testName);
RAllAgTest=zeros(14,NumTest);


for kk=1:NumTest
    load(['Testing_Frequency_' testName{kk} '_General.mat'],'RTotS','RTotA');
    RAllAgTest(:,kk)=(1-pA).*RTotS+pA.*RTotA;
end

RAllAgTest=(R0Sim./R0).*RAllAgTest;

AgBetter=zeros(14,1);
AgRLT1=zeros(14,1);

for qq=1:14
    fAgBetter=find(RAllAgTest(qq,:)<=RRTPCRDelay(qq,1));
    AgBetter(qq)=length(fAgBetter);
    fAgLT1=find(RAllAgTest(qq,:)<1);
    AgRLT1(qq)=length(fAgLT1);
end

subplot('Position',[0.572058823529412,0.093566362715299,0.37,0.4]);
   plot(Freq,AgRLT1,'k-o','LineWidth',2,'MarkerFacecolor','k','color','k','Markersize',10); hold on;

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

set(gca,'LineWidth',2,'Tickdir','out','Fontsize',20,'XTick',[1:14],'Xlim',[1 14],'YTick',[0:NumTest],'Yticklabel',{YTickLabelText{:}},'Yminortick','off','Xminortick','off','xdir', 'reverse' );
xlabel('Duration of quarantine (days)','Fontsize',26);
ylim([0 NumTest]);
ylabel({'Fraction of antigen','tests with {\it R_E} < 1'},'Fontsize',26);
text(17.35,18.315,'D','Fontsize',34,'FontWeight','bold');

xlabel('Days between tests','Fontsize',24);

subplot('Position',[0.097689075630252,0.093566362715299,0.37,0.4]);
   plot(Freq,AgBetter,'k-o','LineWidth',2,'MarkerFacecolor','k','color','k','Markersize',10); hold on;

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

set(gca,'LineWidth',2,'Tickdir','out','Fontsize',20,'XTick',[1:14],'Xlim',[1 14],'YTick',[0:NumTest],'Yticklabel',{YTickLabelText{:}},'Yminortick','off','Xminortick','off','xdir', 'reverse' );
xlabel('Duration of quarantine (days)','Fontsize',26);
ylim([0 NumTest]);
ylabel({'Fraction of antigen tests','outperforming RT-PCR'},'Fontsize',26);
text(17.35,18.315,0,'C','Fontsize',34,'FontWeight','bold');

xlabel('Days between tests','Fontsize',24);

rmpath([pwd '\Non_Delta']);
rmpath([pwd '\Non_Delta\Results']);

print(gcf,'Figure3_Non_Delta','-dpng','-r300');
