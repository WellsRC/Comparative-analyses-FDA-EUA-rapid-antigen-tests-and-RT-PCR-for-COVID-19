%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Plots the R_Eff of the various testing frequencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
clear;
clc;
close all;

addpath([pwd '\Delta_Variant']);
addpath([pwd '\Delta_Variant\Results']);

[pA,~,R0,ts,td] = BaselineParameters;
CTest=[hex2rgb('#231B12');hex2rgb('#375E97');hex2rgb('#486824');hex2rgb('#F9A603');hex2rgb('#CF3721');hex2rgb('#810f7c')];
Freq=[1:14];

RNoTest=integral(@(t)InfectiousnessfromInfection(t,R0,R0,pA,ts,td,1),0,inf);

R=zeros(length(Freq),6);
R_UN=zeros(length(Freq),1000,6);
% RT-PCR with one day delay
load([num2str(1) '-day_Delay_Testing_Frequency_RTPCR_DeltaVOC.mat']);
R(:,1)=(1-pA).*RTotS+pA.*RTotA;
RRTPCRDelay=R(:,1);

load([num2str(1) '-day_Delay_Testing_Frequency_RTPCR_DeltaVOC_Uncertainty.mat']);
R_UN(:,:,1)=(1-pA).*RTotSv+pA.*RTotAv;


load('Testing_Frequency_LumiraDX (Anterior Nasal Swab)_DeltaVOC.mat')
R(:,2)=(1-pA).*RTotS+pA.*RTotA;
load('Testing_Frequency_LumiraDX (Anterior Nasal Swab)_DeltaVOC_Uncertainty.mat')
R_UN(:,:,2)=(1-pA).*RTotSv+pA.*RTotAv;

load('Testing_Frequency_Sofia (FDA)_DeltaVOC.mat')
R(:,3)=(1-pA).*RTotS+pA.*RTotA;
load('Testing_Frequency_Sofia (FDA)_DeltaVOC_Uncertainty.mat')
R_UN(:,:,3)=(1-pA).*RTotSv+pA.*RTotAv;


load('Testing_Frequency_BinaxNOW (FDA)_DeltaVOC.mat')
R(:,4)=(1-pA).*RTotS+pA.*RTotA;
load('Testing_Frequency_BinaxNOW (FDA)_DeltaVOC_Uncertainty.mat')
R_UN(:,:,4)=(1-pA).*RTotSv+pA.*RTotAv;

load('Testing_Frequency_BD Veritor_DeltaVOC.mat')
R(:,5)=(1-pA).*RTotS+pA.*RTotA;
load('Testing_Frequency_BD Veritor_DeltaVOC_Uncertainty.mat')
R_UN(:,:,5)=(1-pA).*RTotSv+pA.*RTotAv;

load('Testing_Frequency_CareStart (Anterior Nasal Swab - FDA)_DeltaVOC.mat')
R(:,6)=(1-pA).*RTotS+pA.*RTotA;
load('Testing_Frequency_CareStart (Anterior Nasal Swab - FDA)_DeltaVOC_Uncertainty.mat')
R_UN(:,:,6)=(1-pA).*RTotSv+pA.*RTotAv;


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

LB_AGTest=zeros(size(R));
UB_AGTest=zeros(size(R));

for jj=1:6
    for ii=1:length(LB_AGTest)
        [~,LB_AGTest(ii,jj),UB_AGTest(ii,jj)]=Credible_Interval_High_Density(R(ii,jj),R_UN(ii,:,jj),0.95,'continuous',[0 RNoTest]);      
    end
end

for ii=1:6
   b(ii).FaceColor=CTest(ii,:); 
   errorbar(b(ii).XEndPoints,b(ii).YEndPoints,b(ii).YEndPoints'-LB_AGTest(:,ii),UB_AGTest(:,ii)-b(ii).YEndPoints','.','Markersize',10^(-16),'LineWidth',2,'color',[0.75 0.75 0.75]);
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
    load([num2str(zz) '-day_Delay_Testing_Frequency_RTPCR_DeltaVOC.mat']);
    R(:,zz+1)=(1-pA).*RTotS+pA.*RTotA;
    
    
    load([num2str(zz) '-day_Delay_Testing_Frequency_RTPCR_DeltaVOC_Uncertainty.mat']);
    RU(:,:,zz+1)=(1-pA).*RTotSv+pA.*RTotAv;
end

LB_RTPCRTest=zeros(size(R));
UB_RTPCRTest=zeros(size(R));

for jj=1:6
    for ii=1:length(LB_AGTest)  
        [~,LB_RTPCRTest(ii,jj),UB_RTPCRTest(ii,jj)]=Credible_Interval_High_Density(R(ii,jj),RU(ii,:,jj),0.95,'continuous',[0 RNoTest]);    
    end
end
subplot('Position',[0.572058823529412,0.586605876393113,0.37,0.4]);


b=bar(Freq,R,'LineStyle','none');
hold on;
plot(linspace(0.5,14.5,1001),ones(1001,1),'-.','color',[0.75 0.75 0.75],'LineWidth',2)

plot(linspace(0.5,14.5,1001),RNoTest.*ones(1001,1),'-','color',[0.75 0.75 0.75],'LineWidth',2)
box off;

text(18.07,1.282/1.3*3.2,'B','Fontsize',34,'FontWeight','bold');
xlim([0.5 14.5]);
fprintf('====================================== \n');
fprintf('RT-PCR \n');
fprintf('====================================== \n');
for ii=1:6
   b(ii).FaceColor=CTest(ii,:); 
   errorbar(b(ii).XEndPoints,b(ii).YEndPoints,b(ii).YEndPoints'-LB_RTPCRTest(:,ii),UB_RTPCRTest(:,ii)-b(ii).YEndPoints','.','Markersize',10^(-16),'LineWidth',2,'color',[0.75 0.75 0.75]);
   if(ii==3)
       testFreq=[b(ii).YEndPoints' LB_RTPCRTest(:,ii) UB_RTPCRTest(:,ii)];
       fprintf('Bounds for a delay of 2 days: %3.2f (%3.2f - %3.2f) \n',testFreq(1,:))
   elseif(ii==2)
       testFreq=[b(ii).YEndPoints' LB_RTPCRTest(:,ii) UB_RTPCRTest(:,ii)];
       fprintf('Bounds for a delay of 24 hrs: %3.2f (%3.2f - %3.2f) \n',testFreq(1,:))       
   elseif(ii==1)
       D=zeros(1000,1);
       for jj=1:1000
           D(jj)=find(RU(:,jj,ii)<1,1,'last');
       end
       MLED=find(R(:,ii)<1,1,'last');
       [~,LB_D,UB_D]=Credible_Interval_High_Density(MLED,D(:),0.95,'discrete',[1 14]);    
       fprintf('Bounds for time between tests no-delay: %d (%d - %d) \n',[MLED LB_D UB_D]);
       
       testFreq=[b(ii).YEndPoints' LB_RTPCRTest(:,ii) UB_RTPCRTest(:,ii)];
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

RAllAgTestU=zeros(14,1000,NumTest);


for kk=1:NumTest
    load(['Testing_Frequency_' testName{kk} '_DeltaVOC.mat'],'RTotS','RTotA');
    RAllAgTest(:,kk)=(1-pA).*RTotS+pA.*RTotA;
    
    load(['Testing_Frequency_' testName{kk} '_DeltaVOC_Uncertainty.mat'],'RTotSv','RTotAv')
    RAllAgTestU(:,:,kk)=(1-pA).*RTotSv+pA.*RTotAv;
end


AgMRedR=zeros(14,1);
PCRMRedR=zeros(14,1);
AgRLT1=zeros(14,1);

AgMRedR_Un=zeros(1000,14);
PCRMRedR_Un=zeros(1000,14);
AgRLT1U=zeros(1000,14);

for qq=1:14
    AgMRedR(qq)=mean(1-(RAllAgTest(qq,:)./R0));
    PCRMRedR(qq)=(1-(RRTPCRDelay(qq)./R0));
    fAgLT1=find(RAllAgTest(qq,:)<1);
    AgRLT1(qq)=length(fAgLT1);
    
    for zz=1:1000
        fAgBetter=find(RAllAgTestU(qq,zz,:)<=RU(qq,zz,2));
        AgMRedR_Un(zz,qq)=mean(1-(RAllAgTestU(qq,zz,:)./R0));
        PCRMRedR_Un(zz,qq)=(1-(RU(qq,zz,2)./R0));
        fAgLT1=find(RAllAgTestU(qq,zz,:)<1);
        AgRLT1U(zz,qq)=length(fAgLT1);
    end
end

fprintf('====================================== \n');
fprintf('Maintain R_E <1  \n');
fprintf('====================================== \n');

MF=zeros(1000,1);
for ii=1:1000
    MF(ii)=max(Freq(AgRLT1U(ii,:)==18));
end
[~,LBPQT,UBPQT]=Credible_Interval_High_Density(max(Freq(AgRLT1==18)),MF,0.95,'discrete',[1 14]); 
testFreq=[max(Freq(AgRLT1==18)) LBPQT UBPQT];
fprintf('Bounds for maximum time between testing for all RA tests such that RE<1: %2.0f (%2.0f - %2.0f) \n',testFreq) 

MF=zeros(1000,1);
for ii=1:1000
    MF(ii)=min(Freq(AgRLT1U(ii,:)==0));
end


[~,LBPQT,UBPQT]=Credible_Interval_High_Density(min(Freq(AgRLT1==0)),MF,0.95,'discrete',[1 14]); 
testFreq=[min(Freq(AgRLT1==0)) LBPQT UBPQT];
fprintf('Bounds for minimum time between testing for all RA tests such that RE>1: %2.0f (%2.0f - %2.0f) \n',testFreq) 

subplot('Position',[0.572058823529412,0.093566362715299,0.37,0.4]);

LB_Ag=zeros(size(AgRLT1));
UB_Ag=zeros(size(AgRLT1));

for ii=1:length(LB_Ag)
    [~,LB_Ag(ii),UB_Ag(ii)]=Credible_Interval_High_Density(AgRLT1(ii),AgRLT1U(:,ii),0.95,'discrete',[0 18]);    
end

errorbar(Freq,AgRLT1,AgRLT1-LB_Ag,UB_Ag-AgRLT1,'-o','color',hex2rgb('#128277'),'MarkerSize',10,'MarkerEdgeColor',hex2rgb('#128277'),'MarkerFaceColor',hex2rgb('#128277'),'LineWidth',2);hold on;
   
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


LB_RTPCRMedR=zeros(size(PCRMRedR));
UB_RTPCRMedR=zeros(size(PCRMRedR));


LB_AgMedR=zeros(size(AgMRedR));
UB_AgMedR=zeros(size(AgMRedR));

for ii=1:length(LB_Ag)
    [~,LB_AgMedR(ii),UB_AgMedR(ii)]=Credible_Interval_High_Density(AgMRedR(ii),AgMRedR_Un(:,ii),0.95,'continuous',[0 1]);    
    [~,LB_RTPCRMedR(ii),UB_RTPCRMedR(ii)]=Credible_Interval_High_Density(PCRMRedR(ii),PCRMRedR_Un(:,ii),0.95,'continuous',[0 1]);  
end

p1=errorbar(Freq,100.*PCRMRedR,100.*PCRMRedR-100.*LB_RTPCRMedR,100.*UB_RTPCRMedR-100.*PCRMRedR,'k-p','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',2);hold on;   
p2=errorbar(Freq,100.*AgMRedR,100.*AgMRedR-100.*LB_AgMedR,100.*UB_AgMedR-100.*AgMRedR,'-o','color',hex2rgb('#128277'),'MarkerSize',10,'MarkerEdgeColor',hex2rgb('#128277'),'MarkerFaceColor',hex2rgb('#128277'),'LineWidth',2);hold on;
plot(linspace(0,6.5,101),100.*(1-RNoTest./R0).*ones(101,1),'-','LineWidth',2,'color',[0.75 0.75 0.75]);
plot(linspace(8.5,14,101),100.*(1-RNoTest./R0).*ones(101,1),'-','LineWidth',2,'color',[0.75 0.75 0.75]);
text(7.5,100.*(1-RNoTest./R0),'No test','horizontalalignment','center','verticalalignment','middle','color',[0.75 0.75 0.75],'Fontsize',20);
box off;
grid on;
legend([p1 p2],{'RT-PCR (1-day delay)','Rapid antigen test (Average)'},'Fontsize',18,'NumColumns',3,'Location','NorthWest');
legend boxoff;

YTickLabelText=cell(NumTest+1,1);
for kk=0:NumTest
    if(rem(kk,2)==0)
        YTickLabelText{kk+1}=[ num2str(kk) '/' num2str(NumTest)];
    else
        YTickLabelText{kk+1}=' ';
    end
end

set(gca,'LineWidth',2,'Tickdir','out','Fontsize',20,'XTick',[1:14],'Xlim',[1 14],'ylim',[0 100],'Ytick',[0:10:100],'Yminortick','on','Xminortick','off','xdir', 'reverse' );
xlabel('Duration of quarantine (days)','Fontsize',26);
ylabel({'Reduction of effective','reproduction number'},'Fontsize',26);
text(17.35,18.315*100/18,0,'C','Fontsize',34,'FontWeight','bold');
ytickformat('percentage')
xlabel('Days between tests','Fontsize',24);


fprintf('====================================== \n');
fprintf('Reduction in R_E \n');
fprintf('====================================== \n');
testFreq=[100.*AgMRedR 100.*LB_AgMedR 100.*UB_AgMedR];
fprintf('Bounds for testing every two days with RA test: %3.1f (%3.1f - %3.1f) \n',testFreq(2,:))  

testFreq=[100.*PCRMRedR 100.*LB_RTPCRMedR 100.*UB_RTPCRMedR];
fprintf('Bounds for testing every two days with RT-PCR test: %3.1f (%3.1f - %3.1f) \n',testFreq(2,:))  
rmpath([pwd '\Delta_Variant']);
rmpath([pwd '\Delta_Variant\Results']);

print(gcf,'Figure3','-dpng','-r300');
print(gcf,'Figure3.tiff','-dtiff','-r600');
