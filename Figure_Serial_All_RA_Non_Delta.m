%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Plots the R_Eff of the various testing frequencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
clear;
clc;
close all;

addpath([pwd '\Non_Delta']);
addpath([pwd '\Non_Delta\Results']);

[pA,~,R0,ts,td] = BaselineParameters;


RNoTest=integral(@(t)InfectiousnessfromInfection(t,R0,R0,pA,ts,td,1),0,inf);

Freq=[1:14];
load('RAgTest_PlotOrder.mat');
NumTest=length(testName);

R=zeros(length(Freq),NumTest+1);
R_UN=zeros(length(Freq),1000,NumTest+1);

load([num2str(1) '-day_Delay_Testing_Frequency_RTPCR_General.mat']);
R(:,1)=(1-pA).*RTotS+pA.*RTotA;

load([num2str(1) '-day_Delay_Testing_Frequency_RTPCR_General_Uncertainty.mat']);
R_UN(:,:,1)=(1-pA).*RTotSv+pA.*RTotAv;

for ii=1:NumTest
    load(['Testing_Frequency_' testName{ii} '_General.mat'],'RTotS','RTotA');
    R(:,ii+1)=(1-pA).*RTotS+pA.*RTotA;
    
    load(['Testing_Frequency_' testName{ii} '_General_Uncertainty.mat'],'RTotSv','RTotAv')
    R_UN(:,:,ii+1)=(1-pA).*RTotSv+pA.*RTotAv;
end

for fr=4:-1:1
    if(rem(fr,2)==0)
        figure('units','normalized','outerposition',[0 0 1 1]);
        subplot('Position',[0.0764,0.5,0.823220338983051/2,0.441236068895644]);
    else
        subplot('Position',[0.584,0.5,0.823220338983051/2,0.441236068895644]);
    end

    b=bar([1:(NumTest+1)],R(fr,:),'LineStyle','none');
    hold on;
    plot(linspace(0.5,(NumTest+1)+0.5,1001),ones(1001,1),'-.','color',[0.75 0.75 0.75],'LineWidth',2)
    box off;
    xlim([0.5 (NumTest+1)+0.5]);

    b.FaceColor = 'flat';

    testname2=cell((NumTest+1),1);
    
    for ii=1:(NumTest+1)
        if(ii==1)
            [C,~]= ColourTests('RTPCR');            
            testname2(1)={'RT-PCR (24-h delay)'};
        else
            [C,~]= ColourTests(testName{ii-1});
            if(length(C(:,1))>1)
                C=C(end,:);
            end
            testname2(ii)={AdjustedNames_Plotting(testName{ii-1})};
        end
       b.CData(ii,:)=C; 
    end
    
    LB_AGTest=zeros(size(R(fr,:)));
    UB_AGTest=zeros(size(R(fr,:)));

    for jj=1:(NumTest+1)
            [~,LB_AGTest(jj),UB_AGTest(jj)]=Credible_Interval_High_Density(R(fr,jj),R_UN(fr,:,jj),0.95,'continuous',[0 RNoTest]); 
    end
    errorbar(b.XEndPoints,b.YEndPoints,b.YEndPoints-LB_AGTest,UB_AGTest-b.YEndPoints,'.','Markersize',10^(-16),'LineWidth',2,'color',[0.75 0.75 0.75]);

    if(fr<=2)
        set(gca,'LineWidth',2,'Tickdir','out','Fontsize',20,'XTick',[1:(NumTest+1)],'Xminortick','off','YTick',[0:0.3:1.9],'ylim',[0 1.9],'Yminortick','on','XTickLabel',testname2);
        xtickangle(90)
        xlabel('Test name','Fontsize',26,'Position',[9 -1.45 0]);
    else
        set(gca,'LineWidth',2,'Tickdir','out','Fontsize',20,'XTick',[1:(NumTest+1)],'Xminortick','off','YTick',[0:0.3:1.9],'ylim',[0 1.9],'Yminortick','on','XTickLabel','');
    end
    ylabel({'Effective','reproduction number'},'Fontsize',26,'Position',[-0.998673431830746,0.95,-1]);
    if(fr==1)
        title('Testing every day','Fontsize',28);
    else
        title(['Testing every ' num2str(fr) ' days'],'Fontsize',28);
    end
    tfr=1+(4-fr);
    text(-2.699383553973718,1.489425287356321./1.4.*1.9,char(64+tfr),'Fontsize',34,'FontWeight','bold');
    if(rem(fr,2)==1)
        print(gcf,['Alltests_Testing_Frequency=' num2str(fr) '_Non_Delta.png'],'-dpng','-r600');
    end
end


rmpath([pwd '\Non_Delta']);
rmpath([pwd '\Non_Delta\Results']);