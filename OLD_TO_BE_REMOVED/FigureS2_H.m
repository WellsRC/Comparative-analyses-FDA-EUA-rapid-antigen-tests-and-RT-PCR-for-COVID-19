%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Plots the R_Eff of the various testing frequencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
close all;
Freq=[1:14];
load('RAgTest_PlotOrder.mat');
NumTest=length(testName);

R=zeros(length(Freq),NumTest+1);
load([num2str(1) '-day_Delay_Testing_Frequency_RTPCR_Hellewell.mat']);
R(:,1)=(1-pA).*RTotS+pA.*RTotA;


for ii=1:NumTest
    load(['Testing_Frequency_' testName{ii} '_Hellewell.mat'],'RTotS','RTotA');
    R(:,ii+1)=(1-pA).*RTotS+pA.*RTotA;
end

for fr=9:-1:6
    if(rem(fr,2)==1)
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
    testname2(1)={'RT-PCR (one-day delay)'};
    testname2(2:end)=testName;
    for ii=1:(NumTest+1)
        if(ii==1)
            [C,~]= ColourTests('RTPCR');
        else
            [C,~]= ColourTests(testName{ii-1});
            if(length(C(:,1))>1)
                C=C(end,:);
            end
        end
       b.CData(ii,:)=C; 
    end

    if(fr<=7)
        set(gca,'LineWidth',2,'Tickdir','out','Fontsize',20,'XTick',[1:(NumTest+1)],'Xminortick','off','YTick',[0:0.2:1.4],'ylim',[0 1.4],'Yminortick','on','XTickLabel',testname2);
        xtickangle(90)
        xlabel('Test name','Fontsize',26,'Position',[9 -1.45 0]);
    else
        set(gca,'LineWidth',2,'Tickdir','out','Fontsize',20,'XTick',[1:(NumTest+1)],'Xminortick','off','YTick',[0:0.2:1.4],'ylim',[0 1.4],'Yminortick','on','XTickLabel','');
    end
    ylabel({'Effective','reproduction number'},'Fontsize',26);
    if(fr==1)
        title('Testing every day','Fontsize',28);
    else
        title(['Testing every ' num2str(fr) ' days'],'Fontsize',28);
    end
    tfr=6+(9-fr);
    text(-2.699383553973718,1.489425287356321,char(59+tfr),'Fontsize',34,'FontWeight','bold');
    if(rem(fr,2)==0)
        print(gcf,['Alltests_Testing_Frequency=' num2str(fr) '.png'],'-dpng','-r600');
    end
end
