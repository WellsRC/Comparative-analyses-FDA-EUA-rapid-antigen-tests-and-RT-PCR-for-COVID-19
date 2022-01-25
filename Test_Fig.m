clear;
close all;
figure('units','normalized','outerposition',[0.1 0 0.65 1]);
subplot('Position',[0.074727272727273,0.725,0.381,0.265]);
set(gca,'LineWidth',2,'tickdir','out','Fontsize',20,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:0.1:1],'Ylim',[0 1]);
ylabel({'Diagnostic sensitivty'},'Fontsize',20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PPA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax1 = axes('Position',[0.227272727272727,0.81661600810537,0.227224980901452,0.172239108409349]);
    set(gca,'LineWidth',1.1,'tickdir','out','Fontsize',12,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:20:100],'Ylim',[0 100]);
    ytickformat('percentage')
    xlabel('Days since symptom onset','Fontsize',12);
    ylabel({'Percent','positive agreement'},'Fontsize',12);
    
subplot('Position',[0.6031,0.725,0.381,0.265]);
yytick=[0 0.005 0.025 0.05 0.1:0.05:0.3];
set(gca,'LineWidth',2,'Tickdir','out','Fontsize',20,'XTick',[1:14],'Xlim',[1 14],'Xminortick','off','YTick',sqrt(yytick),'ylim',[0 sqrt(yytick(end))],'YTickLabel',num2str(yytick'),'Yminortick','on');
ylabel({'Probability of','post-quarantine','transmission'},'Fontsize',20);

subplot('Position',[0.074727272727273,0.407,0.381,0.265]);
set(gca,'LineWidth',2,'tickdir','out','Fontsize',20,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:0.1:1],'Ylim',[0 1]);
ylabel({'Diagnostic sensitivty'},'Fontsize',20);
subplot('Position',[0.6031,0.407,0.381,0.265]);
yytick=[0 0.005 0.025 0.05 0.1:0.05:0.3];
set(gca,'LineWidth',2,'Tickdir','out','Fontsize',20,'XTick',[1:14],'Xlim',[1 14],'Xminortick','off','YTick',sqrt(yytick),'ylim',[0 sqrt(yytick(end))],'YTickLabel',num2str(yytick'),'Yminortick','on');
ylabel({'Probability of','post-quarantine','transmission'},'Fontsize',20);

subplot('Position',[0.074727272727273,0.087,0.381,0.265]);
set(gca,'LineWidth',2,'tickdir','out','Fontsize',20,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:0.1:1],'Ylim',[0 1]);
ylabel({'Diagnostic sensitivty'},'Fontsize',20);
xlabel({'Days post-infection'},'Fontsize',20);
subplot('Position',[0.6031,0.087,0.381,0.265]);
yytick=[0 0.005 0.025 0.05 0.1:0.05:0.3];
set(gca,'LineWidth',2,'Tickdir','out','Fontsize',20,'XTick',[1:14],'Xlim',[1 14],'Xminortick','off','YTick',sqrt(yytick),'ylim',[0 sqrt(yytick(end))],'YTickLabel',num2str(yytick'),'Yminortick','on');
ylabel({'Probability of','post-quarantine','transmission'},'Fontsize',20);
xlabel({'Duration of quarantine (days)'},'Fontsize',20);