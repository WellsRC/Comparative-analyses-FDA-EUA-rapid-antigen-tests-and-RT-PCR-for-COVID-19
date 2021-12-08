clear;
close all;

addpath([pwd '\Delta_Variant']);
addpath([pwd '\Delta_Variant\Results']);
[pAb,~,R0b,~,~] = BaselineParameters;

load('RAgTest_Name.mat','testName');
NumTests=length(testName);
for TestN=1:NumTests
    testNameP=testName{TestN};
    load('ColourMap_Quarantine_Scenario_pA_R0.mat','cc')
    load('ColourMap_Serial_Testing_Scenario_pA_R0.mat','ss')
    hh=0.75;
    figure('units','normalized','outerposition',[0 0.1 1 hh]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Quarantine
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax1=subplot('Position',[0.072226890756302,0.15,0.35,0.576469098277617/hh]);
    [pA,R0v,QD] = Duration_Quarantine_Test(testNameP);
    contourf(R0v,pA,QD,[1:14],'LineStyle','none'); hold on;

    scatter(R0b,pAb,100,'ko','filled')
    caxis([0.5 14.5])
    colormap(ax1,cc);
    box off;
    set(gca,'LineWidth',2,'tickdir','out','YTick',[0.1:0.05:0.95],'XTick',[1:0.5:5],'Xminortick','on','Fontsize',24);
    ylabel('Proportion asymptomatic','Fontsize',26);
    xlabel('Basic reproduction number','Fontsize',26);
    title(['Quarantine: ' testNameP ],'Fontsize',26);
    c=colorbar;
    c.Position=[0.437083332981152,0.15,0.011204481792717,0.576469098277617/hh];
    c.Ticks=[1:14];
    c.Label.String='Duration of quarantine (days)';
    c.Label.FontSize=26;
    c.Label.Rotation=270;
    c.Label.Position=[4.972698438735235,7.500006675720215,0];
    text(0.475375375375376,1.01016333938294,0,'A','Fontsize',32,'FontWeight','bold');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Serial testing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax2=subplot('Position',[0.577731092436973,0.15,0.35,0.576469098277617/hh]);
    [pA,R0v,FS] = Frequency_Serial_Test(testNameP);
    contourf(R0v,pA,FS,[0:14],'LineStyle','none');
    hold on;
    scatter(R0b,pAb,100,'yo','filled')
    caxis([-0.5 14.5])
    colormap(ax2,ss);
    box off;
    set(gca,'LineWidth',2,'tickdir','out','YTick',[0.1:0.05:0.95],'XTick',[1:0.5:5],'Xminortick','on','Fontsize',24);
    ylabel('Proportion asymptomatic','Fontsize',26);
    xlabel('Basic reproduction number','Fontsize',26);
    title(['Serial testing: ' testNameP ],'Fontsize',26);
    c=colorbar;
        c.Position=[0.942335433821487,0.15,0.011204481792717,0.576368866986917/hh];
        c.Ticks=[0:14];
        c.Label.String='Days between tests';
    c.Label.FontSize=26;
    c.Label.Rotation=270;
    c.Label.Position=[4.972698438735235,7.500006675720215,0];
    c.TickLabels{1}='R > 1';
    text(0.475375375375376,1.01016333938294,0,'B','Fontsize',32,'FontWeight','bold');
    print(gcf,['pA_R0_' testNameP '.png'],'-dpng','-r600');
end

rmpath([pwd '\Delta_Variant']);
rmpath([pwd '\Delta_Variant\Results']);