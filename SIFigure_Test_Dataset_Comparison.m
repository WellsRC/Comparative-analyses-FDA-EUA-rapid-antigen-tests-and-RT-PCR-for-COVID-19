clear;
clc;
close all;
load('RAgTest_Name_DatasetCompare.mat');
testNamev=testName;
for ii=1:length(testNamev(:,1))
    
    figure('units','normalized','outerposition',[0.1 0 0.65 1]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    %% Baseline
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    testName=testNamev{ii,1};

    
    [NameP] = AdjustedNames_Plotting(testName);
    addpath([pwd '\Delta_Variant']);
    addpath([pwd '\Delta_Variant\Results']);
    [pA,~,~,ts,~] = BaselineParameters;
    [MLE_RTPCR,MLE_Ag,U_RTPCR,U_Ag,MLE_PPA,U_PPA,Dt,totalpos,truepos,w,t,CCtestRTPCR,SymPRTPCR,CCtest,SymP] = Sensitivity_for_Plotting(testNamev{ii,1},ts);
    [MLE_RTPCR_Ex,MLE_Ag_Ex,U_RTPCR_Ex,U_Ag_Ex,MLE_PPA_Ex,U_PPA_Ex,Dt_Ex,totalpos_Ex,truepos_Ex,w_Ex,t_Ex,CCtestRTPCR_Ex,SymPRTPCR_Ex,CCtest_Ex,SymP_Ex] = Sensitivity_for_Plotting(testNamev{ii,2},ts);
    CCtest=CCtest(1,:);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Diagnostic
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    subplot('Position',[0.074727272727273,0.725,0.381,0.265]);
    
    patch([t flip(t)],[prctile(U_Ag,2.5) flip(prctile(U_Ag,97.5))],CCtest,'LineStyle','none','Facealpha',0.2); hold on;
    plot(t,MLE_Ag,'-','color',CCtest,'LineWidth',2);
    
    patch([t_Ex flip(t_Ex)],[prctile(U_Ag_Ex,2.5) flip(prctile(U_Ag_Ex,97.5))],CCtest_Ex,'LineStyle','none','Facealpha',0.2); 
    plot(t_Ex,MLE_Ag_Ex,'-','color',CCtest_Ex,'LineWidth',2);
    
    plot(ts.*ones(101,1),linspace(0,1,101),'-.','LineWidth',1.5,'color',[0.75 0.75 0.75]);
    set(gca,'LineWidth',2,'tickdir','out','Fontsize',18,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:0.1:1],'Ylim',[0 1]);
    ylabel({'Diagnostic sensitivity'},'Fontsize',18);
    text(-7.140724946695094,0.988,'A','Fontsize',30,'FontWeight','bold');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PPA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax1 = axes('Position',[0.227272727272727,0.81661600810537,0.227224980901452,0.172239108409349]);
    patch([t flip(t)],[prctile(U_PPA,2.5) flip(prctile(U_PPA,97.5))],CCtest,'LineStyle','none','Facealpha',0.2); hold on;
    p1=plot(t,MLE_PPA,'-','color',CCtest,'LineWidth',2); hold on;
    p2=scatter(Dt(w==1),100.*truepos(w==1)./totalpos(w==1),40,SymP{1},'filled','MarkerEdgeColor',CCtest,'MarkerFaceColor',CCtest);
    scatter(Dt(~isnan(w) & w<1),100.*truepos(~isnan(w) & w<1)./totalpos(~isnan(w) & w<1),40,SymP{1},'LineWidth',2,'MarkerEdgeColor',CCtest);
    
    patch([t_Ex flip(t_Ex)],[prctile(U_PPA_Ex,2.5) flip(prctile(U_PPA_Ex,97.5))],CCtest_Ex,'LineStyle','none','Facealpha',0.2); hold on;
    p1=plot(t_Ex,MLE_PPA_Ex,'-','color',CCtest_Ex,'LineWidth',2); hold on;
    p2=scatter(Dt_Ex(w_Ex==1),100.*truepos_Ex(w_Ex==1)./totalpos_Ex(w_Ex==1),40,SymP_Ex{1},'filled','MarkerEdgeColor',CCtest_Ex,'MarkerFaceColor',CCtest_Ex);
    scatter(Dt_Ex(~isnan(w_Ex) & w_Ex<1),100.*truepos_Ex(~isnan(w_Ex) & w_Ex<1)./totalpos_Ex(~isnan(w_Ex) & w_Ex<1),40,SymP_Ex{1},'LineWidth',2,'MarkerEdgeColor',CCtest_Ex);

    set(ax1,'LineWidth',1.1,'tickdir','out','Fontsize',12,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:20:100],'Ylim',[0 100]);
    ytickformat('percentage')
    xlabel('Days since symptom onset','Fontsize',12);
    ylabel({'Percent','positive agreement'},'Fontsize',12);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    % Post-quarantine transmission
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    [q,RRTPCR,RAgX,RAgEX] = PQT_RATests(testNamev{ii,1},'DeltaVOC',pA);
    [q_Ex,RRTPCR_Ex,RAgX_Ex,RAgEX_Ex] = PQT_RATests(testNamev{ii,2},'DeltaVOC',pA);
    subplot('Position',[0.6031,0.725,0.381,0.265]);
    
    plot(q,sqrt(RAgX),SymP{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CCtest,'color',CCtest,'Markersize',10); hold on;
    plot(q,sqrt(RAgEX),SymP{1},'LineStyle','-.','LineWidth',2,'MarkerFacecolor','none','color',CCtest,'Markersize',10);  
    
    plot(q_Ex,sqrt(RAgX_Ex),SymP_Ex{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CCtest_Ex,'color',CCtest_Ex,'Markersize',10); 
    plot(q_Ex,sqrt(RAgEX_Ex),SymP_Ex{1},'LineStyle','-.','LineWidth',2,'MarkerFacecolor','none','color',CCtest_Ex,'Markersize',10);  hold off; 
    
    yytick=[0 0.005 0.025 0.05 0.1:0.05:0.3];
    set(gca,'LineWidth',2,'Tickdir','out','Fontsize',18,'XTick',[1:14],'Xlim',[1 14],'Xminortick','off','YTick',sqrt(yytick),'ylim',[0 sqrt(yytick(end))],'YTickLabel',num2str(yytick'),'Yminortick','on');
    ylabel({'Probability of','post-quarantine','transmission'},'Fontsize',18);
    box off;
    grid on;

    legend({[NameP '(FDA): Exit'],[NameP '(FDA): Entry & Exit'],[NameP '(External): Exit'],[NameP '(External): Entry & Exit']},'Fontsize',16);
    text(-3.605,0.54352544978482,'B','Fontsize',30,'FontWeight','bold');
    rmpath([pwd '\Delta_Variant']);
    rmpath([pwd '\Delta_Variant\Results']);

    addpath([pwd '\Alternative_Curve_Delta_Variant']);
    addpath([pwd '\Alternative_Curve_Delta_Variant\Results']);
    [~,~,~,ts,~] = BaselineParameters;
    [~,MLE_Ag,~,U_Ag,~,~,~,~,~,~,t,~,~,~,~] = Sensitivity_for_Plotting(testNamev{ii,1},ts);
    [~,MLE_Ag_Ex,~,U_Ag_Ex,~,~,~,~,~,~,t_Ex,~,~,~,~] = Sensitivity_for_Plotting(testNamev{ii,2},ts);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Diagnostic
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    subplot('Position',[0.074727272727273,0.407,0.381,0.265]);
    p2=patch([t flip(t)],[prctile(U_Ag,2.5) flip(prctile(U_Ag,97.5))],CCtest,'LineStyle','none','Facealpha',0.2); hold on;
    p4=patch([t_Ex flip(t_Ex)],[prctile(U_Ag_Ex,2.5) flip(prctile(U_Ag_Ex,97.5))],CCtest_Ex,'LineStyle','none','Facealpha',0.2);
    
    p1=plot(t,MLE_Ag,'-','color',CCtest,'LineWidth',2);
    p3=plot(t_Ex,MLE_Ag_Ex,'-','color',CCtest_Ex,'LineWidth',2);
    plot(ts.*ones(101,1),linspace(0,1,101),'-.','LineWidth',1.5,'color',[0.75 0.75 0.75]);
    set(gca,'LineWidth',2,'tickdir','out','Fontsize',18,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:0.1:1],'Ylim',[0 1]);
    ylabel({'Diagnostic sensitivity'},'Fontsize',18);

    legend([p1 p2 p3 p4],{['FDA: MLE'],['FDA: 95% CrI'],['External: MLE'],['External: 95% CrI']},'Fontsize',16);
    text(-7.140724946695094,0.988,'C','Fontsize',30,'FontWeight','bold');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    % Post-quarantine transmission
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    [q,~,RAgX,RAgEX] = PQT_RATests(testNamev{ii,1},'DeltaVOC_Alternative_PCR',pA);
    [q_Ex,~,RAgX_Ex,RAgEX_Ex] = PQT_RATests(testNamev{ii,2},'DeltaVOC_Alternative_PCR',pA);
    subplot('Position',[0.6031,0.407,0.381,0.265]);
    plot(q,sqrt(RAgX),SymP{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CCtest,'color',CCtest,'Markersize',10); hold on;
    plot(q,sqrt(RAgEX),SymP{1},'LineStyle','-.','LineWidth',2,'MarkerFacecolor','none','color',CCtest,'Markersize',10);  
    
    plot(q_Ex,sqrt(RAgX_Ex),SymP_Ex{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CCtest_Ex,'color',CCtest_Ex,'Markersize',10); hold on;
    plot(q_Ex,sqrt(RAgEX_Ex),SymP_Ex{1},'LineStyle','-.','LineWidth',2,'MarkerFacecolor','none','color',CCtest_Ex,'Markersize',10);  hold off; 
    
    yytick=[0 0.005 0.025 0.05 0.1:0.05:0.3];
    set(gca,'LineWidth',2,'Tickdir','out','Fontsize',18,'XTick',[1:14],'Xlim',[1 14],'Xminortick','off','YTick',sqrt(yytick),'ylim',[0 sqrt(yytick(end))],'YTickLabel',num2str(yytick'),'Yminortick','on');
    ylabel({'Probability of','post-quarantine','transmission'},'Fontsize',18);
    box off;
    grid on;    

    text(-3.605,0.54352544978482,'D','Fontsize',30,'FontWeight','bold');
    rmpath([pwd '\Alternative_Curve_Delta_Variant']);
    rmpath([pwd '\Alternative_Curve_Delta_Variant\Results']);


    addpath([pwd '\Non_Delta']);
    addpath([pwd '\Non_Delta\Results']);
    [~,~,~,ts,~] = BaselineParameters;
    [~,MLE_Ag,~,U_Ag,~,~,~,~,~,~,t,~,~,~,~] = Sensitivity_for_Plotting(testNamev{ii,1},ts);
    [~,MLE_Ag_Ex,~,U_Ag_Ex,~,~,~,~,~,~,t_Ex,~,~,~,~] = Sensitivity_for_Plotting(testNamev{ii,2},ts);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Diagnostic
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    subplot('Position',[0.074727272727273,0.087,0.381,0.265]);
    patch([t flip(t)],[prctile(U_Ag,2.5) flip(prctile(U_Ag,97.5))],CCtest,'LineStyle','none','Facealpha',0.2); hold on;
    patch([t_Ex flip(t_Ex)],[prctile(U_Ag_Ex,2.5) flip(prctile(U_Ag_Ex,97.5))],CCtest_Ex,'LineStyle','none','Facealpha',0.2); 
    plot(t,MLE_Ag,'-','color',CCtest,'LineWidth',2);
    plot(t_Ex,MLE_Ag_Ex,'-','color',CCtest_Ex,'LineWidth',2);
    plot(ts.*ones(101,1),linspace(0,1,101),'-.','LineWidth',1.5,'color',[0.75 0.75 0.75]);
    set(gca,'LineWidth',2,'tickdir','out','Fontsize',18,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:0.1:1],'Ylim',[0 1]);
    ylabel({'Diagnostic sensitivity'},'Fontsize',18);
    xlabel({'Days post-infection'},'Fontsize',18);
    text(-7.140724946695094,0.988,'E','Fontsize',30,'FontWeight','bold');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    % Post-quarantine transmission
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    [q,~,RAgX,RAgEX] = PQT_RATests(testNamev{ii,1},'General',pA);
    [q_Ex,~,RAgX_Ex,RAgEX_Ex] = PQT_RATests(testNamev{ii,2},'General',pA);
    subplot('Position',[0.6031,0.087,0.381,0.265]);
    plot(q,sqrt(RAgX),SymP{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CCtest,'color',CCtest,'Markersize',10); hold on;
    plot(q,sqrt(RAgEX),SymP{1},'LineStyle','-.','LineWidth',2,'MarkerFacecolor','none','color',CCtest,'Markersize',10);  
    
    plot(q_Ex,sqrt(RAgX_Ex),SymP_Ex{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CCtest_Ex,'color',CCtest_Ex,'Markersize',10); hold on;
    plot(q_Ex,sqrt(RAgEX_Ex),SymP_Ex{1},'LineStyle','-.','LineWidth',2,'MarkerFacecolor','none','color',CCtest_Ex,'Markersize',10);  hold off; 
    
    yytick=[0 0.005 0.025 0.05 0.1:0.05:0.3];
    set(gca,'LineWidth',2,'Tickdir','out','Fontsize',18,'XTick',[1:14],'Xlim',[1 14],'Xminortick','off','YTick',sqrt(yytick),'ylim',[0 sqrt(yytick(end))],'YTickLabel',num2str(yytick'),'Yminortick','on');
    ylabel({'Probability of','post-quarantine','transmission'},'Fontsize',18);
    box off;
    grid on;    

    text(-3.605,0.54352544978482,'F','Fontsize',30,'FontWeight','bold');
    
    rmpath([pwd '\Non_Delta']);
    rmpath([pwd '\Non_Delta\Results']);
    
    print(gcf,['SI_Figure_' NameP '_Compare_DataSet.png'],'-dpng','-r600');
end