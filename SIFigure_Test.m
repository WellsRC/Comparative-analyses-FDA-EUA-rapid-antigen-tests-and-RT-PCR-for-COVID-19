function SIFigure_Test(testName,ts,tL,testNameAll)

figure('units','normalized','outerposition',[0 0 1 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%% Hellewell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diagnostic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
subplot('Position',[0.052,0.555,0.412,0.39]);

load('MLE-Estimate-RTPCR-Hill_Incubation_8_29_days.mat','beta');
betaRTPCR=beta;
t=linspace(0,40,1001);
if(strcmp(testName,'BinaxNOW')) 
    
    load(['BinaxNOW (FDA)_LR_Parameters.mat'],'beta','Dt','totalpos','truepos','w');
    S = TestSensitivity(t,ts,tL,inf,beta,betaRTPCR);  
    
    
    SR = TestSensitivity(t,ts,tL,inf,[],betaRTPCR);
    [CCtestRTPCR,~]=ColourTests('RTPCR');
    [CCtest,SymP]= ColourTests('BinaxNOW (FDA)'); %
    CCtest=CCtest(2,:);
    plot(t,SR,'-','color',CCtestRTPCR,'LineWidth',2); hold on;
    plot(t,S,'-','color',CCtest,'LineWidth',2); hold on;
    plot(ts.*ones(101,1),linspace(0,1,101),'-.','color',[0.3 0.3 0.3],'LineWidth',2);
    % legend([p1,p2],{'RT-PCR','Antigen test'},'Fontsize',20,'Position',[0.099471830985916,0.805854935548247,0.172241778397155,0.127027023482967]);
    % legend boxoff;
    box off;
    set(gca,'LineWidth',2,'tickdir','out','Fontsize',22,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:0.1:1],'Ylim',[0 1]);
    xlabel('Days post-infection','Fontsize',26);
    ylabel({'Diagnostic sensitivty'},'Fontsize',26);
    text(-4.745,1,'A','Fontsize',32,'FontWeight','bold');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PPA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax1 = axes('Position',[0.250525210084034,0.711246200607903,0.206407563025209,0.232016210739622]);
    p1=plot(t,100.*LR(t,beta),'-','color',CCtest,'LineWidth',2); hold on;
    p2=scatter(Dt(w==1),100.*truepos(w==1)./totalpos(w==1),40,SymP{1},'filled','MarkerEdgeColor',CCtest,'MarkerFaceColor',CCtest);
    scatter(Dt(~isnan(w) & w<1),100.*truepos(~isnan(w) & w<1)./totalpos(~isnan(w) & w<1),40,SymP{1},'LineWidth',2,'MarkerEdgeColor',CCtest);
        
    box off;
    set(gca,'LineWidth',1.1,'tickdir','out','Fontsize',18,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:20:100],'Ylim',[0 100]);
    ytickformat('percentage')
    xlabel('Days since symptom onset','Fontsize',18);
    ylabel({'Percent','positive agreement'},'Fontsize',18);

    if(100.*LR(15.5,beta)<77)
        legend([p1,p2],{'Logistic regression','Percent positive agreement'},'Fontsize',14,'LineWidth',1.1,'Position',[0.330357142857142,0.89225189765564,0.147548351256234,0.062513278830154]);
    else
        legend([p1,p2],{'Logistic regression','Percent positive agreement'},'Fontsize',14,'LineWidth',1.1,'Position',[0.330357142857142,max(LR(40,beta),LR(40,beta2)).*0.141+0.714,0.147548351256234,0.059270515086803]);
    end
elseif(strcmp(testName,'BinaxNOW (Community)'))
    fff=strcmp(testNameAll,'BinaxNOW (FDA)');
    load([testNameAll{fff} '_LR_Parameters.mat'],'beta','Dt','totalpos','truepos','w');    
    S2 = TestSensitivity(t,ts,tL,inf,beta,betaRTPCR);
    [CCtest2,SymP2]= ColourTests(testNameAll{fff});
    CCtest2=CCtest2(1,:);
    beta2=beta;
    Dt2=Dt;
    totalpos2=totalpos;
    truepos2=truepos;
    w2=w;
    
    load([testName '_LR_Parameters.mat'],'beta','Dt','totalpos','truepos','w');
    S = TestSensitivity(t,ts,tL,inf,beta,betaRTPCR);  
    
    
    SR = TestSensitivity(t,ts,tL,inf,[],betaRTPCR);
    [CCtestRTPCR,~]=ColourTests('RTPCR');
    [CCtest,SymP]= ColourTests(testName);
    plot(t,SR,'-','color',CCtestRTPCR,'LineWidth',2); hold on;
    plot(t,S2,'-','color',CCtest2,'LineWidth',2); hold on;
    plot(t,S,'-.','color',CCtest,'LineWidth',2); hold on;
    plot(ts.*ones(101,1),linspace(0,1,101),'-.','color',[0.3 0.3 0.3],'LineWidth',2);
    % legend([p1,p2],{'RT-PCR','Antigen test'},'Fontsize',20,'Position',[0.099471830985916,0.805854935548247,0.172241778397155,0.127027023482967]);
    % legend boxoff;
    box off;
    set(gca,'LineWidth',2,'tickdir','out','Fontsize',22,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:0.1:1],'Ylim',[0 1]);
    xlabel('Days post-infection','Fontsize',26);
    ylabel({'Diagnostic sensitivty'},'Fontsize',26);

    text(-4.745,1,'A','Fontsize',32,'FontWeight','bold');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PPA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax1 = axes('Position',[0.250525210084034,0.711246200607903,0.206407563025209,0.232016210739622]);
    p1=plot(t,100.*LR(t,beta),'-.','color',CCtest,'LineWidth',2); hold on;
    p2=scatter(Dt(w==1),100.*truepos(w==1)./totalpos(w==1),40,SymP{1},'filled','MarkerEdgeColor',CCtest,'MarkerFaceColor',CCtest);
    scatter(Dt(~isnan(w) & w<1),100.*truepos(~isnan(w) & w<1)./totalpos(~isnan(w) & w<1),40,SymP{1},'LineWidth',2,'MarkerEdgeColor',CCtest);
    
    
    plot(t,100.*LR(t,beta2),'-','color',CCtest2,'LineWidth',2); hold on;
    scatter(Dt2(w2==1),100.*truepos2(w2==1)./totalpos2(w2==1),40,SymP2{1},'filled','MarkerEdgeColor',CCtest2,'MarkerFaceColor',CCtest2);
    scatter(Dt2(~isnan(w2) & w2<1),100.*truepos2(~isnan(w2) & w2<1)./totalpos2(~isnan(w2) & w2<1),40,SymP2{1},'LineWidth',2,'MarkerEdgeColor',CCtest2);
    
    box off;
    set(gca,'LineWidth',1.1,'tickdir','out','Fontsize',18,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:20:100],'Ylim',[0 100]);
    ytickformat('percentage')
    xlabel('Days since symptom onset','Fontsize',18);
    ylabel({'Percent','positive agreement'},'Fontsize',18);

        legend([p1,p2],{'Logistic regression','Percent positive agreement'},'Fontsize',14,'LineWidth',1.1,'Position',[0.323529411764705,0.861244417865645,0.147548351256234,0.062513278830154]);
elseif(strcmp(testName,'CareStart (Anterior Nasal Swab)')) 
    fff=strcmp(testNameAll,'CareStart (Nasopharyngeal Swab)');
    load([testNameAll{fff} '_LR_Parameters.mat'],'beta','Dt','totalpos','truepos','w');    
    S2 = TestSensitivity(t,ts,tL,inf,beta,betaRTPCR);
    [CCtest2,SymP2]= ColourTests(testNameAll{fff});
    
    beta2=beta;
    Dt2=Dt;
    totalpos2=totalpos;
    truepos2=truepos;
    w2=w;
    
    load(['CareStart (Anterior Nasal Swab - FDA)_LR_Parameters.mat'],'beta','Dt','totalpos','truepos','w');
    S = TestSensitivity(t,ts,tL,inf,beta,betaRTPCR);  
    
    
    SR = TestSensitivity(t,ts,tL,inf,[],betaRTPCR);
    [CCtestRTPCR,~]=ColourTests('RTPCR');
    [CCtest,SymP]= ColourTests('CareStart (Anterior Nasal Swab - FDA)'); %
    CCtest=CCtest(2,:);
    plot(t,SR,'-','color',CCtestRTPCR,'LineWidth',2); hold on;
    plot(t,S,'-','color',CCtest,'LineWidth',2); hold on;
    plot(t,S2,'-.','color',CCtest2,'LineWidth',2); hold on;
    plot(ts.*ones(101,1),linspace(0,1,101),'-.','color',[0.3 0.3 0.3],'LineWidth',2);
    % legend([p1,p2],{'RT-PCR','Antigen test'},'Fontsize',20,'Position',[0.099471830985916,0.805854935548247,0.172241778397155,0.127027023482967]);
    % legend boxoff;
    box off;
    set(gca,'LineWidth',2,'tickdir','out','Fontsize',22,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:0.1:1],'Ylim',[0 1]);
    xlabel('Days post-infection','Fontsize',26);
    ylabel({'Diagnostic sensitivty'},'Fontsize',26);
    text(-4.745,1,'A','Fontsize',32,'FontWeight','bold');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PPA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax1 = axes('Position',[0.250525210084034,0.711246200607903,0.206407563025209,0.232016210739622]);
    p1=plot(t,100.*LR(t,beta),'-','color',CCtest,'LineWidth',2); hold on;
    p2=scatter(Dt(w==1),100.*truepos(w==1)./totalpos(w==1),40,SymP{1},'filled','MarkerEdgeColor',CCtest,'MarkerFaceColor',CCtest);
    scatter(Dt(~isnan(w) & w<1),100.*truepos(~isnan(w) & w<1)./totalpos(~isnan(w) & w<1),40,SymP{1},'LineWidth',2,'MarkerEdgeColor',CCtest);
    
    
    plot(t,100.*LR(t,beta2),'-.','color',CCtest2,'LineWidth',2); hold on;
    scatter(Dt2(w2==1),100.*truepos2(w2==1)./totalpos2(w2==1),40,SymP2{1},'filled','MarkerEdgeColor',CCtest2,'MarkerFaceColor',CCtest2);
    scatter(Dt2(~isnan(w2) & w2<1),100.*truepos2(~isnan(w2) & w2<1)./totalpos2(~isnan(w2) & w2<1),40,SymP2{1},'LineWidth',2,'MarkerEdgeColor',CCtest2);
    
    box off;
    set(gca,'LineWidth',1.1,'tickdir','out','Fontsize',18,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:20:100],'Ylim',[0 100]);
    ytickformat('percentage')
    xlabel('Days since symptom onset','Fontsize',18);
    ylabel({'Percent','positive agreement'},'Fontsize',18);

    if(100.*max(LR(15.5,beta),LR(15.5,beta2))<77)
        legend([p1,p2],{'Logistic regression','Percent positive agreement'},'Fontsize',14,'LineWidth',1.1,'Position',[0.330357142857142,0.89225189765564,0.147548351256234,0.062513278830154]);
    else
        legend([p1,p2],{'Logistic regression','Percent positive agreement'},'Fontsize',14,'LineWidth',1.1,'Position',[0.330357142857142,max(LR(40,beta),LR(40,beta2)).*0.141+0.714,0.147548351256234,0.059270515086803]);
    end
elseif(strcmp(testName,'CareStart (Anterior Nasal Swab - FDA)')) 
    fff=strcmp(testNameAll,'CareStart (Anterior Nasal Swab - External)');
    load([testNameAll{fff} '_LR_Parameters.mat'],'beta','Dt','totalpos','truepos','w');    
    S2 = TestSensitivity(t,ts,tL,inf,beta,betaRTPCR);
    [CCtest2,SymP2]= ColourTests(testNameAll{fff});
    
    beta2=beta;
    Dt2=Dt;
    totalpos2=totalpos;
    truepos2=truepos;
    w2=w;
    
    load([testName '_LR_Parameters.mat'],'beta','Dt','totalpos','truepos','w');
    S = TestSensitivity(t,ts,tL,inf,beta,betaRTPCR);  
    
    
    SR = TestSensitivity(t,ts,tL,inf,[],betaRTPCR);
    [CCtestRTPCR,~]=ColourTests('RTPCR');
    [CCtest,SymP]= ColourTests(testName);
    CCtest=CCtest(1,:);
    plot(t,SR,'-','color',CCtestRTPCR,'LineWidth',2); hold on;
    plot(t,S,'-','color',CCtest,'LineWidth',2); hold on;
    plot(t,S2,'-.','color',CCtest2,'LineWidth',2); hold on;
    plot(ts.*ones(101,1),linspace(0,1,101),'-.','color',[0.3 0.3 0.3],'LineWidth',2);
    % legend([p1,p2],{'RT-PCR','Antigen test'},'Fontsize',20,'Position',[0.099471830985916,0.805854935548247,0.172241778397155,0.127027023482967]);
    % legend boxoff;
    box off;
    set(gca,'LineWidth',2,'tickdir','out','Fontsize',22,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:0.1:1],'Ylim',[0 1]);
    xlabel('Days post-infection','Fontsize',26);
    ylabel({'Diagnostic sensitivty'},'Fontsize',26);
    text(-4.745,1,'A','Fontsize',32,'FontWeight','bold');


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PPA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax1 = axes('Position',[0.250525210084034,0.711246200607903,0.206407563025209,0.232016210739622]);
    p1=plot(t,100.*LR(t,beta),'-','color',CCtest,'LineWidth',2); hold on;
    p2=scatter(Dt(w==1),100.*truepos(w==1)./totalpos(w==1),40,SymP{1},'filled','MarkerEdgeColor',CCtest,'MarkerFaceColor',CCtest);
    scatter(Dt(~isnan(w) & w<1),100.*truepos(~isnan(w) & w<1)./totalpos(~isnan(w) & w<1),40,SymP{1},'LineWidth',2,'MarkerEdgeColor',CCtest);
    
    
    plot(t,100.*LR(t,beta2),'-.','color',CCtest2,'LineWidth',2); hold on;
    scatter(Dt2(w2==1),100.*truepos2(w2==1)./totalpos2(w2==1),40,SymP2{1},'filled','MarkerEdgeColor',CCtest2,'MarkerFaceColor',CCtest2);
    scatter(Dt2(~isnan(w2) & w2<1),100.*truepos2(~isnan(w2) & w2<1)./totalpos2(~isnan(w2) & w2<1),40,SymP2{1},'LineWidth',2,'MarkerEdgeColor',CCtest2);
    
    box off;
    set(gca,'LineWidth',1.1,'tickdir','out','Fontsize',18,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:20:100],'Ylim',[0 100]);
    ytickformat('percentage')
    xlabel('Days since symptom onset','Fontsize',18);
    ylabel({'Percent','positive agreement'},'Fontsize',18);
    legend([p1,p2],{'Logistic regression','Percent positive agreement'},'Fontsize',14,'LineWidth',1.1,'Position',[0.331407563025209,0.843752648621009,0.147548351256234,0.062513278830154]);
elseif(strcmp(testName,'Sofia (CDC)'))  
    fff=strcmp(testNameAll,'Sofia (FDA)');
    load([testNameAll{fff} '_LR_Parameters.mat'],'beta','Dt','totalpos','truepos','w');    
    S2 = TestSensitivity(t,ts,tL,inf,beta,betaRTPCR);
    [CCtest2,SymP2]= ColourTests(testNameAll{fff});
    CCtest2=CCtest2(1,:);
    beta2=beta;
    Dt2=Dt;
    totalpos2=totalpos;
    truepos2=truepos;
    w2=w;
    
    load([testName '_LR_Parameters.mat'],'beta','Dt','totalpos','truepos','w');
    S = TestSensitivity(t,ts,tL,inf,beta,betaRTPCR);  
    
    
    SR = TestSensitivity(t,ts,tL,inf,[],betaRTPCR);
    [CCtestRTPCR,~]=ColourTests('RTPCR');
    [CCtest,SymP]= ColourTests(testName);
    plot(t,SR,'-','color',CCtestRTPCR,'LineWidth',2); hold on;
    plot(t,S2,'-','color',CCtest2,'LineWidth',2); hold on;
    plot(t,S,'-.','color',CCtest,'LineWidth',2); hold on;
    plot(ts.*ones(101,1),linspace(0,1,101),'-.','color',[0.3 0.3 0.3],'LineWidth',2);
    % legend([p1,p2],{'RT-PCR','Antigen test'},'Fontsize',20,'Position',[0.099471830985916,0.805854935548247,0.172241778397155,0.127027023482967]);
    % legend boxoff;
    box off;
    set(gca,'LineWidth',2,'tickdir','out','Fontsize',22,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:0.1:1],'Ylim',[0 1]);
    xlabel('Days post-infection','Fontsize',26);
    ylabel({'Diagnostic sensitivty'},'Fontsize',26);

    text(-4.745,1,'A','Fontsize',32,'FontWeight','bold');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PPA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax1 = axes('Position',[0.250525210084034,0.711246200607903,0.206407563025209,0.232016210739622]);
    p1=plot(t,100.*LR(t,beta),'-.','color',CCtest,'LineWidth',2); hold on;
    p2=scatter(Dt(w==1),100.*truepos(w==1)./totalpos(w==1),40,SymP{1},'filled','MarkerEdgeColor',CCtest,'MarkerFaceColor',CCtest);
    scatter(Dt(~isnan(w) & w<1),100.*truepos(~isnan(w) & w<1)./totalpos(~isnan(w) & w<1),40,SymP{1},'LineWidth',2,'MarkerEdgeColor',CCtest);
    
    
    plot(t,100.*LR(t,beta2),'-','color',CCtest2,'LineWidth',2); hold on;
    scatter(Dt2(w2==1),100.*truepos2(w2==1)./totalpos2(w2==1),40,SymP2{1},'filled','MarkerEdgeColor',CCtest2,'MarkerFaceColor',CCtest2);
    scatter(Dt2(~isnan(w2) & w2<1),100.*truepos2(~isnan(w2) & w2<1)./totalpos2(~isnan(w2) & w2<1),40,SymP2{1},'LineWidth',2,'MarkerEdgeColor',CCtest2);
    
    box off;
    set(gca,'LineWidth',1.1,'tickdir','out','Fontsize',18,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:20:100],'Ylim',[0 100]);
    ytickformat('percentage')
    xlabel('Days since symptom onset','Fontsize',18);
    ylabel({'Percent','positive agreement'},'Fontsize',18);

    if(100.*max(LR(15.5,beta),LR(15.5,beta2))<77)
        legend([p1,p2],{'Logistic regression','Percent positive agreement'},'Fontsize',14,'LineWidth',1.1,'Position',[0.330357142857142,0.89225189765564,0.147548351256234,0.062513278830154]);
    else
        legend([p1,p2],{'Logistic regression','Percent positive agreement'},'Fontsize',14,'LineWidth',1.1,'Position',[0.330357142857142,min(LR(40,beta),LR(40,beta2)).*0.141+0.714,0.147548351256234,0.059270515086803]);
    end
elseif(strcmp(testName,'Liaison (Anterior Nasal Swab)')) 
    fff=strcmp(testNameAll,'Liaison (Nasalpharyngeal Swab)');
    load([testNameAll{fff} '_LR_Parameters.mat'],'beta','Dt','totalpos','truepos','w');    
    S2 = TestSensitivity(t,ts,tL,inf,beta,betaRTPCR);
    [CCtest2,SymP2]= ColourTests(testNameAll{fff});
    
    beta2=beta;
    Dt2=Dt;
    totalpos2=totalpos;
    truepos2=truepos;
    w2=w;
    
    load([testName '_LR_Parameters.mat'],'beta','Dt','totalpos','truepos','w');
    S = TestSensitivity(t,ts,tL,inf,beta,betaRTPCR);  
    
    
    SR = TestSensitivity(t,ts,tL,inf,[],betaRTPCR);
    [CCtestRTPCR,~]=ColourTests('RTPCR');
    [CCtest,SymP]= ColourTests(testName);
    plot(t,SR,'-','color',CCtestRTPCR,'LineWidth',2); hold on;
    plot(t,S,'-','color',CCtest,'LineWidth',2); hold on;
    plot(t,S2,'-.','color',CCtest2,'LineWidth',2); hold on;
    plot(ts.*ones(101,1),linspace(0,1,101),'-.','color',[0.3 0.3 0.3],'LineWidth',2);
    % legend([p1,p2],{'RT-PCR','Antigen test'},'Fontsize',20,'Position',[0.099471830985916,0.805854935548247,0.172241778397155,0.127027023482967]);
    % legend boxoff;
    box off;
    set(gca,'LineWidth',2,'tickdir','out','Fontsize',22,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:0.1:1],'Ylim',[0 1]);
    xlabel('Days post-infection','Fontsize',26);
    ylabel({'Diagnostic sensitivty'},'Fontsize',26);

    text(-4.745,1,'A','Fontsize',32,'FontWeight','bold');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PPA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax1 = axes('Position',[0.250525210084034,0.711246200607903,0.206407563025209,0.232016210739622]);
    p1=plot(t,100.*LR(t,beta),'-','color',CCtest,'LineWidth',2); hold on;
    p2=scatter(Dt(w==1),100.*truepos(w==1)./totalpos(w==1),40,SymP{1},'filled','MarkerEdgeColor',CCtest,'MarkerFaceColor',CCtest);
    scatter(Dt(~isnan(w) & w<1),100.*truepos(~isnan(w) & w<1)./totalpos(~isnan(w) & w<1),40,SymP{1},'LineWidth',2,'MarkerEdgeColor',CCtest);
    
    
    plot(t,100.*LR(t,beta2),'-.','color',CCtest2,'LineWidth',2); hold on;
    scatter(Dt2(w2==1),100.*truepos2(w2==1)./totalpos2(w2==1),40,SymP2{1},'filled','MarkerEdgeColor',CCtest2,'MarkerFaceColor',CCtest2);
    scatter(Dt2(~isnan(w2) & w2<1),100.*truepos2(~isnan(w2) & w2<1)./totalpos2(~isnan(w2) & w2<1),40,SymP2{1},'LineWidth',2,'MarkerEdgeColor',CCtest2);
    
    box off;
    set(gca,'LineWidth',1.1,'tickdir','out','Fontsize',18,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:20:100],'Ylim',[0 100]);
    ytickformat('percentage')
    xlabel('Days since symptom onset','Fontsize',18);
    ylabel({'Percent','positive agreement'},'Fontsize',18);
    if(100.*max(LR(15.5,beta),LR(15.5,beta2))<77)
        legend([p1,p2],{'Logistic regression','Percent positive agreement'},'Fontsize',14,'LineWidth',1.1,'Position',[0.330357142857142,0.89225189765564,0.147548351256234,0.062513278830154]);
    else
        legend([p1,p2],{'Logistic regression','Percent positive agreement'},'Fontsize',14,'LineWidth',1.1,'Position',[0.330357142857142,min(LR(40,beta),LR(40,beta2)).*0.141+0.714,0.147548351256234,0.059270515086803]);
    end
elseif(strcmp(testName,'LumiraDX (Anterior Nasal Swab)')) 
    fff=strcmp(testNameAll,'LumiraDX (Nasopharyngeal Swabs)');
    load([testNameAll{fff} '_LR_Parameters.mat'],'beta','Dt','totalpos','truepos','w');    
    S2 = TestSensitivity(t,ts,tL,inf,beta,betaRTPCR);
    [CCtest2,SymP2]= ColourTests(testNameAll{fff});
    
    beta2=beta;
    Dt2=Dt;
    totalpos2=totalpos;
    truepos2=truepos;
    w2=w;
    
    load([testName '_LR_Parameters.mat'],'beta','Dt','totalpos','truepos','w');
    S = TestSensitivity(t,ts,tL,inf,beta,betaRTPCR);  
    
    
    SR = TestSensitivity(t,ts,tL,inf,[],betaRTPCR);
    [CCtestRTPCR,~]=ColourTests('RTPCR');
    [CCtest,SymP]= ColourTests(testName);
    plot(t,SR,'-','color',CCtestRTPCR,'LineWidth',2); hold on;
    plot(t,S,'-','color',CCtest,'LineWidth',2); hold on;
    plot(t,S2,'-.','color',CCtest2,'LineWidth',2); hold on;
    plot(ts.*ones(101,1),linspace(0,1,101),'-.','color',[0.3 0.3 0.3],'LineWidth',2);
    % legend([p1,p2],{'RT-PCR','Antigen test'},'Fontsize',20,'Position',[0.099471830985916,0.805854935548247,0.172241778397155,0.127027023482967]);
    % legend boxoff;
    box off;
    set(gca,'LineWidth',2,'tickdir','out','Fontsize',22,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:0.1:1],'Ylim',[0 1]);
    xlabel('Days post-infection','Fontsize',26);
    ylabel({'Diagnostic sensitivty'},'Fontsize',26);

    
    text(-4.745,1,'A','Fontsize',32,'FontWeight','bold');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PPA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax1 = axes('Position',[0.250525210084034,0.711246200607903,0.206407563025209,0.232016210739622]);
    p1=plot(t,100.*LR(t,beta),'-','color',CCtest,'LineWidth',2); hold on;
    p2=scatter(Dt(w==1),100.*truepos(w==1)./totalpos(w==1),40,SymP{1},'filled','MarkerEdgeColor',CCtest,'MarkerFaceColor',CCtest);
    scatter(Dt(~isnan(w) & w<1),100.*truepos(~isnan(w) & w<1)./totalpos(~isnan(w) & w<1),40,SymP{1},'LineWidth',2,'MarkerEdgeColor',CCtest);
    
    
    plot(t,100.*LR(t,beta2),'-.','color',CCtest2,'LineWidth',2); hold on;
    scatter(Dt2(w2==1),100.*truepos2(w2==1)./totalpos2(w2==1),40,SymP2{1},'filled','MarkerEdgeColor',CCtest2,'MarkerFaceColor',CCtest2);
    scatter(Dt2(~isnan(w2) & w2<1),100.*truepos2(~isnan(w2) & w2<1)./totalpos2(~isnan(w2) & w2<1),40,SymP2{1},'LineWidth',2,'MarkerEdgeColor',CCtest2);
    
    box off;
    set(gca,'LineWidth',1.1,'tickdir','out','Fontsize',18,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:20:100],'Ylim',[0 100]);
    ytickformat('percentage')
    xlabel('Days since symptom onset','Fontsize',18);
    ylabel({'Percent','positive agreement'},'Fontsize',18);

    if(100.*max(LR(15.5,beta),LR(15.5,beta2))<77)
        legend([p1,p2],{'Logistic regression','Percent positive agreement'},'Fontsize',14,'LineWidth',1.1,'Position',[0.330357142857142,0.89225189765564,0.147548351256234,0.062513278830154]);
    else        
        legend([p1,p2],{'Logistic regression','Percent positive agreement'},'Fontsize',14,'LineWidth',1.1,'Position',[0.330357142857142,min(LR(40,beta),LR(40,beta2)).*0.141+0.714,0.147548351256234,0.059270515086803]);
    end
elseif(~strcmp(testName,'Sofia (CDC)')&&~strcmp(testName,'CareStart (Anterior Nasal Swab - FDA)')&&~strcmp(testName,'CareStart (Anterior Nasal Swab - External)')&&~strcmp(testName,'LumiraDX (Anterior Nasal Swab)')&&~strcmp(testName,'LumiraDX (Nasopharyngeal Swabs)')&&~strcmp(testName,'Liaison (Anterior Nasal Swab)')&&~strcmp(testName,'Liaison (Nasalpharyngeal Swab)')&&~strcmp(testName,'CareStart (Anterior Nasal Swab)')&&~strcmp(testName,'CareStart (Nasopharyngeal Swab)'))
    
    load([testName '_LR_Parameters.mat'],'beta','Dt','totalpos','truepos','w');

    S = TestSensitivity(t,ts,tL,inf,beta,betaRTPCR);
    SR = TestSensitivity(t,ts,tL,inf,[],betaRTPCR);
    [CCtestRTPCR,~]=ColourTests('RTPCR');
    [CCtest,SymP]= ColourTests(testName);
    CCtest=CCtest(end,:);
    p1=plot(t,SR,'-','color',CCtestRTPCR,'LineWidth',2); hold on;
    p2=plot(t,S,'-','color',CCtest,'LineWidth',2); hold on;
    plot(ts.*ones(101,1),linspace(0,1,101),'-.','color',[0.3 0.3 0.3],'LineWidth',2);
    % legend([p1,p2],{'RT-PCR','Antigen test'},'Fontsize',20,'Position',[0.099471830985916,0.805854935548247,0.172241778397155,0.127027023482967]);
    % legend boxoff;
    box off;
    set(gca,'LineWidth',2,'tickdir','out','Fontsize',22,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:0.1:1],'Ylim',[0 1]);
    xlabel('Days post-infection','Fontsize',26);
    ylabel({'Diagnostic sensitivty'},'Fontsize',26);

    
    text(-4.745,1,'A','Fontsize',32,'FontWeight','bold');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PPA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax1 = axes('Position',[0.250525210084034,0.711246200607903,0.206407563025209,0.232016210739622]);
    p1=plot(t,100.*LR(t,beta),'-','color',CCtest,'LineWidth',2); hold on;
    p2=scatter(Dt(w==1),100.*truepos(w==1)./totalpos(w==1),40,SymP{1},'filled','MarkerEdgeColor',CCtest,'MarkerFaceColor',CCtest);
    scatter(Dt(~isnan(w) & w<1),100.*truepos(~isnan(w) & w<1)./totalpos(~isnan(w) & w<1),40,SymP{1},'LineWidth',2,'MarkerEdgeColor',CCtest);
    box off;
    set(gca,'LineWidth',1.1,'tickdir','out','Fontsize',18,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:20:100],'Ylim',[0 100]);
    ytickformat('percentage')
    xlabel('Days since symptom onset','Fontsize',18);
    ylabel({'Percent','positive agreement'},'Fontsize',18);
    if(100.*LR(15.5,beta)<77)
        legend([p1,p2],{'Logistic regression','Percent positive agreement'},'Fontsize',14,'LineWidth',1.1,'Position',[0.330357142857142,0.89225189765564,0.147548351256234,0.062513278830154]);
    else
        legend([p1,p2],{'Logistic regression','Percent positive agreement'},'Fontsize',14,'LineWidth',1.1,'Position',[0.330357142857142,LR(40,beta).*0.141+0.714,0.147548351256234,0.059270515086803]);
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Probabilility of PQT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot('Position',[0.5745,0.555,0.412,0.39]);
Risk=1;
load('TestingonExit_RTPCR_24hrDelay_Hellewell.mat','IDSLA','IDSLS','pA','q')
RRTPCR=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);
[CRTPCR,MFRTPCR]= ColourTests('RTPCR');
if(strcmp(testName,'BinaxNOW')) 
    [CAg,MFAg]= ColourTests('BinaxNOW (FDA)'); % As FDA is baseline not the average
    CAg=CAg(2,:);
    load(['TestingonExit_BinaxNOW (FDA)_NoDelay_Hellewell.mat'],'IDSLA','IDSLS');
    RAgX=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    load(['TestingonEntryExit_BinaxNOW (FDA)_NoDelay_Hellewell.mat'],'IDSLA','IDSLS');
    RAgEX=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    plot(q,sqrt(RRTPCR),MFRTPCR{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CRTPCR,'color',CRTPCR,'Markersize',10); hold on;
    plot(q,sqrt(RAgX),MFAg{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CAg,'color',CAg,'Markersize',10); 
    plot(q,sqrt(RAgEX),MFAg{1},'LineStyle','-.','LineWidth',2,'MarkerFacecolor','none','color',CAg,'Markersize',10);  hold off; 

    yytick=[0 0.001 0.01 0.025 0.05 0.075 0.1:0.05:0.3];
    set(gca,'LineWidth',2,'Tickdir','out','Fontsize',20,'XTick',[1:14],'Xlim',[1 14],'Xminortick','off','YTick',sqrt(yytick),'ylim',[0 sqrt(yytick(end))],'YTickLabel',num2str(yytick'),'Yminortick','on');
    xlabel('Duration of quarantine (days)','Fontsize',26);
    ylabel({'Probability of','post-quarantine transmission'},'Fontsize',24);
    legend({'RT-PCR on exit',['BinaxNOW on exit'],['BinaxNOW on entry and exit']},'Fontsize',18);
    legend boxoff;
    box off;
    grid on;
    
    text(-1.205,sqrt(0.36),'BinaxNOW','Fontsize',32,'HorizontalAlignment','center') 
    

    text(-1.764,sqrt(0.3),'B','Fontsize',32,'FontWeight','bold');
elseif(strcmp(testName,'BinaxNOW (Community)')) 
    [CAg,MFAg]= ColourTests(testName);
    load(['TestingonExit_' testName '_NoDelay_Hellewell.mat'],'IDSLA','IDSLS');
    RAgX=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    load(['TestingonEntryExit_' testName '_NoDelay_Hellewell.mat'],'IDSLA','IDSLS');
    RAgEX=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    fff=strcmp(testNameAll,'BinaxNOW (FDA)');
    [CAg2,MFAg2]= ColourTests(testNameAll{fff});
    CAg2=CAg2(1,:);
    load(['TestingonExit_' testNameAll{fff} '_NoDelay_Hellewell.mat'],'IDSLA','IDSLS');
    RAgX2=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    load(['TestingonEntryExit_' testNameAll{fff} '_NoDelay_Hellewell.mat'],'IDSLA','IDSLS');
    RAgEX2=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    plot(q,sqrt(RRTPCR),MFRTPCR{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CRTPCR,'color',CRTPCR,'Markersize',10); hold on;
    plot(q,sqrt(RAgX2),MFAg2{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CAg2,'color',CAg2,'Markersize',10); 
    plot(q,sqrt(RAgEX2),MFAg2{1},'LineStyle','-.','LineWidth',2,'MarkerFacecolor','none','color',CAg2,'Markersize',10);
    plot(q,sqrt(RAgX),MFAg{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CAg,'color',CAg,'Markersize',10); 
    plot(q,sqrt(RAgEX),MFAg{1},'LineStyle','-.','LineWidth',2,'MarkerFacecolor','none','color',CAg,'Markersize',10);  hold off; 

    yytick=[0 0.001 0.01 0.025 0.05 0.075 0.1:0.05:0.3];
    set(gca,'LineWidth',2,'Tickdir','out','Fontsize',20,'XTick',[1:14],'Xlim',[1 14],'Xminortick','off','YTick',sqrt(yytick),'ylim',[0 sqrt(yytick(end))],'YTickLabel',num2str(yytick'),'Yminortick','on');
    xlabel('Duration of quarantine (days)','Fontsize',26);
    ylabel({'Probability of','post-quarantine transmission'},'Fontsize',24);
    legend({'RT-PCR on exit',['BinaxNOW (Internal) on exit'],['BinaxNOW (Internal) on entry and exit'],['BinaxNOW (External) on exit'],['BinaxNOW (External) on entry and exit']},'Fontsize',18,'Position',[0.671918776381958,0.78858494220131,0.324579823057817,0.163120562695684]);
    legend boxoff;
    box off;
    grid on;
    
    text(-1.205,sqrt(0.36),'BinaxNOW','Fontsize',32,'HorizontalAlignment','center') 
elseif(strcmp(testName,'CareStart (Anterior Nasal Swab)')) 
    [CAg,MFAg]= ColourTests('CareStart (Anterior Nasal Swab - FDA)'); % As FDA is baseline not the average
    CAg=CAg(2,:);
    load(['TestingonExit_CareStart (Anterior Nasal Swab - FDA)_NoDelay_Hellewell.mat'],'IDSLA','IDSLS');
    RAgX=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    load(['TestingonEntryExit_CareStart (Anterior Nasal Swab - FDA)_NoDelay_Hellewell.mat'],'IDSLA','IDSLS');
    RAgEX=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    fff=strcmp(testNameAll,'CareStart (Nasopharyngeal Swab)');
    [CAg2,MFAg2]= ColourTests(testNameAll{fff});
    load(['TestingonExit_' testNameAll{fff} '_NoDelay_Hellewell.mat'],'IDSLA','IDSLS');
    RAgX2=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    load(['TestingonEntryExit_' testNameAll{fff} '_NoDelay_Hellewell.mat'],'IDSLA','IDSLS');
    RAgEX2=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    plot(q,sqrt(RRTPCR),MFRTPCR{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CRTPCR,'color',CRTPCR,'Markersize',10); hold on;
    plot(q,sqrt(RAgX),MFAg{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CAg,'color',CAg,'Markersize',10); 
    plot(q,sqrt(RAgEX),MFAg{1},'LineStyle','-.','LineWidth',2,'MarkerFacecolor','none','color',CAg,'Markersize',10); 
    plot(q,sqrt(RAgX2),MFAg2{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CAg2,'color',CAg2,'Markersize',10); 
    plot(q,sqrt(RAgEX2),MFAg2{1},'LineStyle','-.','LineWidth',2,'MarkerFacecolor','none','color',CAg2,'Markersize',10); hold off; 

    yytick=[0 0.001 0.01 0.025 0.05 0.075 0.1:0.05:0.3];
    set(gca,'LineWidth',2,'Tickdir','out','Fontsize',20,'XTick',[1:14],'Xlim',[1 14],'Xminortick','off','YTick',sqrt(yytick),'ylim',[0 sqrt(yytick(end))],'YTickLabel',num2str(yytick'),'Yminortick','on');
    xlabel('Duration of quarantine (days)','Fontsize',26);
    ylabel({'Probability of','post-quarantine transmission'},'Fontsize',24);
    legend({'RT-PCR on exit',[testName ' on exit'],[testName ' on entry and exit'],[testNameAll{fff} ' on exit'],[testNameAll{fff} ' on entry and exit']},'Fontsize',18,'Position',[0.671918776381958,0.78858494220131,0.324579823057817,0.163120562695684]);
    legend boxoff;
    box off;
    grid on;
    
    text(-1.205,sqrt(0.36),'CareStart','Fontsize',32,'HorizontalAlignment','center') 
    

    text(-1.764,sqrt(0.3),'B','Fontsize',32,'FontWeight','bold');
elseif(strcmp(testName,'CareStart (Anterior Nasal Swab - FDA)')) 
    [CAg,MFAg]= ColourTests(testName);
    CAg=CAg(1,:);
    load(['TestingonExit_' testName '_NoDelay_Hellewell.mat'],'IDSLA','IDSLS');
    RAgX=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    load(['TestingonEntryExit_' testName '_NoDelay_Hellewell.mat'],'IDSLA','IDSLS');
    RAgEX=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    fff=strcmp(testNameAll,'CareStart (Anterior Nasal Swab - External)');
    [CAg2,MFAg2]= ColourTests(testNameAll{fff});
    load(['TestingonExit_' testNameAll{fff} '_NoDelay_Hellewell.mat'],'IDSLA','IDSLS');
    RAgX2=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    load(['TestingonEntryExit_' testNameAll{fff} '_NoDelay_Hellewell.mat'],'IDSLA','IDSLS');
    RAgEX2=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    plot(q,sqrt(RRTPCR),MFRTPCR{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CRTPCR,'color',CRTPCR,'Markersize',10); hold on;
    plot(q,sqrt(RAgX),MFAg{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CAg,'color',CAg,'Markersize',10); 
    plot(q,sqrt(RAgEX),MFAg{1},'LineStyle','-.','LineWidth',2,'MarkerFacecolor','none','color',CAg,'Markersize',10); 
    plot(q,sqrt(RAgX2),MFAg2{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CAg2,'color',CAg2,'Markersize',10); 
    plot(q,sqrt(RAgEX2),MFAg2{1},'LineStyle','-.','LineWidth',2,'MarkerFacecolor','none','color',CAg2,'Markersize',10); hold off; 

    yytick=[0 0.001 0.01 0.025 0.05 0.075 0.1:0.05:0.3];
    set(gca,'LineWidth',2,'Tickdir','out','Fontsize',20,'XTick',[1:14],'Xlim',[1 14],'Xminortick','off','YTick',sqrt(yytick),'ylim',[0 sqrt(yytick(end))],'YTickLabel',num2str(yytick'),'Yminortick','on');
    xlabel('Duration of quarantine (days)','Fontsize',26);
    ylabel({'Probability of','post-quarantine transmission'},'Fontsize',24);
    legend({'RT-PCR on exit',['CareStart (Internal) on exit'],['CareStart (Internal) on entry and exit'],['CareStart (External) on exit'],['CareStart (External) on entry and exit']},'Fontsize',18,'Position',[0.671918776381958,0.78858494220131,0.324579823057817,0.163120562695684]);
    legend boxoff;
    box off;
    grid on;
    
    text(-1.205,sqrt(0.36),'CareStart (Anterior Nasal Swab)','Fontsize',32,'HorizontalAlignment','center') 
    
    text(-1.764,sqrt(0.3),'B','Fontsize',32,'FontWeight','bold');
 elseif(strcmp(testName,'Sofia (CDC)')) 
    [CAg,MFAg]= ColourTests(testName);
    load(['TestingonExit_' testName '_NoDelay_Hellewell.mat'],'IDSLA','IDSLS');
    RAgX=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    load(['TestingonEntryExit_' testName '_NoDelay_Hellewell.mat'],'IDSLA','IDSLS');
    RAgEX=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    fff=strcmp(testNameAll,'Sofia (FDA)');
    [CAg2,MFAg2]= ColourTests(testNameAll{fff});
    CAg2=CAg2(1,:);
    load(['TestingonExit_' testNameAll{fff} '_NoDelay_Hellewell.mat'],'IDSLA','IDSLS');
    RAgX2=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    load(['TestingonEntryExit_' testNameAll{fff} '_NoDelay_Hellewell.mat'],'IDSLA','IDSLS');
    RAgEX2=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    plot(q,sqrt(RRTPCR),MFRTPCR{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CRTPCR,'color',CRTPCR,'Markersize',10); hold on;
    plot(q,sqrt(RAgX2),MFAg2{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CAg2,'color',CAg2,'Markersize',10); 
    plot(q,sqrt(RAgEX2),MFAg2{1},'LineStyle','-.','LineWidth',2,'MarkerFacecolor','none','color',CAg2,'Markersize',10);
    plot(q,sqrt(RAgX),MFAg{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CAg,'color',CAg,'Markersize',10); 
    plot(q,sqrt(RAgEX),MFAg{1},'LineStyle','-.','LineWidth',2,'MarkerFacecolor','none','color',CAg,'Markersize',10);  hold off; 

    yytick=[0 0.001 0.01 0.025 0.05 0.075 0.1:0.05:0.3];
    set(gca,'LineWidth',2,'Tickdir','out','Fontsize',20,'XTick',[1:14],'Xlim',[1 14],'Xminortick','off','YTick',sqrt(yytick),'ylim',[0 sqrt(yytick(end))],'YTickLabel',num2str(yytick'),'Yminortick','on');
    xlabel('Duration of quarantine (days)','Fontsize',26);
    ylabel({'Probability of','post-quarantine transmission'},'Fontsize',24);
    legend({'RT-PCR on exit',['Sofia (Internal) on exit'],['Sofia (Internal) on entry and exit'],['Sofia (External) on exit'],['Sofia (External) on entry and exit']},'Fontsize',18,'Position',[0.671918776381958,0.78858494220131,0.324579823057817,0.163120562695684]);
    legend boxoff;
    box off;
    grid on;
    
    text(-1.205,sqrt(0.36),'Sofia','Fontsize',32,'HorizontalAlignment','center') 
elseif(strcmp(testName,'Liaison (Anterior Nasal Swab)')) 
    [CAg,MFAg]= ColourTests(testName);
    load(['TestingonExit_' testName '_NoDelay_Hellewell.mat'],'IDSLA','IDSLS');
    RAgX=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    load(['TestingonEntryExit_' testName '_NoDelay_Hellewell.mat'],'IDSLA','IDSLS');
    RAgEX=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    fff=strcmp(testNameAll,'Liaison (Nasalpharyngeal Swab)');
    [CAg2,MFAg2]= ColourTests(testNameAll{fff});
    load(['TestingonExit_' testNameAll{fff} '_NoDelay_Hellewell.mat'],'IDSLA','IDSLS');
    RAgX2=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    load(['TestingonEntryExit_' testNameAll{fff} '_NoDelay_Hellewell.mat'],'IDSLA','IDSLS');
    RAgEX2=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    plot(q,sqrt(RRTPCR),MFRTPCR{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CRTPCR,'color',CRTPCR,'Markersize',10); hold on;
    plot(q,sqrt(RAgX),MFAg{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CAg,'color',CAg,'Markersize',10); 
    plot(q,sqrt(RAgEX),MFAg{1},'LineStyle','-.','LineWidth',2,'MarkerFacecolor','none','color',CAg,'Markersize',10); 
    plot(q,sqrt(RAgX2),MFAg2{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CAg2,'color',CAg2,'Markersize',10); 
    plot(q,sqrt(RAgEX2),MFAg2{1},'LineStyle','-.','LineWidth',2,'MarkerFacecolor','none','color',CAg2,'Markersize',10); hold off; 

    yytick=[0 0.001 0.01 0.025 0.05 0.075 0.1:0.05:0.3];
    set(gca,'LineWidth',2,'Tickdir','out','Fontsize',20,'XTick',[1:14],'Xlim',[1 14],'Xminortick','off','YTick',sqrt(yytick),'ylim',[0 sqrt(yytick(end))],'YTickLabel',num2str(yytick'),'Yminortick','on');
    xlabel('Duration of quarantine (days)','Fontsize',26);
    ylabel({'Probability of','post-quarantine transmission'},'Fontsize',24);
    legend({'RT-PCR on exit',[testName ' on exit'],[testName ' on entry and exit'],[testNameAll{fff} ' on exit'],[testNameAll{fff} ' on entry and exit']},'Fontsize',18,'Position',[0.671918776381958,0.78858494220131,0.324579823057817,0.163120562695684]);
    legend boxoff;
    box off;
    grid on;
    text(-1.205,sqrt(0.36),'Liaison','Fontsize',32,'HorizontalAlignment','center') 
    
    text(-1.764,sqrt(0.3),'B','Fontsize',32,'FontWeight','bold');
elseif(strcmp(testName,'LumiraDX (Anterior Nasal Swab)')) 
    [CAg,MFAg]= ColourTests(testName);
    load(['TestingonExit_' testName '_NoDelay_Hellewell.mat'],'IDSLA','IDSLS');
    RAgX=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    load(['TestingonEntryExit_' testName '_NoDelay_Hellewell.mat'],'IDSLA','IDSLS');
    RAgEX=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    fff=strcmp(testNameAll,'LumiraDX (Nasopharyngeal Swabs)');
    [CAg2,MFAg2]= ColourTests(testNameAll{fff});
    load(['TestingonExit_' testNameAll{fff} '_NoDelay_Hellewell.mat'],'IDSLA','IDSLS');
    RAgX2=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    load(['TestingonEntryExit_' testNameAll{fff} '_NoDelay_Hellewell.mat'],'IDSLA','IDSLS');
    RAgEX2=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    plot(q,sqrt(RRTPCR),MFRTPCR{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CRTPCR,'color',CRTPCR,'Markersize',10); hold on;
    plot(q,sqrt(RAgX),MFAg{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CAg,'color',CAg,'Markersize',10); 
    plot(q,sqrt(RAgEX),MFAg{1},'LineStyle','-.','LineWidth',2,'MarkerFacecolor','none','color',CAg,'Markersize',10); 
    plot(q,sqrt(RAgX2),MFAg2{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CAg2,'color',CAg2,'Markersize',10); 
    plot(q,sqrt(RAgEX2),MFAg2{1},'LineStyle','-.','LineWidth',2,'MarkerFacecolor','none','color',CAg2,'Markersize',10); hold off; 

    yytick=[0 0.001 0.01 0.025 0.05 0.075 0.1:0.05:0.3];
    set(gca,'LineWidth',2,'Tickdir','out','Fontsize',20,'XTick',[1:14],'Xlim',[1 14],'Xminortick','off','YTick',sqrt(yytick),'ylim',[0 sqrt(yytick(end))],'YTickLabel',num2str(yytick'),'Yminortick','on');
    xlabel('Duration of quarantine (days)','Fontsize',26);
    ylabel({'Probability of','post-quarantine transmission'},'Fontsize',24);
    legend({'RT-PCR on exit',[testName ' on exit'],[testName ' on entry and exit'],[testNameAll{fff} ' on exit'],[testNameAll{fff} ' on entry and exit']},'Fontsize',18,'Position',[0.671918776381958,0.78858494220131,0.324579823057817,0.163120562695684]);
    legend boxoff;
    box off;
    grid on;
    text(-1.205,sqrt(0.36),'LumiraDX','Fontsize',32,'HorizontalAlignment','center') 
    
    text(-1.764,sqrt(0.3),'B','Fontsize',32,'FontWeight','bold');
elseif(~strcmp(testName,'Sofia (CDC)')&&~strcmp(testName,'CareStart (Anterior Nasal Swab - FDA)')&&~strcmp(testName,'CareStart (Anterior Nasal Swab - External)')&&~strcmp(testName,'LumiraDX (Anterior Nasal Swab)')&&~strcmp(testName,'LumiraDX (Nasopharyngeal Swabs)')&&~strcmp(testName,'Liaison (Anterior Nasal Swab)')&&~strcmp(testName,'Liaison (Nasalpharyngeal Swab)')&&~strcmp(testName,'CareStart (Anterior Nasal Swab)')&&~strcmp(testName,'CareStart (Nasopharyngeal Swab)'))

    [CAg,MFAg]= ColourTests(testName);
    CAg=CAg(end,:);
    load(['TestingonExit_' testName '_NoDelay_Hellewell.mat'],'IDSLA','IDSLS');
    RAgX=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    load(['TestingonEntryExit_' testName '_NoDelay_Hellewell.mat'],'IDSLA','IDSLS');
    RAgEX=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    plot(q,sqrt(RRTPCR),MFRTPCR{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CRTPCR,'color',CRTPCR,'Markersize',10); hold on;
    plot(q,sqrt(RAgX),MFAg{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CAg,'color',CAg,'Markersize',10); 
    plot(q,sqrt(RAgEX),MFAg{1},'LineStyle','-.','LineWidth',2,'MarkerFacecolor','none','color',CAg,'Markersize',10); hold off;

    yytick=[0 0.001 0.01 0.025 0.05 0.075 0.1:0.05:0.3];
    set(gca,'LineWidth',2,'Tickdir','out','Fontsize',20,'XTick',[1:14],'Xlim',[1 14],'Xminortick','off','YTick',sqrt(yytick),'ylim',[0 sqrt(yytick(end))],'YTickLabel',num2str(yytick'),'Yminortick','on');
    xlabel('Duration of quarantine (days)','Fontsize',26);
    ylabel({'Probability of','post-quarantine transmission'},'Fontsize',24);
    legend({'RT-PCR on exit',[testName ' on exit'],[testName ' on entry and exit']},'Fontsize',18);
    legend boxoff;
    box off;
    grid on;
    if(strcmp(testName,'Sofia (FDA)'))
        text(-1.205,sqrt(0.36),'Sofia','Fontsize',32,'HorizontalAlignment','center') 
    else
        text(-1.205,sqrt(0.36),testName,'Fontsize',32,'HorizontalAlignment','center') 
    end
    
    text(-1.764,sqrt(0.3),'B','Fontsize',32,'FontWeight','bold');
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%% Nat Comm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diagnostic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
subplot('Position',[0.052,0.105,0.412,0.39]);

t=linspace(0,40,1001);
if(strcmp(testName,'BinaxNOW')) 
        
    load(['BinaxNOW (FDA)_LR_Parameters.mat'],'beta');
    S = TestSensitivityOLD(t,ts,tL,inf,beta,0);  
    
    
    SR = TestSensitivityOLD(t,ts,tL,inf,[],0);
    [CCtestRTPCR,~]=ColourTests('RTPCR');
    [CCtest,SymP]= ColourTests('BinaxNOW (FDA)'); % As we have this as the baseline
    CCtest=CCtest(2,:);
    p1=plot(t,SR,'-','color',CCtestRTPCR,'LineWidth',2); hold on;
    p2=plot(t,S,'-','color',CCtest,'LineWidth',2); hold on;
    plot(ts.*ones(101,1),linspace(0,1,101),'-.','color',[0.3 0.3 0.3],'LineWidth',2);
    legend([p1,p2],{'RT-PCR','Antigen test'},'Fontsize',18);
    legend boxoff;
    box off;
    set(gca,'LineWidth',2,'tickdir','out','Fontsize',22,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:0.1:1],'Ylim',[0 1]);
    xlabel('Days post-infection','Fontsize',26);
    ylabel({'Diagnostic sensitivty'},'Fontsize',26);
    
    text(-4.745,1,'C','Fontsize',32,'FontWeight','bold');
 elseif(strcmp(testName,'BinaxNOW (Community)')) 
    fff=strcmp(testNameAll,'BinaxNOW (FDA)');
    load([testNameAll{fff} '_LR_Parameters.mat'],'beta');    
    S2 = TestSensitivityOLD(t,ts,tL,inf,beta,0);
    [CCtest2,SymP2]= ColourTests(testNameAll{fff});
    CCtest2=CCtest2(1,:);    
    load([testName '_LR_Parameters.mat'],'beta');
    S = TestSensitivityOLD(t,ts,tL,inf,beta,0);  
    
    
    SR = TestSensitivityOLD(t,ts,tL,inf,[],0);
    [CCtestRTPCR,~]=ColourTests('RTPCR');
    [CCtest,SymP]= ColourTests(testName);
    p1=plot(t,SR,'-','color',CCtestRTPCR,'LineWidth',2); hold on;
    p2=plot(t,S2,'-','color',CCtest2,'LineWidth',2); hold on;
    p3=plot(t,S,'-.','color',CCtest,'LineWidth',2); hold on;
    plot(ts.*ones(101,1),linspace(0,1,101),'-.','color',[0.3 0.3 0.3],'LineWidth',2);
    legend([p1,p2,p3],{'RT-PCR','Antigen test (Internal)','Antigen test (External)'},'Fontsize',18);
    legend boxoff;
    box off;
    set(gca,'LineWidth',2,'tickdir','out','Fontsize',22,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:0.1:1],'Ylim',[0 1]);
    xlabel('Days post-infection','Fontsize',26);
    ylabel({'Diagnostic sensitivty'},'Fontsize',26);
    
    text(-4.745,1,'C','Fontsize',32,'FontWeight','bold');
elseif(strcmp(testName,'CareStart (Anterior Nasal Swab)')) 
    fff=strcmp(testNameAll,'CareStart (Nasopharyngeal Swab)');
    load([testNameAll{fff} '_LR_Parameters.mat'],'beta');    
    S2 = TestSensitivityOLD(t,ts,tL,inf,beta,0);
    [CCtest2,SymP2]= ColourTests(testNameAll{fff});
        
    load([testName '_LR_Parameters.mat'],'beta');
    S = TestSensitivityOLD(t,ts,tL,inf,beta,0);  
    
    
    SR = TestSensitivityOLD(t,ts,tL,inf,[],0);
    [CCtestRTPCR,~]=ColourTests('RTPCR');
    [CCtest,SymP]= ColourTests('CareStart (Anterior Nasal Swab - FDA)'); % As we have this as the baseline
    CCtest=CCtest(2,:);
    p1=plot(t,SR,'-','color',CCtestRTPCR,'LineWidth',2); hold on;
    p2=plot(t,S,'-','color',CCtest,'LineWidth',2); hold on;
    p3=plot(t,S2,'-.','color',CCtest2,'LineWidth',2); hold on;
    plot(ts.*ones(101,1),linspace(0,1,101),'-.','color',[0.3 0.3 0.3],'LineWidth',2);
    legend([p1,p2,p3],{'RT-PCR','Antigen test (Anterior)','Antigen test (Nasopharyngeal)'},'Fontsize',18);
    legend boxoff;
    box off;
    set(gca,'LineWidth',2,'tickdir','out','Fontsize',22,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:0.1:1],'Ylim',[0 1]);
    xlabel('Days post-infection','Fontsize',26);
    ylabel({'Diagnostic sensitivty'},'Fontsize',26);
    
    text(-4.745,1,'C','Fontsize',32,'FontWeight','bold');
elseif(strcmp(testName,'CareStart (Anterior Nasal Swab - FDA)')) 
    fff=strcmp(testNameAll,'CareStart (Anterior Nasal Swab - External)');
    load([testNameAll{fff} '_LR_Parameters.mat'],'beta');    
    S2 = TestSensitivityOLD(t,ts,tL,inf,beta,0);
    [CCtest2,SymP2]= ColourTests(testNameAll{fff});
        
    load([testName '_LR_Parameters.mat'],'beta');
    S = TestSensitivityOLD(t,ts,tL,inf,beta,0);  
    
    
    SR = TestSensitivityOLD(t,ts,tL,inf,[],0);
    [CCtestRTPCR,~]=ColourTests('RTPCR');
    [CCtest,SymP]= ColourTests(testName);
    CCtest=CCtest(1,:);
    p1=plot(t,SR,'-','color',CCtestRTPCR,'LineWidth',2); hold on;
    p2=plot(t,S,'-','color',CCtest,'LineWidth',2); hold on;
    p3=plot(t,S2,'-.','color',CCtest2,'LineWidth',2); hold on;
    plot(ts.*ones(101,1),linspace(0,1,101),'-.','color',[0.3 0.3 0.3],'LineWidth',2);
    legend([p1,p2,p3],{'RT-PCR','Antigen test (Internal)','Antigen test (External)'},'Fontsize',18);
    legend boxoff;
    box off;
    set(gca,'LineWidth',2,'tickdir','out','Fontsize',22,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:0.1:1],'Ylim',[0 1]);
    xlabel('Days post-infection','Fontsize',26);
    ylabel({'Diagnostic sensitivty'},'Fontsize',26);
    
    text(-4.745,1,'C','Fontsize',32,'FontWeight','bold');
 elseif(strcmp(testName,'Sofia (CDC)')) 
    fff=strcmp(testNameAll,'Sofia (FDA)');
    load([testNameAll{fff} '_LR_Parameters.mat'],'beta');    
    S2 = TestSensitivityOLD(t,ts,tL,inf,beta,0);
    [CCtest2,SymP2]= ColourTests(testNameAll{fff});
    CCtest2=CCtest2(1,:);    
    load([testName '_LR_Parameters.mat'],'beta');
    S = TestSensitivityOLD(t,ts,tL,inf,beta,0);  
    
    
    SR = TestSensitivityOLD(t,ts,tL,inf,[],0);
    [CCtestRTPCR,~]=ColourTests('RTPCR');
    [CCtest,SymP]= ColourTests(testName);
    p1=plot(t,SR,'-','color',CCtestRTPCR,'LineWidth',2); hold on;
    p2=plot(t,S2,'-','color',CCtest2,'LineWidth',2); hold on;
    p3=plot(t,S,'-.','color',CCtest,'LineWidth',2); hold on;
    plot(ts.*ones(101,1),linspace(0,1,101),'-.','color',[0.3 0.3 0.3],'LineWidth',2);
    legend([p1,p2,p3],{'RT-PCR','Antigen test (Internal)','Antigen test (External)'},'Fontsize',18);
    legend boxoff;
    box off;
    set(gca,'LineWidth',2,'tickdir','out','Fontsize',22,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:0.1:1],'Ylim',[0 1]);
    xlabel('Days post-infection','Fontsize',26);
    ylabel({'Diagnostic sensitivty'},'Fontsize',26);
    
    text(-4.745,1,'C','Fontsize',32,'FontWeight','bold');
elseif(strcmp(testName,'Liaison (Anterior Nasal Swab)')) 
    fff=strcmp(testNameAll,'Liaison (Nasalpharyngeal Swab)');
    load([testNameAll{fff} '_LR_Parameters.mat'],'beta');    
    S2 = TestSensitivityOLD(t,ts,tL,inf,beta,0);
    [CCtest2,SymP2]= ColourTests(testNameAll{fff});
        
    load([testName '_LR_Parameters.mat'],'beta');
    S = TestSensitivityOLD(t,ts,tL,inf,beta,0);  
    
    
    SR = TestSensitivityOLD(t,ts,tL,inf,[],0);
    [CCtestRTPCR,~]=ColourTests('RTPCR');
    [CCtest,SymP]= ColourTests(testName);
    p1=plot(t,SR,'-','color',CCtestRTPCR,'LineWidth',2); hold on;
    p2=plot(t,S,'-','color',CCtest,'LineWidth',2); hold on;
    p3=plot(t,S2,'-.','color',CCtest2,'LineWidth',2); hold on;
    plot(ts.*ones(101,1),linspace(0,1,101),'-.','color',[0.3 0.3 0.3],'LineWidth',2);
    legend([p1,p2,p3],{'RT-PCR','Antigen test (Anterior)','Antigen test (Nasopharyngeal)'},'Fontsize',18);
    legend boxoff;
    box off;
    set(gca,'LineWidth',2,'tickdir','out','Fontsize',22,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:0.1:1],'Ylim',[0 1]);
    xlabel('Days post-infection','Fontsize',26);
    ylabel({'Diagnostic sensitivty'},'Fontsize',26);
    
    text(-4.745,1,'C','Fontsize',32,'FontWeight','bold');
elseif(strcmp(testName,'LumiraDX (Anterior Nasal Swab)')) 
    fff=strcmp(testNameAll,'LumiraDX (Nasopharyngeal Swabs)');
    load([testNameAll{fff} '_LR_Parameters.mat'],'beta');    
    S2 = TestSensitivityOLD(t,ts,tL,inf,beta,0);
    [CCtest2,SymP2]= ColourTests(testNameAll{fff});
        
    load([testName '_LR_Parameters.mat'],'beta');
    S = TestSensitivityOLD(t,ts,tL,inf,beta,0);  
    
    
    SR = TestSensitivityOLD(t,ts,tL,inf,[],0);
    [CCtestRTPCR,~]=ColourTests('RTPCR');
    [CCtest,SymP]= ColourTests(testName);
    p1=plot(t,SR,'-','color',CCtestRTPCR,'LineWidth',2); hold on;
    p2=plot(t,S,'-','color',CCtest,'LineWidth',2); hold on;
    p3=plot(t,S2,'-.','color',CCtest2,'LineWidth',2); hold on;
    plot(ts.*ones(101,1),linspace(0,1,101),'-.','color',[0.3 0.3 0.3],'LineWidth',2);
    legend([p1,p2,p3],{'RT-PCR','Antigen test (Anterior)','Antigen test (Nasopharyngeal)'},'Fontsize',18);
    legend boxoff;
    box off;
    set(gca,'LineWidth',2,'tickdir','out','Fontsize',22,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:0.1:1],'Ylim',[0 1]);
    xlabel('Days post-infection','Fontsize',26);
    ylabel({'Diagnostic sensitivty'},'Fontsize',26);
    
    text(-4.745,1,'C','Fontsize',32,'FontWeight','bold');
elseif(~strcmp(testName,'Sofia (CDC)')&&~strcmp(testName,'CareStart (Anterior Nasal Swab - FDA)')&&~strcmp(testName,'CareStart (Anterior Nasal Swab - External)')&&~strcmp(testName,'LumiraDX (Anterior Nasal Swab)')&&~strcmp(testName,'LumiraDX (Nasopharyngeal Swabs)')&&~strcmp(testName,'Liaison (Anterior Nasal Swab)')&&~strcmp(testName,'Liaison (Nasalpharyngeal Swab)')&&~strcmp(testName,'CareStart (Anterior Nasal Swab)')&&~strcmp(testName,'CareStart (Nasopharyngeal Swab)'))
    load([testName '_LR_Parameters.mat'],'beta');

    S = TestSensitivityOLD(t,ts,tL,inf,beta,0);
    SR = TestSensitivityOLD(t,ts,tL,inf,[],0);
    [CCtestRTPCR,~]=ColourTests('RTPCR');
    [CCtest,SymP]= ColourTests(testName);
    CCtest=CCtest(end,:);
    p1=plot(t,SR,'-','color',CCtestRTPCR,'LineWidth',2); hold on;
    p2=plot(t,S,'-','color',CCtest,'LineWidth',2); hold on;
    plot(ts.*ones(101,1),linspace(0,1,101),'-.','color',[0.3 0.3 0.3],'LineWidth',2);
    legend([p1,p2],{'RT-PCR','Antigen test'},'Fontsize',18);
    legend boxoff;
    box off;
    set(gca,'LineWidth',2,'tickdir','out','Fontsize',22,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:0.1:1],'Ylim',[0 1]);
    xlabel('Days post-infection','Fontsize',26);
    ylabel({'Diagnostic sensitivty'},'Fontsize',26);
    
    text(-4.745,1,'C','Fontsize',32,'FontWeight','bold');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Probabilility of PQT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot('Position',[0.5745,0.105,0.412,0.39]);
Risk=1;
load('TestingonExit_RTPCR_24hrDelay_NatComm.mat','IDSLA','IDSLS','pA','q')
RRTPCR=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);
[CRTPCR,MFRTPCR]= ColourTests('RTPCR');
if(strcmp(testName,'BinaxNOW')) 
    [CAg,MFAg]= ColourTests('BinaxNOW (FDA)'); % as we have this as the baseline
    CAg=CAg(2,:);
    load(['TestingonExit_BinaxNOW (FDA)_NoDelay_NatComm.mat'],'IDSLA','IDSLS');
    RAgX=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    load(['TestingonEntryExit_BinaxNOW (FDA)_NoDelay_NatComm.mat'],'IDSLA','IDSLS');
    RAgEX=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    plot(q,sqrt(RRTPCR),MFRTPCR{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CRTPCR,'color',CRTPCR,'Markersize',10); hold on;
    plot(q,sqrt(RAgX),MFAg{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CAg,'color',CAg,'Markersize',10); 
    plot(q,sqrt(RAgEX),MFAg{1},'LineStyle','-.','LineWidth',2,'MarkerFacecolor','none','color',CAg,'Markersize',10);  hold off; 

    yytick=[0 0.001 0.01 0.025 0.05 0.075 0.1:0.05:0.3];
    set(gca,'LineWidth',2,'Tickdir','out','Fontsize',20,'XTick',[1:14],'Xlim',[1 14],'Xminortick','off','YTick',sqrt(yytick),'ylim',[0 sqrt(yytick(end))],'YTickLabel',num2str(yytick'),'Yminortick','on');
    xlabel('Duration of quarantine (days)','Fontsize',26);
    ylabel({'Probability of','post-quarantine transmission'},'Fontsize',24);
    text(-1.764,sqrt(0.3),'D','Fontsize',32,'FontWeight','bold');
%     legend({'RT-PCR on exit',[testName ' on exit'],[testName ' on entry and exit'],[testNameAll{fff} ' on exit'],[testNameAll{fff} ' on entry and exit']},'Fontsize',18,'Position',[0.671918776381958,0.78858494220131,0.324579823057817,0.163120562695684]);
%     legend boxoff;
    box off;
    grid on;
print(gcf,['SIFigure_' testName '.png'],'-dpng','-r300');
elseif(strcmp(testName,'BinaxNOW (Community)')) 
    [CAg,MFAg]= ColourTests(testName);
    load(['TestingonExit_' testName '_NoDelay_NatComm.mat'],'IDSLA','IDSLS');
    RAgX=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    load(['TestingonEntryExit_' testName '_NoDelay_NatComm.mat'],'IDSLA','IDSLS');
    RAgEX=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    fff=strcmp(testNameAll,'BinaxNOW (FDA)');
    [CAg2,MFAg2]= ColourTests(testNameAll{fff});
    CAg2=CAg2(1,:);
    load(['TestingonExit_' testNameAll{fff} '_NoDelay_NatComm.mat'],'IDSLA','IDSLS');
    RAgX2=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    load(['TestingonEntryExit_' testNameAll{fff} '_NoDelay_NatComm.mat'],'IDSLA','IDSLS');
    RAgEX2=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    plot(q,sqrt(RRTPCR),MFRTPCR{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CRTPCR,'color',CRTPCR,'Markersize',10); hold on;
    plot(q,sqrt(RAgX2),MFAg2{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CAg2,'color',CAg2,'Markersize',10); 
    plot(q,sqrt(RAgEX2),MFAg2{1},'LineStyle','-.','LineWidth',2,'MarkerFacecolor','none','color',CAg2,'Markersize',10);
    plot(q,sqrt(RAgX),MFAg{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CAg,'color',CAg,'Markersize',10); 
    plot(q,sqrt(RAgEX),MFAg{1},'LineStyle','-.','LineWidth',2,'MarkerFacecolor','none','color',CAg,'Markersize',10);  hold off; 

    yytick=[0 0.001 0.01 0.025 0.05 0.075 0.1:0.05:0.3];
    set(gca,'LineWidth',2,'Tickdir','out','Fontsize',20,'XTick',[1:14],'Xlim',[1 14],'Xminortick','off','YTick',sqrt(yytick),'ylim',[0 sqrt(yytick(end))],'YTickLabel',num2str(yytick'),'Yminortick','on');
    xlabel('Duration of quarantine (days)','Fontsize',26);
    ylabel({'Probability of','post-quarantine transmission'},'Fontsize',24);
    text(-1.764,sqrt(0.3),'D','Fontsize',32,'FontWeight','bold');
%     legend({'RT-PCR on exit',[testName ' on exit'],[testName ' on entry and exit'],[testNameAll{fff} ' on exit'],[testNameAll{fff} ' on entry and exit']},'Fontsize',18,'Position',[0.671918776381958,0.78858494220131,0.324579823057817,0.163120562695684]);
%     legend boxoff;
    box off;
    grid on;
print(gcf,['SIFigure_' testName '.png'],'-dpng','-r300');
elseif(strcmp(testName,'CareStart (Anterior Nasal Swab)')) 
    [CAg,MFAg]= ColourTests('CareStart (Anterior Nasal Swab - FDA)'); % as we have this as the baseline
    CAg=CAg(2,:);
    load(['TestingonExit_CareStart (Anterior Nasal Swab - FDA)_NoDelay_NatComm.mat'],'IDSLA','IDSLS');
    RAgX=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    load(['TestingonEntryExit_CareStart (Anterior Nasal Swab - FDA)_NoDelay_NatComm.mat'],'IDSLA','IDSLS');
    RAgEX=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    fff=strcmp(testNameAll,'CareStart (Nasopharyngeal Swab)');
    [CAg2,MFAg2]= ColourTests(testNameAll{fff});
    load(['TestingonExit_' testNameAll{fff} '_NoDelay_NatComm.mat'],'IDSLA','IDSLS');
    RAgX2=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    load(['TestingonEntryExit_' testNameAll{fff} '_NoDelay_NatComm.mat'],'IDSLA','IDSLS');
    RAgEX2=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    plot(q,sqrt(RRTPCR),MFRTPCR{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CRTPCR,'color',CRTPCR,'Markersize',10); hold on;
    plot(q,sqrt(RAgX),MFAg{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CAg,'color',CAg,'Markersize',10); 
    plot(q,sqrt(RAgEX),MFAg{1},'LineStyle','-.','LineWidth',2,'MarkerFacecolor','none','color',CAg,'Markersize',10); 
    plot(q,sqrt(RAgX2),MFAg2{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CAg2,'color',CAg2,'Markersize',10); 
    plot(q,sqrt(RAgEX2),MFAg2{1},'LineStyle','-.','LineWidth',2,'MarkerFacecolor','none','color',CAg2,'Markersize',10); hold off; 

    yytick=[0 0.001 0.01 0.025 0.05 0.075 0.1:0.05:0.3];
    set(gca,'LineWidth',2,'Tickdir','out','Fontsize',20,'XTick',[1:14],'Xlim',[1 14],'Xminortick','off','YTick',sqrt(yytick),'ylim',[0 sqrt(yytick(end))],'YTickLabel',num2str(yytick'),'Yminortick','on');
    xlabel('Duration of quarantine (days)','Fontsize',26);
    ylabel({'Probability of','post-quarantine transmission'},'Fontsize',24);
    text(-1.764,sqrt(0.3),'D','Fontsize',32,'FontWeight','bold');
%     legend({'RT-PCR on exit',[testName ' on exit'],[testName ' on entry and exit'],[testNameAll{fff} ' on exit'],[testNameAll{fff} ' on entry and exit']},'Fontsize',18,'Position',[0.671918776381958,0.78858494220131,0.324579823057817,0.163120562695684]);
%     legend boxoff;
    box off;
    grid on;
print(gcf,['SIFigure_' testName '.png'],'-dpng','-r300');
elseif(strcmp(testName,'CareStart (Anterior Nasal Swab - FDA)')) 
    [CAg,MFAg]= ColourTests(testName);
    CAg=CAg(1,:);
    load(['TestingonExit_' testName '_NoDelay_NatComm.mat'],'IDSLA','IDSLS');
    RAgX=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    load(['TestingonEntryExit_' testName '_NoDelay_NatComm.mat'],'IDSLA','IDSLS');
    RAgEX=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    fff=strcmp(testNameAll,'CareStart (Anterior Nasal Swab - External)');
    [CAg2,MFAg2]= ColourTests(testNameAll{fff});
    load(['TestingonExit_' testNameAll{fff} '_NoDelay_NatComm.mat'],'IDSLA','IDSLS');
    RAgX2=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    load(['TestingonEntryExit_' testNameAll{fff} '_NoDelay_NatComm.mat'],'IDSLA','IDSLS');
    RAgEX2=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    plot(q,sqrt(RRTPCR),MFRTPCR{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CRTPCR,'color',CRTPCR,'Markersize',10); hold on;
    plot(q,sqrt(RAgX),MFAg{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CAg,'color',CAg,'Markersize',10); 
    plot(q,sqrt(RAgEX),MFAg{1},'LineStyle','-.','LineWidth',2,'MarkerFacecolor','none','color',CAg,'Markersize',10); 
    plot(q,sqrt(RAgX2),MFAg2{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CAg2,'color',CAg2,'Markersize',10); 
    plot(q,sqrt(RAgEX2),MFAg2{1},'LineStyle','-.','LineWidth',2,'MarkerFacecolor','none','color',CAg2,'Markersize',10); hold off; 

    yytick=[0 0.001 0.01 0.025 0.05 0.075 0.1:0.05:0.3];
    set(gca,'LineWidth',2,'Tickdir','out','Fontsize',20,'XTick',[1:14],'Xlim',[1 14],'Xminortick','off','YTick',sqrt(yytick),'ylim',[0 sqrt(yytick(end))],'YTickLabel',num2str(yytick'),'Yminortick','on');
    xlabel('Duration of quarantine (days)','Fontsize',26);
    ylabel({'Probability of','post-quarantine transmission'},'Fontsize',24);
    text(-1.764,sqrt(0.3),'D','Fontsize',32,'FontWeight','bold');
%     legend({'RT-PCR on exit',[testName ' on exit'],[testName ' on entry and exit'],[testNameAll{fff} ' on exit'],[testNameAll{fff} ' on entry and exit']},'Fontsize',18,'Position',[0.671918776381958,0.78858494220131,0.324579823057817,0.163120562695684]);
%     legend boxoff;
    box off;
    grid on;
print(gcf,['SIFigure_' testName '.png'],'-dpng','-r300');
elseif(strcmp(testName,'Sofia (CDC)')) 
    [CAg,MFAg]= ColourTests(testName);
    load(['TestingonExit_' testName '_NoDelay_NatComm.mat'],'IDSLA','IDSLS');
    RAgX=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    load(['TestingonEntryExit_' testName '_NoDelay_NatComm.mat'],'IDSLA','IDSLS');
    RAgEX=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    fff=strcmp(testNameAll,'Sofia (FDA)');
    [CAg2,MFAg2]= ColourTests(testNameAll{fff});
    CAg2=CAg2(1,:);
    load(['TestingonExit_' testNameAll{fff} '_NoDelay_NatComm.mat'],'IDSLA','IDSLS');
    RAgX2=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    load(['TestingonEntryExit_' testNameAll{fff} '_NoDelay_NatComm.mat'],'IDSLA','IDSLS');
    RAgEX2=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    plot(q,sqrt(RRTPCR),MFRTPCR{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CRTPCR,'color',CRTPCR,'Markersize',10); hold on;
    plot(q,sqrt(RAgX2),MFAg2{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CAg2,'color',CAg2,'Markersize',10); 
    plot(q,sqrt(RAgEX2),MFAg2{1},'LineStyle','-.','LineWidth',2,'MarkerFacecolor','none','color',CAg2,'Markersize',10);
    plot(q,sqrt(RAgX),MFAg{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CAg,'color',CAg,'Markersize',10); 
    plot(q,sqrt(RAgEX),MFAg{1},'LineStyle','-.','LineWidth',2,'MarkerFacecolor','none','color',CAg,'Markersize',10);  hold off; 

    yytick=[0 0.001 0.01 0.025 0.05 0.075 0.1:0.05:0.3];
    set(gca,'LineWidth',2,'Tickdir','out','Fontsize',20,'XTick',[1:14],'Xlim',[1 14],'Xminortick','off','YTick',sqrt(yytick),'ylim',[0 sqrt(yytick(end))],'YTickLabel',num2str(yytick'),'Yminortick','on');
    xlabel('Duration of quarantine (days)','Fontsize',26);
    ylabel({'Probability of','post-quarantine transmission'},'Fontsize',24);
    text(-1.764,sqrt(0.3),'D','Fontsize',32,'FontWeight','bold');
%     legend({'RT-PCR on exit',[testName ' on exit'],[testName ' on entry and exit'],[testNameAll{fff} ' on exit'],[testNameAll{fff} ' on entry and exit']},'Fontsize',18,'Position',[0.671918776381958,0.78858494220131,0.324579823057817,0.163120562695684]);
%     legend boxoff;
    box off;
    grid on;
print(gcf,['SIFigure_' testName '.png'],'-dpng','-r300');
elseif(strcmp(testName,'Liaison (Anterior Nasal Swab)')) 
    [CAg,MFAg]= ColourTests(testName);
    load(['TestingonExit_' testName '_NoDelay_NatComm.mat'],'IDSLA','IDSLS');
    RAgX=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    load(['TestingonEntryExit_' testName '_NoDelay_NatComm.mat'],'IDSLA','IDSLS');
    RAgEX=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    fff=strcmp(testNameAll,'Liaison (Nasalpharyngeal Swab)');
    [CAg2,MFAg2]= ColourTests(testNameAll{fff});
    load(['TestingonExit_' testNameAll{fff} '_NoDelay_NatComm.mat'],'IDSLA','IDSLS');
    RAgX2=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    load(['TestingonEntryExit_' testNameAll{fff} '_NoDelay_NatComm.mat'],'IDSLA','IDSLS');
    RAgEX2=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    plot(q,sqrt(RRTPCR),MFRTPCR{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CRTPCR,'color',CRTPCR,'Markersize',10); hold on;
    plot(q,sqrt(RAgX),MFAg{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CAg,'color',CAg,'Markersize',10); 
    plot(q,sqrt(RAgEX),MFAg{1},'LineStyle','-.','LineWidth',2,'MarkerFacecolor','none','color',CAg,'Markersize',10); 
    plot(q,sqrt(RAgX2),MFAg2{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CAg2,'color',CAg2,'Markersize',10); 
    plot(q,sqrt(RAgEX2),MFAg2{1},'LineStyle','-.','LineWidth',2,'MarkerFacecolor','none','color',CAg2,'Markersize',10); hold off; 

    yytick=[0 0.001 0.01 0.025 0.05 0.075 0.1:0.05:0.3];
    set(gca,'LineWidth',2,'Tickdir','out','Fontsize',20,'XTick',[1:14],'Xlim',[1 14],'Xminortick','off','YTick',sqrt(yytick),'ylim',[0 sqrt(yytick(end))],'YTickLabel',num2str(yytick'),'Yminortick','on');
    xlabel('Duration of quarantine (days)','Fontsize',26);
    ylabel({'Probability of','post-quarantine transmission'},'Fontsize',24);
    text(-1.764,sqrt(0.3),'D','Fontsize',32,'FontWeight','bold');
%     legend({'RT-PCR on exit',[testName ' on exit'],[testName ' on entry and exit'],[testNameAll{fff} ' on exit'],[testNameAll{fff} ' on entry and exit']},'Fontsize',18,'Position',[0.671918776381958,0.78858494220131,0.324579823057817,0.163120562695684]);
%     legend boxoff;
    box off;
    grid on;
print(gcf,['SIFigure_' testName '.png'],'-dpng','-r300');
elseif(strcmp(testName,'LumiraDX (Anterior Nasal Swab)')) 
    [CAg,MFAg]= ColourTests(testName);
    load(['TestingonExit_' testName '_NoDelay_NatComm.mat'],'IDSLA','IDSLS');
    RAgX=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    load(['TestingonEntryExit_' testName '_NoDelay_NatComm.mat'],'IDSLA','IDSLS');
    RAgEX=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    fff=strcmp(testNameAll,'LumiraDX (Nasopharyngeal Swabs)');
    [CAg2,MFAg2]= ColourTests(testNameAll{fff});
    load(['TestingonExit_' testNameAll{fff} '_NoDelay_NatComm.mat'],'IDSLA','IDSLS');
    RAgX2=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    load(['TestingonEntryExit_' testNameAll{fff} '_NoDelay_NatComm.mat'],'IDSLA','IDSLS');
    RAgEX2=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    plot(q,sqrt(RRTPCR),MFRTPCR{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CRTPCR,'color',CRTPCR,'Markersize',10); hold on;
    plot(q,sqrt(RAgX),MFAg{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CAg,'color',CAg,'Markersize',10); 
    plot(q,sqrt(RAgEX),MFAg{1},'LineStyle','-.','LineWidth',2,'MarkerFacecolor','none','color',CAg,'Markersize',10); 
    plot(q,sqrt(RAgX2),MFAg2{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CAg2,'color',CAg2,'Markersize',10); 
    plot(q,sqrt(RAgEX2),MFAg2{1},'LineStyle','-.','LineWidth',2,'MarkerFacecolor','none','color',CAg2,'Markersize',10); hold off; 

    yytick=[0 0.001 0.01 0.025 0.05 0.075 0.1:0.05:0.3];
    set(gca,'LineWidth',2,'Tickdir','out','Fontsize',20,'XTick',[1:14],'Xlim',[1 14],'Xminortick','off','YTick',sqrt(yytick),'ylim',[0 sqrt(yytick(end))],'YTickLabel',num2str(yytick'),'Yminortick','on');
    xlabel('Duration of quarantine (days)','Fontsize',26);
    ylabel({'Probability of','post-quarantine transmission'},'Fontsize',24);
    text(-1.764,sqrt(0.3),'D','Fontsize',32,'FontWeight','bold');
%     legend({'RT-PCR on exit',[testName ' on exit'],[testName ' on entry and exit'],[testNameAll{fff} ' on exit'],[testNameAll{fff} ' on entry and exit']},'Fontsize',18,'Position',[0.671918776381958,0.78858494220131,0.324579823057817,0.163120562695684]);
%     legend boxoff;
    box off;
    grid on;
print(gcf,['SIFigure_' testName '.png'],'-dpng','-r300');
elseif(~strcmp(testName,'Sofia (CDC)')&&~strcmp(testName,'CareStart (Anterior Nasal Swab - FDA)')&&~strcmp(testName,'CareStart (Anterior Nasal Swab - External)')&&~strcmp(testName,'LumiraDX (Anterior Nasal Swab)')&&~strcmp(testName,'LumiraDX (Nasopharyngeal Swabs)')&&~strcmp(testName,'Liaison (Anterior Nasal Swab)')&&~strcmp(testName,'Liaison (Nasalpharyngeal Swab)')&&~strcmp(testName,'CareStart (Anterior Nasal Swab)')&&~strcmp(testName,'CareStart (Nasopharyngeal Swab)'))   
    [CAg,MFAg]= ColourTests(testName);
    CAg=CAg(end,:);
    load(['TestingonExit_' testName '_NoDelay_NatComm.mat'],'IDSLA','IDSLS');
    RAgX=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    load(['TestingonEntryExit_' testName '_NoDelay_NatComm.mat'],'IDSLA','IDSLS');
    RAgEX=Probability_Onward((1-pA).*IDSLS+pA.*IDSLA,Risk);

    plot(q,sqrt(RRTPCR),MFRTPCR{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CRTPCR,'color',CRTPCR,'Markersize',10); hold on;
    plot(q,sqrt(RAgX),MFAg{1},'LineStyle','-','LineWidth',2,'MarkerFacecolor',CAg,'color',CAg,'Markersize',10); 
    plot(q,sqrt(RAgEX),MFAg{1},'LineStyle','-.','LineWidth',2,'MarkerFacecolor','none','color',CAg,'Markersize',10); hold off;

    yytick=[0 0.001 0.01 0.025 0.05 0.075 0.1:0.05:0.3];
    set(gca,'LineWidth',2,'Tickdir','out','Fontsize',20,'XTick',[1:14],'Xlim',[1 14],'Xminortick','off','YTick',sqrt(yytick),'ylim',[0 sqrt(yytick(end))],'YTickLabel',num2str(yytick'),'Yminortick','on');
    xlabel('Duration of quarantine (days)','Fontsize',26);
    ylabel({'Probability of','post-quarantine transmission'},'Fontsize',24);
    
    
    text(-1.764,sqrt(0.3),'D','Fontsize',32,'FontWeight','bold');
%     legend({'RT-PCR on exit',[testName ' on exit'],[testName ' on entry and exit']},'Fontsize',18);
%     legend boxoff;
    box off;
    grid on;
print(gcf,['SIFigure_' testName '.png'],'-dpng','-r300');
end
end