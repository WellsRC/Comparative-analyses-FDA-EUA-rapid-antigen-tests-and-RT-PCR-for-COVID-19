clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Sensitivity for a RA test with and without the cutoff relative to RT-PCR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
addpath([pwd '\Delta_Variant']);
addpath([pwd '\Delta_Variant\Results']);
[pA,~,~,ts,~] = BaselineParameters;
CTest=[hex2rgb('#231B12');hex2rgb('#375E97');hex2rgb('#486824');hex2rgb('#F9A603');hex2rgb('#CF3721');hex2rgb('#810f7c')];

figure('units','normalized','outerposition',[0.2 0.2 0.6 0.6]);
subplot('Position',[0.085387323943662,0.174774774774775,0.899295774647887,0.803603603603605]);
tt=1;

t=linspace(0,40,1001);

 
[betaRTPCR,betaAg]=ParameterCOVIDTest('BD Veritor',1);
 S = TestSensitivity(t,ts,[],betaRTPCR); 
 
p1=plot(t,S,'-','color',CTest(1,:),'LineWidth',2);
 hold on
 
Sf = TestSensitivity(t,ts,betaAg,betaRTPCR);


rmpath([pwd '\Delta_Variant']);
rmpath([pwd '\Delta_Variant\Results']);

addpath([pwd '\Ag_Cut']);

    [~,~,~,ts,td] = BaselineParameters;
AgCutoffPSO=fminbnd(@(x)(integral(@(t)ViralShedding_Symptomatic(t,td),0,ts+x)-0.99).^2,0,10);
S = TestSensitivity(t,ts,betaAg,betaRTPCR,AgCutoffPSO);

p4=plot(t,S,'-','color',hex2rgb('D9B44A'),'LineWidth',2);

AgCutoffPSO=10;
S = TestSensitivity(t,ts,betaAg,betaRTPCR,AgCutoffPSO);


p3=plot(t,S,'-','color',hex2rgb('8D230F'),'LineWidth',2);



p2=plot(t,Sf,':','color',CTest(5,:),'LineWidth',1.5);
rmpath([pwd '\Ag_Cut']);
plot(ts.*ones(101,1),linspace(0,1,101),'-.','color',[0.7 0.7 0.7],'LineWidth',2);
set(gca,'LineWidth',2,'tickdir','out','Fontsize',20,'XTick',[0:2:40],'xlim',[0 22],'XMinorTick','on','Yminortick','on','YTick',[0:0.2:1],'Ylim',[0 1]);
xlabel('Days post-infection','Fontsize',24);
ylabel('Diagnostic sensitivity','Fontsize',24);
box off;
legend([p1 p2 p3 p4],{'RT-PCR','BD Veritor','BD Veritor (10 day post-symptom onset based cut-off)','BD Veritor (5.6 day post-symptom onset based cut-off)'},'Fontsize',18);
legend boxoff;
xlim([0 18]);
print(gcf,['Figure_Sensitivity_Ag_Cut.png'],'-dpng','-r600');