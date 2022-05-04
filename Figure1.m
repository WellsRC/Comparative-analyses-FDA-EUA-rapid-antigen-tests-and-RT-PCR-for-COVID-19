clear;
% close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Hellewell et al
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
addpath([pwd '\Delta_Variant']);
addpath([pwd '\Delta_Variant\Results']);
[pA,~,~,ts,~] = BaselineParameters;
CTest=[hex2rgb('#231B12');hex2rgb('#bfd3e6');hex2rgb('#41ae76');hex2rgb('#F9A603');hex2rgb('#993404');hex2rgb('#88419d')];

figure('units','normalized','outerposition',[0.2 0.2 0.6 0.6]);
subplot('Position',[0.085387323943662,0.174774774774775,0.899295774647887,0.803603603603605]);
tt=1;

t=linspace(0,40,1001);

 [betaRTPCR,betaAg]=ParameterCOVIDTest('LumiraDX (Anterior Nasal Swab)',1);
 S = TestSensitivity(t,ts,[],betaRTPCR); 
 
plot(t,S,'-p','color',CTest(1,:),'LineWidth',2,'MarkerFaceColor',CTest(1,:),'MarkerIndices',1:25:length(t));
 hold on
 
 
 S = TestSensitivity(t,ts,betaAg,betaRTPCR);
 plot(t,S,'-s','color',CTest(2,:),'LineWidth',2,'MarkerFaceColor',CTest(2,:),'MarkerIndices',1:25:length(t));
 
[betaRTPCR,betaAg]=ParameterCOVIDTest('Sofia (FDA)',1);
S = TestSensitivity(t,ts,betaAg,betaRTPCR);
plot(t,S,'-d','color',CTest(3,:),'LineWidth',2,'MarkerFaceColor',CTest(3,:),'MarkerIndices',1:25:length(t));

[betaRTPCR,betaAg]=ParameterCOVIDTest('BinaxNOW (FDA)',1);
S = TestSensitivity(t,ts,betaAg,betaRTPCR);

plot(t,S,'-^','color',CTest(4,:),'LineWidth',2,'MarkerFaceColor',CTest(4,:),'MarkerIndices',1:25:length(t));

[betaRTPCR,betaAg]=ParameterCOVIDTest('BD Veritor',1);
S = TestSensitivity(t,ts,betaAg,betaRTPCR);

plot(t,S,'-o','color',CTest(5,:),'LineWidth',2,'MarkerFaceColor',CTest(5,:),'MarkerIndices',1:25:length(t));

[betaRTPCR,betaAg]=ParameterCOVIDTest('CareStart (Anterior Nasal Swab - FDA)',1);
S = TestSensitivity(t,ts,betaAg,betaRTPCR);

plot(t,S,'-h','color',CTest(6,:),'LineWidth',2,'MarkerFaceColor',CTest(5,:),'MarkerIndices',1:25:length(t));

plot(ts.*ones(101,1),linspace(0,1,101),'-.','color',[0.7 0.7 0.7],'LineWidth',2);
set(gca,'LineWidth',2,'tickdir','out','Fontsize',20,'XTick',[0:2:40],'xlim',[0 22],'XMinorTick','on','Yminortick','on','YTick',[0:0.2:1],'Ylim',[0 1]);
xlabel('Days post-infection','Fontsize',24);
ylabel('Diagnostic sensitivity','Fontsize',24);
box off;
legend({'RT-PCR','LumiraDx','Sofia','BinaxNOW','BD Veritor','CareStart'},'Fontsize',20,'NumColumns',3,'Location','NorthEast');
legend boxoff;
rmpath([pwd '\Delta_Variant']);
addpath([pwd '\Delta_Variant\Results']);
print(gcf,['Figure1.png'],'-dpng','-r600');
print(gcf,['Figure1.tiff'],'-dtiff','-r600');