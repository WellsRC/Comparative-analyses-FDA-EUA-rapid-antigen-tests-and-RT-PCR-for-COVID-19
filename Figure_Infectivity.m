clear;
% close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Hellewell et al
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


figure('units','normalized','outerposition',[0.2 0.2 0.6 0.6]);
subplot('Position',[0.097711267605634,0.174774774774775,0.886971830985915,0.803603603603605]);

addpath([pwd '\Delta_Variant']);
addpath([pwd '\Delta_Variant\Results']);
[~,~,~,ts,td] = BaselineParameters;
t=linspace(-10,25,1001);
plot(t-ts,ViralShedding_Symptomatic(t,td),'k','LineWidth',2); hold on;
rmpath([pwd '\Delta_Variant']);
addpath([pwd '\Delta_Variant\Results']);

addpath([pwd '\Non_Delta']);
addpath([pwd '\Non_Delta\Results']);
[~,~,~,ts,td] = BaselineParameters;
plot(t-ts,ViralShedding_Symptomatic(t,td),'r','LineWidth',2);
rmpath([pwd '\Non_Delta']);
addpath([pwd '\Non_Delta\Results']);


plot(zeros(101,1),linspace(0,1,101),':','color',[0.75 0.75 0.75],'LineWidth',2);
set(gca,'LineWidth',2,'tickdir','out','Fontsize',20,'XTick',[-6:10],'xlim',[-6 8],'XMinorTick','on','Yminortick','on','YTick',[0:0.05:0.25],'Ylim',[0 0.25]);
xlabel('Days since symptom onset','Fontsize',24);
ylabel('Relative infectivity','Fontsize',24);
box off;
legend({'Delta Variant','Non-Delta Variant'},'Fontsize',20,'Location','NorthEast');
legend boxoff;

print(gcf,['Figure_Infectivity.png'],'-dpng','-r600');