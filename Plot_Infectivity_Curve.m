clear;
close all;
figure('units','normalized','outerposition',[0 0 0.5 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Baseline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
subplot('Position',[0.1134,0.592256332320174,0.85,0.390000000000001]);
addpath([pwd '\Delta_Variant']);
addpath([pwd '\Delta_Variant\Results']);
[pA,IncubationI,R0,ts,td] = BaselineParameters;
tt=linspace(0,35,100001);
SelfIsolate=1;
VS1=InfectiousnessfromInfection(tt,R0,R0,pA,ts,td,SelfIsolate);
VNS1=InfectiousnessfromInfection(tt,R0,R0,pA,ts,td,0);

RS=integral(@(x)InfectiousnessfromInfection(x,R0,R0,pA,ts,td,SelfIsolate),0,inf);


patch([tt flip(tt)],[VS1 zeros(size(tt))],hex2rgb('#F5BE41'),'LineStyle','none','Facealpha',0.3); hold on
text(-1.534999999999991+ts,0.2, num2str(round(RS,1),'%2.1f'),'Fontsize',16,'HorizontalAlignment','center');

p3=plot(tt,VS1,'color',hex2rgb('#F5BE41'),'LineWidth',2); hold on
p1=plot(tt,VNS1,'color',hex2rgb('#2A3132'),'LineWidth',2); hold on


legend([p1 p3],{'No self-isolation',['Self-isolation and ' num2str(100*pA) '% asymptomatic']},'Fontsize',14);
legend boxoff;
box off;
xlabel('Days post-infection','Fontsize',18);
ylabel('Infectivity','Fontsize',18);

ylim([0 0.8]);
xlim([0 16]);
set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'XTick',[0:25],'XMinorTick','on','Yminortick','on','YTick',[0:0.1:0.9]);
text(-1.779975278121137,0.8,'a','Fontsize',32,'Fontweight','bold');

rmpath([pwd '\Delta_Variant']);
rmpath([pwd '\Delta_Variant\Results']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Non-Delta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
subplot('Position',[0.1134,0.114,0.85,0.390000000000001]);

addpath([pwd '\Non_Delta']);
addpath([pwd '\Non_Delta\Results']);
[pA,IncubationI,R0,ts,td] = BaselineParameters;
tt=linspace(0,35,100001);
SelfIsolate=1;
VS2=InfectiousnessfromInfection(tt,R0,R0,pA,ts,td,SelfIsolate);
VNS2=InfectiousnessfromInfection(tt,R0,R0,pA,ts,td,0);

RS=integral(@(x)InfectiousnessfromInfection(x,R0,R0,pA,ts,td,SelfIsolate),0,inf);


patch([tt flip(tt)],[VS2 zeros(size(tt))],hex2rgb('#F5BE41'),'LineStyle','none','Facealpha',0.3); hold on
text(-1.534999999999991+ts,0.2, num2str(round(RS,1),'%2.1f'),'Fontsize',16,'HorizontalAlignment','center');

p3=plot(tt,VS2,'color',hex2rgb('#F5BE41'),'LineWidth',2); hold on
p1=plot(tt,VNS2,'color',hex2rgb('#2A3132'),'LineWidth',2); hold on


legend([p1 p3],{'No self-isolation',['Self-isolation and ' num2str(100*pA) '% asymptomatic']},'Fontsize',14);
legend boxoff;
box off;
xlabel('Days post-infection','Fontsize',18);
ylabel('Infectivity','Fontsize',18);

ylim([0 0.8]);
xlim([0 16]);
set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'XTick',[0:25],'XMinorTick','on','Yminortick','on','YTick',[0:0.1:0.9]);
text(-1.779975278121137,0.8,'b','Fontsize',32,'Fontweight','bold');

rmpath([pwd '\Non_Delta']);
rmpath([pwd '\Non_Delta\Results']);

print(gcf,['Infectivity_Curves.png'],'-dpng','-r600');