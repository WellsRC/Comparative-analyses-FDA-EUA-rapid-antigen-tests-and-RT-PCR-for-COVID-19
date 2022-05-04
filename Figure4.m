%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Calculations for the false positive
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
close all;
clear;
figure('units','normalized','outerposition',[0.05 0.05 1 1]);
subplot('Position',[0.091,0.12,0.847689075630252,0.861631205673759]);

addpath([pwd '\Delta_Variant']);
addpath([pwd '\Delta_Variant\Results']);
[pA,~,~,~,~] = BaselineParameters;

f=[1:14];
R=zeros(length(f),11);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Expected transmission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([num2str(0) '-day_Delay_Testing_Frequency_RTPCR_DeltaVOC.mat']);
R(:,1)=(1-pA).*RTotS+pA.*RTotA;

load('Testing_Frequency_LumiraDX (Anterior Nasal Swab)_DeltaVOC.mat')
R(:,2)=(1-pA).*RTotS+pA.*RTotA;

load('Testing_Frequency_Sofia (FDA)_DeltaVOC.mat')
R(:,3)=(1-pA).*RTotS+pA.*RTotA;

load('Testing_Frequency_BinaxNOW (FDA)_DeltaVOC.mat')
R(:,4)=(1-pA).*RTotS+pA.*RTotA;

load('Testing_Frequency_BD Veritor_DeltaVOC.mat')
R(:,5)=(1-pA).*RTotS+pA.*RTotA;

load('Testing_Frequency_CareStart (Anterior Nasal Swab - FDA)_DeltaVOC.mat')
R(:,6)=(1-pA).*RTotS+pA.*RTotA;


MF={'-p','-.s','-.d','-.^','-.o','-.h'};

for dd=1:5
    load([num2str(dd) '-day_Delay_Testing_Frequency_RTPCR_DeltaVOC.mat']);
    R(:,dd+6)=(1-pA).*RTotS+pA.*RTotA;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Expected transmission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ProbFP=zeros(length(f),11);
DurT=14; % Duration in which the testing is being conducted

NumTestSerial=ceil(DurT./f); % Number of tests conducted over the time period

for jj=1:length(f)
    [ProbFP(jj,1)]= CalcFalsePositive('RTPCR',DurT,f(jj),1);
    [ProbFP(jj,2)]=CalcFalsePositive('LumiraDX (Anterior Nasal Swab)',DurT,f(jj),1);
    [ProbFP(jj,3)]=CalcFalsePositive('Sofia (FDA)',DurT,f(jj),1);
    [ProbFP(jj,4)]=CalcFalsePositive('BinaxNOW (FDA)',DurT,f(jj),1);
    [ProbFP(jj,5)]=CalcFalsePositive('BD Veritor',DurT,f(jj),1);
    [ProbFP(jj,6)]=CalcFalsePositive('CareStart (Anterior Nasal Swab - FDA)',DurT,f(jj),1);

    for dd=1:5
        [ProbFP(jj,dd+6)]=CalcFalsePositive('RTPCR',DurT,f(jj),1);
    end
end


CTest=[hex2rgb('#231B12');hex2rgb('#bfd3e6');hex2rgb('#41ae76');hex2rgb('#F9A603');hex2rgb('#993404');hex2rgb('#88419d'); hex2rgb('525252');hex2rgb('737373');hex2rgb('969696');hex2rgb('bdbdbd');hex2rgb('d9d9d9')]; 
hold on;
pbt=plot(ProbFP(1,1),R(1,1),ProbFP(1,2),R(1,2),ProbFP(1,3),R(1,3),ProbFP(1,4),R(1,4),ProbFP(1,5),R(1,5),ProbFP(1,6),R(1,6),ProbFP(1,7),R(1,7),ProbFP(1,8),R(1,8),ProbFP(1,9),R(1,9),ProbFP(1,10),R(1,10),ProbFP(1,11),R(1,11),'LineWidth',2);
MF={'p','s','d','^','o','h','p','p','p','p','p'};
for ii=1:11
   pbt(ii).Color=CTest(ii,:); 
   pbt(ii).Marker=MF{ii};
   pbt(ii).MarkerFaceColor=CTest(ii,:); 
   pbt(ii).MarkerSize=9;
end
pbt2=plot(ProbFP,R,'LineWidth',2);
for ii=1:11
   pbt2(ii).Color=CTest(ii,:); 
end
scatter(ProbFP(:,1),R(:,1),250.*exp(-f./5),CTest(1,:),'filled','p')
scatter(ProbFP(:,2),R(:,2),250.*exp(-f./5),CTest(2,:),'filled','s')
scatter(ProbFP(:,3),R(:,3),250.*exp(-f./5),CTest(3,:),'filled','d')
scatter(ProbFP(:,4),R(:,4),250.*exp(-f./5),CTest(4,:),'filled','^')
scatter(ProbFP(:,5),R(:,5),250.*exp(-f./5),CTest(5,:),'filled')
scatter(ProbFP(:,6),R(:,6),250.*exp(-f./6),CTest(6,:),'filled','h')
for dd=1:5
    scatter(ProbFP(:,6+dd),R(:,6+dd),250.*exp(-f./5),CTest(6+dd,:),'filled','p')
end

XX=[0.145 0.25 0.145 0.25 0.145 0.25 0.145 0.25 0.145 0.25 0.145 0.25 0.145 0.25];
YY=(2.8/1.7).*([1.2 1.2 1.15 1.15 1.1 1.1 1.05 1.05 1 1 0.95 0.95 0.9 0.9]+0.3);
scatter((XX),YY,250.*exp(-f./5),[0.5 0.5 0.5],'filled')
for ii=1:14
   text((XX(ii)+0.01),YY(ii),[num2str(ii) '-day'],'color',[0.5 0.5 0.5],'fontsize',22);
end
xxtick=10.^[-4:0];
semilogx(10.^linspace(-6,0,1001),ones(1001,1),'k-.','LineWidth',1.5);
grid on;
ylabel({'Effective reproduction number'},'Fontsize',28);
xlabel('Probability of at least one false-positive over two weeks','Fontsize',28);
box off;
legend({'RT-PCR (No delay)','LumiraDx','Sofia','BinaxNOW','BD Veritor','CareStart','RT-PCR 1-day delay','RT-PCR 2-day delay','RT-PCR 3-day delay','RT-PCR 4-day delay','RT-PCR 5-day delay'},'Fontsize',22,'NumColumns',4,'Location','NorthEast');
legend boxoff;
set(gca,'LineWidth',2,'Tickdir','out','Fontsize',26,'xtick',(xxtick),'xticklabel',{num2str(xxtick')},'xscale','log','YTick',[0:0.2:2.8],'Yminortick','on');
xlim([7.5*10^(-4) 0.5]);
ylim([0 2.8]);


rmpath([pwd '\Delta_Variant']);
rmpath([pwd '\Delta_Variant\Results']);

print(gcf,['Figure4.png'],'-dpng','-r600');
print(gcf,['Figure4.tiff'],'-dtiff','-r600');