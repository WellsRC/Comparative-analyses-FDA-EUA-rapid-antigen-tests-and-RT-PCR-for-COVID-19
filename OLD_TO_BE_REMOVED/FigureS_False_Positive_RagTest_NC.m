%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Calculations for the false positive
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
close all;
clear;
figure('units','normalized','outerposition',[0.05 0.05 1 1]);
subplot('Position',[0.091,0.12,0.847689075630252,0.861631205673759]);
pA=0.308;


load('RAgTest_PlotOrder.mat');
NumTest=length(testName);

f=[1:14];
R=zeros(length(f),NumTest+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Expected transmission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RT-PCR 1 day delay
 load(['1-day_Delay_Testing_Frequency_RTPCR_NatComm.mat']);
 R(:,1)=(1-pA).*RTotS+pA.*RTotA;
for ii=1:NumTest
    load(['Testing_Frequency_' testName{ii} '_NatComm.mat'],'RTotS','RTotA');
    R(:,ii+1)=(1-pA).*RTotS+pA.*RTotA;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Expected transmission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ProbFP=zeros(length(f),NumTest+1);
DurT=14; % Duration in which the testing is being conducted

for jj=1:length(f)
    [ProbFP(jj,1)]= CalcFalsePositive('RTPCR',DurT,f(jj));
end

for tt=1:NumTest
    for jj=1:length(f)
        [ProbFP(jj,tt+1)]= CalcFalsePositive(testName{tt},DurT,f(jj));
    end
end


pbt=plot(ProbFP(1:2,:),R(1:2,:),'LineWidth',2); hold on
pbt2=plot(ProbFP,R,'LineWidth',2); 

ii=0;

[C,MF]= ColourTests('RTPCR');
pbt(ii+1).Color=C;    
pbt(ii+1).Marker=MF{1}; 
pbt(ii+1).MarkerFaceColor=C; 
pbt(ii+1).MarkerSize=8; 
pbt2(ii+1).Color=C;  
scatter(ProbFP(:,ii+1),R(:,ii+1),250.*exp(-f./5),C,'filled',MF{1})

for ii=1:NumTest
    [C,MF]= ColourTests(testName{ii});
    
    if(length(C(:,1))>1)
        C=C(end,:);
    end
    pbt(ii+1).Color=C;    
    pbt(ii+1).Marker=MF{1}; 
    pbt(ii+1).MarkerFaceColor=C; 
    pbt(ii+1).MarkerSize=8; 
    pbt2(ii+1).Color=C;  
    scatter(ProbFP(:,ii+1),R(:,ii+1),250.*exp(-f./5),C,'filled',MF{1})
end

XX=[0.145 0.25 0.145 0.25 0.145 0.25 0.145 0.25 0.145 0.25 0.145 0.25 0.145 0.25]-0.038;
YY=[1.2 1.2 1.15 1.15 1.1 1.1 1.05 1.05 1 1 0.95 0.95 0.9 0.9]+0.225;
scatter((XX),YY,250.*exp(-f./5),[0.5 0.5 0.5],'filled')
for ii=1:14
   text((XX(ii)+0.01),YY(ii),[num2str(ii) '-day'],'color',[0.5 0.5 0.5],'fontsize',22);
end

xxtick=10.^[-4:0];
semilogx(10.^linspace(-6,0,1001),ones(1001,1),'k-.','LineWidth',1.5);
grid on;
ylabel({'Effective reproduction number'},'Fontsize',28);
xlabel('Probability of at least one false-postive over two weeks','Fontsize',28);
box off;
set(gca,'LineWidth',2,'Tickdir','out','Fontsize',26,'xtick',(xxtick),'xticklabel',{num2str(xxtick')},'xscale','log','YTick',[0:0.1:1.7]);
xlim([5*10^(-4) 0.5]);
ylim([0 1.7]);

legendtestname=cell(NumTest+1,1);
legendtestname(2:end)=testName;
legendtestname(1)={'RT-PCR (one-day delay)'};

legend(pbt,legendtestname,'Fontsize',16,'NumColumns',4,'Location','NorthEast');
legend boxoff

print(gcf,['False_positive_All_NatComm.png'],'-dpng','-r600');