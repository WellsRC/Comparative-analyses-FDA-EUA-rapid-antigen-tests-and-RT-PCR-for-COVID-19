%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Calculations for the false positive
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
close all;
clear;


addpath([pwd '\Delta_Variant']);
addpath([pwd '\Delta_Variant\Results']);

[pA,~,R0,ts,td] = BaselineParameters;

figure('units','normalized','outerposition',[0.05 0.05 1 1]);
subplot('Position',[0.091,0.12,0.847689075630252,0.861631205673759]);


load('RAgTest_PlotOrder.mat');
NumTest=length(testName);

f=[1:14];
R=zeros(length(f),NumTest+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Expected transmission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RT-PCR 1 day delay
 load(['1-day_Delay_Testing_Frequency_RTPCR_DeltaVOC.mat']);
 R(:,1)=(1-pA).*(RTotS)+pA.*(RTotA);
 
testname2=cell((NumTest+1),1);
testname2(1)={'RT-PCR (24-h delay)'};
for ii=1:NumTest
    load(['Testing_Frequency_' testName{ii} '_DeltaVOC.mat'],'RTotS','RTotA');
    R(:,ii+1)=(1-pA).*(RTotS)+pA.*(RTotA);
    testname2(ii+1)={AdjustedNames_Plotting(testName{ii})};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Expected transmission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ProbFP=zeros(length(f),NumTest+1);
DurT=14; % Duration in which the testing is being conducted

for jj=1:length(f)
    [ProbFP(jj,1)]= CalcFalsePositive('RTPCR',DurT,f(jj),1);
end

for tt=1:NumTest
    for jj=1:length(f)
        [ProbFP(jj,tt+1)]= CalcFalsePositive(testName{tt},DurT,f(jj),1);
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

dx=0.025;
Mxx=0.027;
XX=[Mxx Mxx+dx Mxx Mxx+dx Mxx Mxx+dx Mxx Mxx+dx Mxx Mxx+dx Mxx Mxx+dx Mxx Mxx+dx];
dy=0.075;
Mxy=2.7;
YY=[Mxy Mxy Mxy-dy Mxy-dy Mxy-2.*dy Mxy-2.*dy Mxy-3.*dy Mxy-3.*dy Mxy-4.*dy Mxy-4.*dy Mxy-5.*dy Mxy-5.*dy Mxy-6.*dy Mxy-6.*dy];
scatter((XX),YY,250.*exp(-f./5),[0.5 0.5 0.5],'filled')
for ii=1:14
   text((XX(ii)+dx./10),YY(ii),[num2str(ii) '-day'],'color',[0.5 0.5 0.5],'fontsize',22);
end

xxtick=10.^[-4:0];
semilogx(10.^linspace(-6,0,1001),ones(1001,1),'k-.','LineWidth',1.5);
grid on;
ylabel({'Effective reproduction number'},'Fontsize',28);
xlabel('Probability of at least one false-positive over two weeks','Fontsize',28);
box off;
set(gca,'LineWidth',2,'Tickdir','out','Fontsize',26,'xtick',(xxtick),'xticklabel',{num2str(xxtick')},'xscale','log','YTick',[0:0.25:2.25]);
xlim([5*10^(-4) 0.5]);
ylim([0 2.75]);


legend(pbt,testname2,'Fontsize',16,'NumColumns',4,'Location','NorthWest');
legend boxoff
print(gcf,['False_positive_All_Baseline.png'],'-dpng','-r600');
[r,p]=corr(ProbFP(:),R(:));
fprintf('Pearson Correlation between False positve and RE: r= %4.3f p=%3.2e \n',[r,p]);
rmpath([pwd '\Delta_Variant']);
rmpath([pwd '\Delta_Variant\Results']);