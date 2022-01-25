close all;
clear;
clc;
addpath([pwd '\Non_Delta']);
addpath([pwd '\Non_Delta\Results']);

It=[0:25];
[I1,I2]=meshgrid(It,It);
I1=I1(:);
I2=I2(:);
[betaRTPCR,betaAg]=ParameterCOVIDTest('BD Veritor',1);
[SAg] = Test_Specificity ('BD Veritor',1);
[SR] = Test_Specificity ('RT-PCR',1);
    testtype=cell(3,1);
    testtype{1}=[];
    testtype{2}=betaAg;
    testtype{3}=[];
% TestTimes=[2 5 8];
% FPE=[0 0 0];
% TPE=[0 0 0];
% N=[124 121 96];

TestTimes=[3 6 9];
FPE=[2 0 0];
TPE=[0 3 0];
N=[458 458 457];

q=1;


[~,~,R0,ts,td] = BaselineParameters;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% SPECIFY WHAT R0 Should be
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55




C1=integral(@(x)(1-TestSensitivity_Case_Estimate(x,ts,testtype{1},betaRTPCR,td,SR)),0,td-q); % Inividual needs to be infectious after their release from quarantine (i.e. the 24 h delay from obtaining test results)
% False positve
FI1=(1./C1).*integral(@(x)((1-TestSensitivity_Case_Estimate(x,ts,testtype{1},betaRTPCR,td,SR)).*(TestSensitivity_Case_Estimate(x+q+TestTimes(1),ts,testtype{2},betaRTPCR,td,SAg)).*(1-TestSensitivity_Case_Estimate(x+q+TestTimes(1),ts,testtype{3},betaRTPCR,td,SR))),0,td-q);
FU1=(1-SAg).*SR;

FI2=(1./C1).*integral(@(x)((1-TestSensitivity_Case_Estimate(x,ts,testtype{1},betaRTPCR,td,SR)).*(TestSensitivity_Case_Estimate(x+q+TestTimes(2),ts,testtype{2},betaRTPCR,td,SAg)).*(1-TestSensitivity_Case_Estimate(x+q+TestTimes(2),ts,testtype{3},betaRTPCR,td,SR))),0,td-q);
FU2=(1-SAg).*SR;

FI3=(1./C1).*integral(@(x)((1-TestSensitivity_Case_Estimate(x,ts,testtype{1},betaRTPCR,td,SR)).*(TestSensitivity_Case_Estimate(x+q+TestTimes(3),ts,testtype{2},betaRTPCR,td,SAg)).*(1-TestSensitivity_Case_Estimate(x+q+TestTimes(3),ts,testtype{3},betaRTPCR,td,SR))),0,td-q);
FU3=(1-SAg).*SR;
% Confirmed positve

PU1=(1-SAg).*(1-SR);
PI1=(1./C1).*integral(@(x)((1-TestSensitivity_Case_Estimate(x,ts,testtype{1},betaRTPCR,td,SR)).*(TestSensitivity_Case_Estimate(x+q+TestTimes(1),ts,testtype{2},betaRTPCR,td,SAg)).*(TestSensitivity_Case_Estimate(x+q+TestTimes(1),ts,testtype{3},betaRTPCR,td,SR))),0,td-q);


PU2=(1-SAg).*(1-SR);
PI2=(1./C1).*integral(@(x)((1-TestSensitivity_Case_Estimate(x,ts,testtype{1},betaRTPCR,td,SR)).*(TestSensitivity_Case_Estimate(x+q+TestTimes(2),ts,testtype{2},betaRTPCR,td,SAg)).*(TestSensitivity_Case_Estimate(x+q+TestTimes(2),ts,testtype{3},betaRTPCR,td,SR))),0,td-q);


PU3=(1-SAg).*(1-SR);
PI3=(1./C1).*integral(@(x)((1-TestSensitivity_Case_Estimate(x,ts,testtype{1},betaRTPCR,td,SR)).*(TestSensitivity_Case_Estimate(x+q+TestTimes(3),ts,testtype{2},betaRTPCR,td,SAg)).*(TestSensitivity_Case_Estimate(x+q+TestTimes(3),ts,testtype{3},betaRTPCR,td,SR))),0,td-q);

% False negative from RA test

RNI1=(1./C1).*integral(@(x)((1-TestSensitivity_Case_Estimate(x,ts,testtype{1},betaRTPCR,td,SR)).*(1-TestSensitivity_Case_Estimate(x+q+TestTimes(1),ts,testtype{2},betaRTPCR,td,SAg))),0,td-q);
RNI2=(1./C1).*integral(@(x)((1-TestSensitivity_Case_Estimate(x,ts,testtype{1},betaRTPCR,td,SR)).*(1-TestSensitivity_Case_Estimate(x+q+TestTimes(2),ts,testtype{2},betaRTPCR,td,SAg))),0,td-q);
RNI3=(1./C1).*integral(@(x)((1-TestSensitivity_Case_Estimate(x,ts,testtype{1},betaRTPCR,td,SR)).*(1-TestSensitivity_Case_Estimate(x+q+TestTimes(3),ts,testtype{2},betaRTPCR,td,SAg))),0,td-q);

% Prob PQT

testtype=cell(1,1);
testtype{1}=[];

RTPCR_Folllow=[0];

PQT=zeros(4,1);
PQT(1)=(1./C1).*integral2(@(x,t)(InfectiousnessfromInfectionTesting_EstimateOffshore(t+x,x,[0],testtype,RTPCR_Folllow,R0,R0,1,ts,td,0,betaRTPCR)),0,td-q,q,q+TestTimes(1));

testtype=cell(2,1);
testtype{1}=[];
testtype{2}=[betaAg];

RTPCR_Folllow=[0 1];

PQT(2)=(1./C1).*integral2(@(x,t)(InfectiousnessfromInfectionTesting_EstimateOffshore(t+x,x,[0 q+TestTimes(1)],testtype,RTPCR_Folllow,R0,R0,1,ts,td,0,betaRTPCR)),0,td-q,q+TestTimes(1),q+TestTimes(2));

testtype=cell(3,1);
testtype{1}=[];
testtype{2}=[betaAg];
testtype{3}=[betaAg];

RTPCR_Folllow=[0 1 1];

PQT(3)=(1./C1).*integral2(@(x,t)(InfectiousnessfromInfectionTesting_EstimateOffshore(t+x,x,[0 q+TestTimes(1) q+TestTimes(2)],testtype,RTPCR_Folllow,R0,R0,1,ts,td,0,betaRTPCR)),0,td-q,q+TestTimes(2),q+TestTimes(3));

testtype=cell(4,1);
testtype{1}=[];
testtype{2}=[betaAg];
testtype{3}=[betaAg];
testtype{4}=[betaAg];

RTPCR_Folllow=[0 1 1 1];

PQT(4)=(1./C1).*integral2(@(x,t)(InfectiousnessfromInfectionTesting_EstimateOffshore(t+x,x,[0 q+TestTimes(1) q+TestTimes(2) q+TestTimes(3)],testtype,RTPCR_Folllow,R0,R0,1,ts,td,0,betaRTPCR)),0,td-q,q+TestTimes(3),inf);

PQT=Probability_Onward(sum(PQT),1);


I=I1;
L1=poisspdf(FPE(1),FU1.*(N(1)-I)+FI1.*I).*poisspdf(TPE(1),PU1.*(N(1)-I)+PI1.*I).*(RNI1.^(I-(FI1+PI1).*I));

L2=poisspdf(FPE(2),FU2.*(N(2)-I.*N(2)/N(1))+FI2.*I.*N(2)/N(1)).*poisspdf(TPE(2),PU2.*(N(2)-I.*N(2)/N(1))+PI2.*I.*N(2)/N(1)).*(RNI2.^(I.*N(2)/N(1)-(FI2+PI2).*I.*N(2)/N(1)));

L3=poisspdf(FPE(3),FU3.*(N(3)-I.*N(3)/N(1))+FI3.*I.*N(3)/N(1)).*poisspdf(TPE(3),PU3.*(N(3)-I.*N(3)/N(1))+PI3.*I.*N(3)/N(1)).*(RNI3.^(I.*N(3)/N(1)-(FI3+PI3).*I.*N(3)/N(1)));

LT1=L1.*L2.*L3.*(1-PQT).^I;

FPE=[0 0 0];
TPE=[0 0 0];
N=[124 121 96];


q=1;


[~,~,R0,ts,td] = BaselineParameters;
C1=integral(@(x)(1-TestSensitivity_Case_Estimate(x,ts,testtype{1},betaRTPCR,td,SR)),0,td-q); % Inividual needs to be infectious after their release from quarantine (i.e. the 24 h delay from obtaining test results)
% False positve
FI1=(1./C1).*integral(@(x)((1-TestSensitivity_Case_Estimate(x,ts,testtype{1},betaRTPCR,td,SR)).*(TestSensitivity_Case_Estimate(x+q+TestTimes(1),ts,testtype{2},betaRTPCR,td,SAg)).*(1-TestSensitivity_Case_Estimate(x+q+TestTimes(1),ts,testtype{3},betaRTPCR,td,SR))),0,td-q);
FU1=(1-SAg).*SR;

FI2=(1./C1).*integral(@(x)((1-TestSensitivity_Case_Estimate(x,ts,testtype{1},betaRTPCR,td,SR)).*(TestSensitivity_Case_Estimate(x+q+TestTimes(2),ts,testtype{2},betaRTPCR,td,SAg)).*(1-TestSensitivity_Case_Estimate(x+q+TestTimes(2),ts,testtype{3},betaRTPCR,td,SR))),0,td-q);
FU2=(1-SAg).*SR;

FI3=(1./C1).*integral(@(x)((1-TestSensitivity_Case_Estimate(x,ts,testtype{1},betaRTPCR,td,SR)).*(TestSensitivity_Case_Estimate(x+q+TestTimes(3),ts,testtype{2},betaRTPCR,td,SAg)).*(1-TestSensitivity_Case_Estimate(x+q+TestTimes(3),ts,testtype{3},betaRTPCR,td,SR))),0,td-q);
FU3=(1-SAg).*SR;
% Confirmed positve

PU1=(1-SAg).*(1-SR);
PI1=(1./C1).*integral(@(x)((1-TestSensitivity_Case_Estimate(x,ts,testtype{1},betaRTPCR,td,SR)).*(TestSensitivity_Case_Estimate(x+q+TestTimes(1),ts,testtype{2},betaRTPCR,td,SAg)).*(TestSensitivity_Case_Estimate(x+q+TestTimes(1),ts,testtype{3},betaRTPCR,td,SR))),0,td-q);


PU2=(1-SAg).*(1-SR);
PI2=(1./C1).*integral(@(x)((1-TestSensitivity_Case_Estimate(x,ts,testtype{1},betaRTPCR,td,SR)).*(TestSensitivity_Case_Estimate(x+q+TestTimes(2),ts,testtype{2},betaRTPCR,td,SAg)).*(TestSensitivity_Case_Estimate(x+q+TestTimes(2),ts,testtype{3},betaRTPCR,td,SR))),0,td-q);


PU3=(1-SAg).*(1-SR);
PI3=(1./C1).*integral(@(x)((1-TestSensitivity_Case_Estimate(x,ts,testtype{1},betaRTPCR,td,SR)).*(TestSensitivity_Case_Estimate(x+q+TestTimes(3),ts,testtype{2},betaRTPCR,td,SAg)).*(TestSensitivity_Case_Estimate(x+q+TestTimes(3),ts,testtype{3},betaRTPCR,td,SR))),0,td-q);

% False negative from RA test

RNI1=(1./C1).*integral(@(x)((1-TestSensitivity_Case_Estimate(x,ts,testtype{1},betaRTPCR,td,SR)).*(1-TestSensitivity_Case_Estimate(x+q+TestTimes(1),ts,testtype{2},betaRTPCR,td,SAg))),0,td-q);
RNI2=(1./C1).*integral(@(x)((1-TestSensitivity_Case_Estimate(x,ts,testtype{1},betaRTPCR,td,SR)).*(1-TestSensitivity_Case_Estimate(x+q+TestTimes(2),ts,testtype{2},betaRTPCR,td,SAg))),0,td-q);
RNI3=(1./C1).*integral(@(x)((1-TestSensitivity_Case_Estimate(x,ts,testtype{1},betaRTPCR,td,SR)).*(1-TestSensitivity_Case_Estimate(x+q+TestTimes(3),ts,testtype{2},betaRTPCR,td,SAg))),0,td-q);

% Prob PQT

testtype=cell(1,1);
testtype{1}=[];

RTPCR_Folllow=[0];

PQT=zeros(4,1);
PQT(1)=(1./C1).*integral2(@(x,t)(InfectiousnessfromInfectionTesting_EstimateOffshore(t+x,x,[0],testtype,RTPCR_Folllow,R0,R0,1,ts,td,0,betaRTPCR)),0,td-q,q,q+TestTimes(1));

testtype=cell(2,1);
testtype{1}=[];
testtype{2}=[betaAg];

RTPCR_Folllow=[0 1];

PQT(2)=(1./C1).*integral2(@(x,t)(InfectiousnessfromInfectionTesting_EstimateOffshore(t+x,x,[0 q+TestTimes(1)],testtype,RTPCR_Folllow,R0,R0,1,ts,td,0,betaRTPCR)),0,td-q,q+TestTimes(1),q+TestTimes(2));

testtype=cell(3,1);
testtype{1}=[];
testtype{2}=[betaAg];
testtype{3}=[betaAg];

RTPCR_Folllow=[0 1 1];

PQT(3)=(1./C1).*integral2(@(x,t)(InfectiousnessfromInfectionTesting_EstimateOffshore(t+x,x,[0 q+TestTimes(1) q+TestTimes(2)],testtype,RTPCR_Folllow,R0,R0,1,ts,td,0,betaRTPCR)),0,td-q,q+TestTimes(2),q+TestTimes(3));

testtype=cell(4,1);
testtype{1}=[];
testtype{2}=[betaAg];
testtype{3}=[betaAg];
testtype{4}=[betaAg];

RTPCR_Folllow=[0 1 1 1];

PQT(4)=(1./C1).*integral2(@(x,t)(InfectiousnessfromInfectionTesting_EstimateOffshore(t+x,x,[0 q+TestTimes(1) q+TestTimes(2) q+TestTimes(3)],testtype,RTPCR_Folllow,R0,R0,1,ts,td,0,betaRTPCR)),0,td-q,q+TestTimes(3),inf);

PQT=Probability_Onward(sum(PQT),1);


I=I2;
L1=poisspdf(FPE(1),FU1.*(N(1)-I)+FI1.*I).*poisspdf(TPE(1),PU1.*(N(1)-I)+PI1.*I).*(RNI1.^(I-(FI1+PI1).*I));

L2=poisspdf(FPE(2),FU2.*(N(2)-I.*N(2)/N(1))+FI2.*I.*N(2)/N(1)).*poisspdf(TPE(2),PU2.*(N(2)-I.*N(2)/N(1))+PI2.*I.*N(2)/N(1)).*(RNI2.^(I.*N(2)/N(1)-(FI2+PI2).*I.*N(2)/N(1)));

L3=poisspdf(FPE(3),FU3.*(N(3)-I.*N(3)/N(1))+FI3.*I.*N(3)/N(1)).*poisspdf(TPE(3),PU3.*(N(3)-I.*N(3)/N(1))+PI3.*I.*N(3)/N(1)).*(RNI3.^(I.*N(3)/N(1)-(FI3+PI3).*I.*N(3)/N(1)));

LT2=L1.*L2.*L3.*(1-PQT).^I;


% figure('units','normalized','outerposition',[0 0 1 1]);
% 
 L=LT1.*LT2;
% 
% 
% contourf(It,It,reshape(log(L),length(It),length(It)),'LineStyle','none');

RWE=3./(I1+I2);
RWE(I1+I2==0)=1;
RWE(RWE>1)=1;
PPr=L./sum(L);
RWEu=unique(RWE);
PP=zeros(size(RWEu));
for ii=1:length(RWEu)
    f=find(RWE==RWEu(ii));
    PP(ii)=sum(PPr(f));
end

figure('units','normalized','outerposition',[0 0.1 1 0.65]);


ITot=I1+I2;
Iu=unique(ITot);
PPr=L./sum(L);
PPI=zeros(size(Iu));
for ii=1:length(Iu)
    f=find(ITot==Iu(ii));
    PPI(ii)=sum(PPr(f));
end

subplot('Position',[0.07 0.16 0.4 0.8])
b=bar(Iu,PPI,'k','LineStyle','none'); hold on
p=plot(3.*ones(101,1),linspace(0,0.2,101),'r-.','LineWidth',2);
legend([b p],{'Estimated distribution','Observed RT-PCR confirmed'},'Fontsize',24);
legend boxoff;
xlabel('Number of expected cases offshore');
ylabel('Density','Fontsize',24)
box off;
xlim([0 20])
ylim([0 0.2]);
text(-2.755,0.2,'A','Fontsize',32,'FontWeight','bold');
set(gca,'LineWidth',2,'tickdir','out','Fontsize',22,'YTick',[0:0.05:0.2],'YMinortick','on','XTick',[0:5:20],'Xminortick','on');
%%%% 
% Determine the credible interval
%%%
alphaZ=0.95;
    z=flip(unique(PPI));
    z=z(z>0);
    AUC=zeros(size(z));
    for kk=1:length(z)
       AUC(kk)=sum(PPI(PPI>=z(kk)))-alphaZ;
    end
    zmin=z(AUC==min(AUC(AUC>=0)));
    lb=min(Iu(PPI>=zmin));
    ub=max(Iu(PPI>=zmin));
%%%%%%%%%%%%%%%%%%%%%

fprintf('Number of infected individuals inferred from anaylsys: %2.0f (%2.0f-%2.0f) \n',[Iu(PPI==max(PPI)) lb ub]);


subplot('Position',[0.58 0.16 0.4 0.8])
scatter(RWEu,PP./integral(@(t)pchip(RWEu,PP,t),0,1),40','k','filled'); hold on;
plot(linspace(0,1,101),pchip(RWEu,PP,linspace(0,1,101))./integral(@(t)pchip(RWEu,PP,t),0,1),'color',[0.7 0.7 0.7],'LineWidth',2);
legend({'Discrete','Continuous approximation'},'Fontsize',24,'Location','NorthWest');
legend boxoff;
xlabel('Real world effectiveness','Fontsize',24);
ylabel('Density','Fontsize',24)

set(gca,'LineWidth',2,'tickdir','out','Fontsize',22,'YTick',[0:0.25:2.5],'YMinortick','on','XTick',[0:0.1:1],'Xminortick','on');
ylim([0 2.5]);
text(-0.1457,2.5,'B','Fontsize',32,'FontWeight','bold');
box off;
%%%%%%%%%%%%%%%%%%%%%%
% Determine credible interval
%%%%%%%%%%%%%%%%%%%%%%%
mle=3/Iu(PPI==max(PPI));
bpf=fminbnd(@(bp)(integral(@(t)pchip(RWEu,PP,t),max(mle-bp,0),min(mle+bp,1))./integral(@(t)pchip(RWEu,PP,t),0,1)-alphaZ).^2,0,mle);
lb=max(mle-bpf,0);
ub=min(mle+bpf,1);
% [bnd]=fminbnd(@(bp)(pchip(RWEu,PP,bp)-pdf(pd,icdf(pd,min(cdf(pd,bp)+alphaZ,1)))).^2+((cdf(pd,icdf(pd,min(cdf(pd,bp)+alphaZ,1)))-cdf(pd,bp))-alphaZ).^2,databnd(1),mle,options);
% 
%         lb=bnd;
%         ub=icdf(pd,min(cdf(pd,bnd)+alphaZ,1));
%         

fprintf('Effectiveness of the serial testing strategy: %3.1f%% (%3.1f%%-%3.1f%%) \n',100.*[3/Iu(PPI==max(PPI)) lb ub]);

rmpath([pwd '\Non_Delta']);
rmpath([pwd '\Non_Delta\Results']);
print(gcf,['OffShore_Infected_Cases.png'],'-dpng','-r600');