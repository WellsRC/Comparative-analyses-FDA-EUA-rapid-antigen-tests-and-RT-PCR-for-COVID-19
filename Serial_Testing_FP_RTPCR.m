
addpath([pwd '\Delta_Variant']);
addpath([pwd '\Delta_Variant\Results']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Determine the minimal testing frequency
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    ft=[1:14];
    pA=0.351;
    DurT=14;
    load([num2str(1) '-day_Delay_Testing_Frequency_RTPCR_DeltaVOC.mat'],'RTotA','RTotS')
    R=(1-pA).*RTotS+pA.*RTotA;
    FreqT(1)=max(ft(R<1)); % take the maximum time between tests
    
    load([num2str(1) '-day_Delay_Testing_Frequency_RTPCR_DeltaVOC_Uncertainty.mat'],'RTotAv','RTotSv')
    R=(1-pA).*RTotSv+pA.*RTotAv;
    TestF=zeros(length(R(1,:)),1);
    FPU=zeros(length(R(1,:)),250);
    for ii=1:length(TestF)
        if(~isempty(ft(R(:,ii)<1)))
            TestF(ii)=max(ft(R(:,ii)<1)); 
        else
            TestF(ii)=1;
        end
       parfor jj=1:250
            FPU(ii,jj)=CalcFalsePositive('RT-PCR',DurT,TestF(ii),0);
       end
    end
    
    [~,FreqT(2),FreqT(3)]=Credible_Interval_High_Density(FreqT(1),TestF,0.95,'discrete',[1 14]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Risk of false positve
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    FPProb(1)=CalcFalsePositive('RT-PCR',DurT,FreqT(1),1); % Calacute the probability of a 
    [~,FPProb(2),FPProb(3)]=Credible_Interval_High_Density(FPProb(1),FPU(:),0.95,'continuous',[0 1]);
    
    rmpath([pwd '\Delta_Variant']);
rmpath([pwd '\Delta_Variant\Results']);