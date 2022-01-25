addpath([pwd '\Non_Delta']);
addpath([pwd '\Non_Delta\Results']);

load('RAgTest_PlotOrder.mat');
testName=sort(testName);
[~,~,~,ts,~] = BaselineParameters;

Ntest=length(testName);
Day=[1:60]';
Integration_Start=Day-1;
Integration_End=Day;

Test_Name=cell(Ntest+1,1);

Sensitivity_Day=zeros(Ntest+1,length(Day));

ii=1;
[betaRTPCR,betaAg]=ParameterCOVIDTest([],1);

Test_Name{ii}='RT-PCR';


for dd=1:length(Day)
    Sensitivity_Day(ii,dd)=integral(@(x)TestSensitivity(x,ts,[],betaRTPCR),Integration_Start(dd),Integration_End(dd));
end

for ii=1:Ntest
    
    [Test_Name{ii+1}] = AdjustedNames_Plotting(testName{ii});
    [betaRTPCR,betaAg]=ParameterCOVIDTest(testName{ii},1);
    for dd=1:length(Day)
        Sensitivity_Day(ii+1,dd)=integral(@(x)TestSensitivity(x,ts,betaAg,betaRTPCR),Integration_Start(dd),Integration_End(dd));
    end
end

T=table(Test_Name,Sensitivity_Day);

writetable(T,'Discrete_Sensitivity_Ashcroft_et_al.csv');