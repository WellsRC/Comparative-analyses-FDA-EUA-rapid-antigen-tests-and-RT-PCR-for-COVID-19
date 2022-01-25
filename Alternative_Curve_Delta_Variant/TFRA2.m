%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency testing: BD Veritor 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parpool(16); % Parallel pool
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEW Sensitivity curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('RAgTest_Name.mat','testName');
NumTests=length(testName);
load('RTPCR_Parameter_Uncertainty_Serial.mat','betaRTPCRv');

NSS=length(betaRTPCRv(:,1));

for TestN=5:7
    [~,~,R0,ts,td] = BaselineParameters;

    SelfIsolate=1; %If sympmatics self-isolate
    R0S=R0;
    R0A=R0;
    RTotA=zeros(1,NSS);
    RTotS=zeros(1,NSS);

    RTotAv=zeros(14,NSS);
    RTotSv=zeros(14,NSS);
    
    

    for dT=1:14  
        parfor ns=1:NSS
            [~,betaAg]=ParameterCOVIDTest(testName{TestN},1);

                timetoff=[0:dT:(floor(td))];
                NT=length([timetoff]);
                temp=cell(NT,1);

                for ii=1:length([timetoff])
                    temp{ii}=betaAg;
                end
                testtype=temp;
            [RS,RA] = SerialTesting(testtype,[timetoff],R0S,R0A,ts,td,SelfIsolate,NT,dT,betaRTPCRv(ns,:));

            RTotS(ns)=sum(RS);
            RTotA(ns)=sum(RA);
        end
        RTotAv(dT,:)=RTotA;
        RTotSv(dT,:)=RTotS;
    end
    save(['Testing_Frequency_' testName{TestN} '_DeltaVOC_Alternative_PCR_Uncertainty.mat']);
end

