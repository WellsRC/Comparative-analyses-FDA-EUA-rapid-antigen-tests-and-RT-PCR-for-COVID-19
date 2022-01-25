%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency testing: BD Veritor 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parpool(16); % Parallel pool
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEW Sensitivity curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('RAgTest_Name.mat','testName');
NumTests=length(testName);
% for TestN=1:NumTests
%     [~,~,R0,ts,td] = BaselineParameters;
% 
%     SelfIsolate=1; %If sympmatics self-isolate
%     R0S=R0;
%     R0A=R0;
%     RTotA=zeros(14,1);
%     RTotS=zeros(14,1);
% 
%     
%     [betaRTPCR,betaAg]=ParameterCOVIDTest(testName{TestN},1);
% 
%     testtype=cell(14,1);
%     timetoff=cell(14,1);
%     NT=zeros(14,1);
%     for dT=1:14
%         timetoff{dT}=[0:dT:(floor(td))];
%         NT(dT)=length([timetoff{dT}]);
%         temp=cell(NT(dT),1);
% 
%         for ii=1:length([timetoff{dT}])
%             temp{ii}=betaAg;
%         end
%         testtype{dT}=temp;
%     end
% 
%     parfor dT=1:14  
% 
%         [RS,RA] = SerialTesting(testtype{dT},[timetoff{dT}],R0S,R0A,ts,td,SelfIsolate,NT(dT),dT,betaRTPCR);
% 
%         RTotS(dT)=sum(RS);
%         RTotA(dT)=sum(RA);
%     end
%     save(['Testing_Frequency_' testName{TestN} '_DeltaVOC.mat']);
% end 

load('RTPCR_Parameter_Uncertainty_Serial.mat','betaRTPCRv');

NSS=length(betaRTPCRv(:,1));

for TestN=NumTests:-1:1
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
    save(['Testing_Frequency_' testName{TestN} '_DeltaVOC_Uncertainty.mat']);
end

