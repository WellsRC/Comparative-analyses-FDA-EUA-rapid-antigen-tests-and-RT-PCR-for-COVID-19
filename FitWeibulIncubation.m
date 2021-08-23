clear;
%Estimation of incubation period distribution of COVID-19 using disease onset forward time: a novel cross-sectional and forward follow-up study
%7·76 days (95% CI: 7·02-8·53), and mean is 8·29 days (95% CI: 7·67-8·9), the 90th percentile is 14·28 days (95% CI: 13·64-14·90), and the 99th percentile is 20·31 days (95% CI: 19·15- 21·47)

Day=[2.07 4.97 7.76 11.04 14.28 16.32 20.31 24.95];
Prc=[5 25 50 75 90 95 99 99.9];

opts= optimset('MaxIter',10^4,'MaxFunEvals',10^4,'TolFun',10^(-12),'TolX',10^(-12),'UseParallel',false,'Display','off'); 
p=fmincon(@(x)sum((wblinv(Prc./100,x(1),x(2))-Day).^2),[2 2],[],[],[],[],[0 0],[100 100],[],opts);
