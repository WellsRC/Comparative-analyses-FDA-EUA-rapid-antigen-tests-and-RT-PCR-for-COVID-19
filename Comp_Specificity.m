clear;
clc;
x=[2196 212; 3 1];

% Column 1 is the emperical testing data
% Column 2 is the FDA testing data
% row 1 is number of Negative agreement (RA=- and RTPCR -)
% row 2 is the number of disagreements (RA=+ and RT-PCR=-);

[h p]=fishertest(x)