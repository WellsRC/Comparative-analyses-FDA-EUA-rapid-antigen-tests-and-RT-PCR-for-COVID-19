# Comparative analyses of eighteen rapid antigen tests and RT-PCR for COVID-19 quarantine and surveillance-based isolation

Chad R. Wells, Abhishek Pandey, Seyed M. Moghadas, Burton H. Singer, Gary Krieger, Richard J.L. Herron, David E. Turner , Justin P. Abshire , Kimberly M. Phillips , A. Michael Donoghue, Alison P. Galvani and Jeffrey P. Townsend

Copyright (C) <2022>, Chad R. Wells et. al. All rights reserved. Released under the GNU General Public License (GPL)

This repository contains codes and data used to simulate and analyze COVID- testing strategies in the scenarios of 
1. Quarnatine
2. Serial testing

The model code is written in MATLAB and results are saved as MATLAB data files (extension .mat), with plots also being constructed in MATLAB. 

## OS System requirements
The codes developed here are tested on Windows operating system (Windows 10 Home: 64-bit). 

## Installation guide
### MATLAB
Installation instruction for MATLAB can be found at https://www.mathworks.com/help/install/install-products.html. Typical install time for MATLAB on a "normal" desktop is around 30-40 minutes. The current codes were developed and tested on MATLAB R2019b.

## Demo
Figure1 generates Figure 1 in the main text. Corresponding figures and analysis can be generated by the provided scripts. All mat files to run this script are available. 

## Instructions for use
To generate the Figures and output of the calculations, select a script from Figures section to run in MATLAB and enter the name in the command line. All mat file are available to generate figures and conduct the calculations. To run analysis on a different set of parameters, adjust the parameters in the script and enter the name of the script in the command line to run. 

There is a separate folder for each of the different scenarios considered in the manuscript:
Baseline: Delta_Variant
5.72 day incubation period: Non_Delta
Alternative RT-PCR curve: Alternative_Curve_Delta_Variant
Diagnostic sensitivty cutoff for antigen test: Ag_Cut

In these folders there are scripts to run the analysis for quarantine and serial testing

DTX0 - RT-PCR testing on exit from quarantine with a 24-h delay
DTFR0 - RT-PCR serial testing with a 24-h delay
TEX1 - Rapid antigen testing on entry and exit from quarantine
TX1 - Rapid antigen testing on from quarantine
TFR1 - Rapid antigen serial testing
