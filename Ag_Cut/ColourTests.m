function [C,MF]= ColourTests(CalltestName)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


if(strcmp(CalltestName,'RTPCR'))
    C=hex2rgb('#080706');
    MF={'p'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gradual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif(strcmp(CalltestName,'BinaxNOW'))
    C=hex2rgb('#810f7c');    
    MF={'s'};
elseif(strcmp(CalltestName,'BinaxNOW (FDA)'))
    C=[hex2rgb('4D85BD');hex2rgb('#810f7c')];    
    MF={'s'};
elseif(strcmp(CalltestName,'BinaxNOW (Community)'))
    C=hex2rgb('#CB0000');    
    MF={'o'};
elseif(strcmp(CalltestName,'Sofia (CDC)')) 
    C=hex2rgb('#ce1256');  
    MF={'d'};
elseif(strcmp(CalltestName,'Simoa')) 
    C=hex2rgb('#9ebcda');  
    MF={'o'};
elseif(strcmp(CalltestName,'Status COVID+Flu')) 
    C=hex2rgb('#8c96c6');  
    MF={'^'};
elseif(strcmp(CalltestName,'SCoV-2')) 
    C=hex2rgb('#8c6bb1');  
    MF={'>'};
elseif(strcmp(CalltestName,'CareStart (Anterior Nasal Swab - FDA)')) 
    C=[hex2rgb('#9ecae1'); hex2rgb('#88419d');];     
    MF={'h'};
elseif(strcmp(CalltestName,'LumiraDX (Anterior Nasal Swab)'))
    C=hex2rgb('#bfd3e6');   
    MF={'v'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif(strcmp(CalltestName,'Sofia')) 
    C=hex2rgb('#41ae76');  
    MF={'o'};
elseif(strcmp(CalltestName,'Sofia (FDA)')) 
    C=[hex2rgb('#4292c6');hex2rgb('#41ae76')];   
    MF={'s'};
elseif(strcmp(CalltestName,'Vitros')) 
    C=hex2rgb('#00441b');  
    MF={'d'};
elseif(strcmp(CalltestName,'Omnia')) 
    C=hex2rgb('#006d2c');  
    MF={'^'};
elseif(strcmp(CalltestName,'Liaison (Anterior Nasal Swab)')) 
    C=hex2rgb('#66c2a4');  
    MF={'>'};
elseif(strcmp(CalltestName,'Liaison (Nasalpharyngeal Swab)')) 
    C=hex2rgb('#238b45');  
    MF={'h'};
elseif(strcmp(CalltestName,'Clip COVID')) 
    C=hex2rgb('#99d8c9');  
    MF={'v'};
elseif(strcmp(CalltestName,'CareStart (Anterior Nasal Swab)')) 
    C=hex2rgb('#006d2c');      
    MF={'<'};
elseif(strcmp(CalltestName,'CareStart (Anterior Nasal Swab - External)')) 
    C=hex2rgb('#980043');      
    MF={'s'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rapid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif(strcmp(CalltestName,'Ellume')) 
    C=hex2rgb('#fe9929');  
    MF={'s'};
elseif(strcmp(CalltestName,'Sofia 2 Flu+SARS')) 
    C=hex2rgb('#ec7014');    
    MF={'d'};
elseif(strcmp(CalltestName,'CareStart (Nasopharyngeal Swab)')) 
    C=hex2rgb('#cc4c02');
    MF={'o'};
elseif(strcmp(CalltestName,'Celltrion DiaTrust')) 
    C=hex2rgb('#fec44f');  
    MF={'^'};   
elseif(strcmp(CalltestName,'BD Veritor'))
    C=hex2rgb('#993404');  
    MF={'>'};    
elseif(strcmp(CalltestName,'LumiraDX (Nasopharyngeal Swabs)'))
    C=hex2rgb('#fee391');   
    MF={'h'};
end
end