function S = Test_Specificity (CalltestName)

% Returns the test specificity for the different tests


% Specificity for the RT-PCR test
S=12392/12404; 

if(strcmp(CalltestName,'LumiraDX (Anterior Nasal Swab)'))
    S= S*(168/174);
elseif(strcmp(CalltestName,'LumiraDX (Nasopharyngeal Swabs)'))
    S= S*(210/215);
elseif(strcmp(CalltestName,'BD Veritor'))
    S= S*(212/213);
elseif(strcmp(CalltestName,'Celltrion DiaTrust'))
    S= S*(102/103);
elseif(strcmp(CalltestName,'CareStart (Anterior Nasal Swab - FDA)'))
    S= S*(53/53);
elseif(strcmp(CalltestName,'CareStart (Anterior Nasal Swab - External)'))
    S= S*(1243/1264);
elseif(strcmp(CalltestName,'CareStart (Anterior Nasal Swab)'))
    S= S*(1296/1317);
elseif(strcmp(CalltestName,'CareStart (Nasopharyngeal Swab)')) 
    S= S*(147/148);
elseif(strcmp(CalltestName,'Clip COVID')) 
    S= S*(134/134);
elseif(strcmp(CalltestName,'Liaison (Anterior Nasal Swab)'))
    S= S*(108/108);
elseif(strcmp(CalltestName,'Liaison (Nasalpharyngeal Swab)')) 
    S= S*(133/134);
elseif(strcmp(CalltestName,'Omnia')) 
    S= S*(32/32);
elseif(strcmp(CalltestName,'SCoV-2')) 
    S= S*(257/257);
elseif(strcmp(CalltestName,'Status COVID+Flu')) 
    S= S*(76/76);
elseif(strcmp(CalltestName,'Vitros')) 
    S= S*(75/75);
elseif(strcmp(CalltestName,'Sofia 2 Flu+SARS')) 
    S= S*(122/122);
elseif(strcmp(CalltestName,'Simoa')) 
    S= S*(38/38);
elseif(strcmp(CalltestName,'Ellume')) 
    S= S*(156/161);
elseif(strcmp(CalltestName,'Sofia (FDA)')) 
    S= S*(179/179);
elseif(strcmp(CalltestName,'Sofia (CDC)')) 
    S= S*(1025/1041);
elseif(strcmp(CalltestName,'Sofia')) 
    S= S*(1204/1220);
elseif(strcmp(CalltestName,'BinaxNOW (FDA)'))
    S= S*(338/343);
elseif(strcmp(CalltestName,'BinaxNOW (Community)'))
    S= S*(2004/2016);
elseif(strcmp(CalltestName,'BinaxNOW'))
    S= S*(2342/2359);
end
end