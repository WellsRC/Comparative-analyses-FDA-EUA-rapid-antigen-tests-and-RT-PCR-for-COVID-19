function [Name] = AdjustedNames_Plotting(testName)

if(strcmp('LumiraDX (Anterior Nasal Swab)',testName))
    Name='LumiraDX (AN)';
elseif(strcmp('CareStart (Anterior Nasal Swab - FDA)',testName))    
    Name='CareStart (AN)';
elseif(strcmp('BinaxNOW (FDA)',testName))  
    Name='BinaxNOW';
elseif(strcmp('Liaison (Anterior Nasal Swab)',testName))
    Name='Liaison (AN)';
elseif(strcmp('Sofia (FDA)',testName))
    Name='Sofia';
elseif(strcmp('Liaison (Nasalpharyngeal Swab)',testName))
    Name='Liaison (NS)';
elseif(strcmp('LumiraDX (Nasopharyngeal Swabs)',testName))
    Name='LumiraDX (NS)';
elseif(strcmp('CareStart (Nasopharyngeal Swab)',testName))
    Name='CareStart (NS)';
elseif(strcmp('SCoV-2',testName))
    Name='SCoV-2 Ag Detect';
elseif(strcmp('Status COVID+Flu',testName))
    Name='Status COVID-19/Flu';    
else
    Name=testName;
end

