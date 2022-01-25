function [S,AG_TP,RTPCR_RTP,RTPCR_TP,TP] = Test_Specificity (CalltestName,MLEv)

% Returns the test specificity for the different tests

    % Specificity for the RT-PCR test
   RTPCR_TP=12392;
   TP=12404; 
    if(strcmp(CalltestName,'LumiraDX (Anterior Nasal Swab)'))
        AG_TP=168; 
        RTPCR_RTP=174;
    elseif(strcmp(CalltestName,'LumiraDX (Nasopharyngeal Swabs)'))
        AG_TP=210; 
        RTPCR_RTP=215;
    elseif(strcmp(CalltestName,'BD Veritor'))
        AG_TP=212; 
        RTPCR_RTP=213;
    elseif(strcmp(CalltestName,'Celltrion DiaTrust'))
        AG_TP=102; 
        RTPCR_RTP=103;
    elseif(strcmp(CalltestName,'CareStart (Anterior Nasal Swab - FDA)'))
        AG_TP=53; 
        RTPCR_RTP=53;
    elseif(strcmp(CalltestName,'CareStart (Anterior Nasal Swab - External)'))
        AG_TP=1243; 
        RTPCR_RTP=1264;
    elseif(strcmp(CalltestName,'CareStart (Anterior Nasal Swab)'))
        AG_TP=1296; 
        RTPCR_RTP=1317;
    elseif(strcmp(CalltestName,'CareStart (Nasopharyngeal Swab)')) 
        AG_TP=147; 
        RTPCR_RTP=148;
    elseif(strcmp(CalltestName,'Clip COVID')) 
        AG_TP=134; 
        RTPCR_RTP=134;
    elseif(strcmp(CalltestName,'Liaison (Anterior Nasal Swab)'))
        AG_TP=108; 
        RTPCR_RTP=108;
    elseif(strcmp(CalltestName,'Liaison (Nasalpharyngeal Swab)')) 
        AG_TP=133; 
        RTPCR_RTP=134;
    elseif(strcmp(CalltestName,'Omnia')) 
        AG_TP=32; 
        RTPCR_RTP=32;
    elseif(strcmp(CalltestName,'SCoV-2')) 
        AG_TP=257; 
        RTPCR_RTP=257;
    elseif(strcmp(CalltestName,'Status COVID+Flu')) 
        AG_TP=76; 
        RTPCR_RTP=76;
    elseif(strcmp(CalltestName,'Vitros')) 
        AG_TP=75; 
        RTPCR_RTP=75;
    elseif(strcmp(CalltestName,'Sofia 2 Flu+SARS')) 
        AG_TP=122; 
        RTPCR_RTP=122;
    elseif(strcmp(CalltestName,'Simoa')) 
        AG_TP=38; 
        RTPCR_RTP=38;
    elseif(strcmp(CalltestName,'Ellume')) 
        AG_TP=156; 
        RTPCR_RTP=161;
    elseif(strcmp(CalltestName,'Sofia (FDA)')) 
        AG_TP=179; 
        RTPCR_RTP=179;
    elseif(strcmp(CalltestName,'Sofia (CDC)')) 
        AG_TP=1025; 
        RTPCR_RTP=1041;
    elseif(strcmp(CalltestName,'Sofia')) 
        AG_TP=1204;
        RTPCR_RTP=1220;
    elseif(strcmp(CalltestName,'BinaxNOW (FDA)'))
        AG_TP=338; 
        RTPCR_RTP=343;
    elseif(strcmp(CalltestName,'BinaxNOW (Community)'))
        AG_TP=2004;
        RTPCR_RTP=2016;
    elseif(strcmp(CalltestName,'BinaxNOW'))
        AG_TP=2342; 
        RTPCR_RTP=2359;
    else        
        AG_TP=[];
        RTPCR_RTP=[];
    end

if(MLEv~=0)
    S=RTPCR_TP/TP;
    if(~isempty(AG_TP))
        S=S.*(AG_TP/RTPCR_RTP);
    end
else
    p=linspace(0,1,10001);
    L=binopdf(RTPCR_TP,TP,p);
    w=cumsum(L./sum(L));
    r=rand(1);
    findx=find(r<=w,1);
    S=p(findx);
    if(~isempty(AG_TP))
        L=binopdf(AG_TP,RTPCR_RTP,p);
        w=cumsum(L./sum(L));
        r=rand(1);
        findx=find(r<=w,1);
        S=S.*p(findx);
    end
end

end