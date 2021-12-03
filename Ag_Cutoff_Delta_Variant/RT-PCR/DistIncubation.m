function [f,F] = DistIncubation(Inc)
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8392962/

 p=[1.540253105635185   0.434826206700732];
 
 f=lognpdf(Inc,p(1),p(2));
 F=logncdf(Inc,p(1),p(2));

end

