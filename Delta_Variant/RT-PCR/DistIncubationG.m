function [f,F] = DistIncubationG(Inc)
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8392962/
% lognx=fmincon(@(x)mean(abs(X(:,2)-lognpdf(X(:,1),x(1),x(2)))),[1 1],[],[],[],[],[ 0 0],[5 5],[])
% If done through least squares, changes slightly; but does not have a
% drastic impact on the parameter estimation of the RT-PCR curve. 

% log-normal
%  p=[1.540253105635185   0.434826206700732];
%  
%  f=lognpdf(Inc,p(1),p(2));
%  F=logncdf(Inc,p(1),p(2));

gx=[5.70477179566407,0.852822456014538];

f=gampdf(Inc,gx(1),gx(2));
F=gamcdf(Inc,gx(1),gx(2));
end

