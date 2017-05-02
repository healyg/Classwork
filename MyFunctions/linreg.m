function [a,b,r] = linreg(explanatory, response) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This function plots a linear regression
% 
% INPUTS:
% explanatory - the explanatory (independent) variable
% 
% response - the response (dependent) variable
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nr = length(response); 
ne = length(explanatory);
rm = mean(response);
em = mean(explanatory);
stdr = std(response);
stde = std(explanatory);
xmin = min(explanatory);
xmax = max(explanatory);
ymin = min(response);
ymax = max(response);

if ne == nr
    varstats = cov(explanatory, response); %covariance matrix
    covar = varstats(1,2); %covariance(x,y)
    varx = varstats(1,1); %variance of x alone
%   vare = var(explanatory); %explanatory variance
%   varr = var(response); %response variance
    b = covar/(varx); %slope
    a = rm - b*(em);  %initial value
    sve = (explanatory - em)/stde; %standardized value 
    svr = (response - rm)/stdr; %standardized value
    r = (sum(sve.*svr) / (nr-1)); %correlation coefficient
    
    t = 0:ne; 
    reg = a + b*t; %regression formula

    figure 
    hold all; 
    scatter(explanatory, response);
    plot(t,reg);
    xlabel('Explanatory');ylabel('Response');
    axis([(xmin-1) (xmax+1) (ymin-1) (ymax+1)]);
else 
    disp('Input variables need to be the same size')
end 

end