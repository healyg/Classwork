clear all
close all

[el1 el2 el3]=Simulated_Signal(3,0.01,0.3);
save test el1 el2 el3;
pr14_2

% el1 and el3 are much more highly coupled than el1 and el2. When looking
% at the signals themselves, it is hard to see any real differnce. However,
% looking at their spectrum shows the similarities between el1 & el3, as
% well as the differences in el2. el1 & el3 have nearly indentical spectra,
% whereas el2 has a vastly different one. 

% correlations 

[c12,l12]=xcorr(el1,el2,'coeff'); 
[c23,l23]=xcorr(el2,el3,'coeff');
[c13,l13]=xcorr(el1,el3,'coeff');

figure 
plot(l12,c12);
title('Correlation against Lag for el1 & el2');

figure 
plot(l23,c23);
title('Correlation against Lag for el2 & el3');

figure 
plot(l13,c13);
title('Correlation against Lag for el1 & el3');

% By looking at these graphs, we can see that the xcorr is much more
% regularly structured and centered around 0 lag for the relationship
% between el1 & el3. While there is obviously a relationship between the
% less coupled el1 & el2, it is not centered extremely well around 0 and is
% not symmetrical. 
% Solely from these graphs, I can establish Granger causality between el1 &
% el3, since they exist with 0 lag, I would be able to know more about el3
% by knowing el1. The same could be said about el1 & el2, but to a far
% lesser degree. 


