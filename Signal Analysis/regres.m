function [m,b,r]=regres(x,y);
% Calculates Regression Line 
% and correlation based on 
% input variables
% x and y

% Calculate Regression Help variables
n=length(x);
Sx=sum(x);
Sy=sum(y);
Sxy=sum(x.*y);
Sx2=sum(x.^2);
Sy2=sum(y.^2);
% Direction (m), Intercept (b) and Correlation (r)
m=(n*Sxy-Sx*Sy)/(n*Sx2-Sx^2);
b=(Sy-m*Sx)/n;
r=(n*Sxy-Sx*Sy)/(sqrt((n*Sx2-Sx^2)*(n*Sy2-Sy^2)));
R2=r^2;


