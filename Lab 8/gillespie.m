function [xplot,tplot] = gillespie(c, x, tmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This function solves the Gillespie algorithm. It requires the rates and
% the initial states, along with a maximum time. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = length(c);
N= length(x); 

t = 0; %sets initial time 
tplot = zeros(1); 
k = 2; 
xplot = ones(1); 
tau=1; 

while (t<tmax) && tau~=0
    r1 = rand;
    r2 = rand; 
    a = c*x;  
    A = sum(a); 
    tau = -log(r1)/A;
    t = t+tau; 
    X = X+a; 
    tplot(k) = t; 
    xplot(k) = X; 
    k = k+1; 
end
end