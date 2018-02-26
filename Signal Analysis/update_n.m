function [n]=update_n(mv,v_rest, dt,n);
% Function that calculates the n coefficient
% in the H&H equations based on trans-membrane potential mv
% and resting potential v_rest
% Dimensions:   potentials are in mV
%               dn is a dimensionless coefficient update
%               dt time interval in ms

v=-mv-v_rest;
	bn=.125*exp(v/80);
	% avoid division by 0
	if v==-10; v=10.1;end
	an=.01*(v+10)/(exp((v+10)/10)-1);
	% Updates for the coefficients n
	% ------------------------------
	dn=(an*(1-n)-bn*n)*dt; 
    n=n+dn;
	
   