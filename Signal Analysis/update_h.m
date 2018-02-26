function [h]=update_h(mv,v_rest, dt, h);
% Function that calculates the  h coefficient
% in the H&H equations based on trans-membrane potential mv
% and resting potential v_rest
% Dimensions:   potentials are in mV
%               dh is the dimensionless coefficient updates
%               dt time interval in ms

v=-mv-v_rest;
	bh=1/(exp((v+35)/10)+1);
	ah=.07*exp(v/20);
	% Updates for the coefficients h
	% -------------------------------
	dh=(ah*(1-h)-bh*h)*dt; 
    h=h+dh;
	
   