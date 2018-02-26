function [m]=update_m(mv,v_rest, dt, m);
% Function that calculates the  m  coefficient
% in the H&H equations based on trans-membrane potential mv
% and resting potential v_rest
% Dimensions:   potentials are in mV
%               dmis a dimensionless coefficient update
%               dt time interval in ms

v=-mv-v_rest;
	bm=4*exp(v/18);
	% avoid division by 0
	if v==-23; v=23.1;end
	am=.1*(v+23)/(exp((v+23)/10)-1);
	% Updates for the coefficient m
	% -----------------------------
	dm=(am*(1-m)-bm*m)*dt; 
    m=m+dm;

	
   