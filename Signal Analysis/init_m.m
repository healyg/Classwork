function [m,mt]=init_m(mv,v_rest);
% Function that calculates the  m  coefficient
% in the H&H equations based on trans-membrane potential mv
% and resting potential v_rest
% Dimensions:   potentials are in mV
%               m is a dimensionless coefficient
%               mt is the time constant

v =-mv-v_rest;  % calculate departure from rest
	bm=4*exp(v/18);
                % avoid division by 0
	if v==-23; v=-23.01;end;
	am=.1*(v+23)/(exp((v+23)/10)-1);	
	% Coefficients
	% ------------	
	m=am/(am+bm);
    mt=1/(am+bm);
	