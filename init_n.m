function [n,nt]=init_n(mv,v_rest);
% Function that calculates the n coefficient
% in the H&H equations based on trans-membrane potential mv
% and resting potential v_rest
% Dimensions:   potentials are in mV
%               n is a dimensionless coefficient
%               nt is the time constant

v =-mv-v_rest;  % calculate departure from rest
	bn=.125*exp(v/80);
                % avoid division by 0
	if v==-10; v=-10.1;end
 	an=.01*(v+10)/(exp((v+10)/10)-1);
	% Coefficients
	% ------------
	n=an/(an+bn);
    nt=1/(an+bn);
	