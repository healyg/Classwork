function [h,ht]=init_h(mv,v_rest);
% Function that calculates the  h coefficient
% in the H&H equations based on trans-membrane potential mv
% and resting potential v_rest
% Dimensions:   potentials are in mV
%               h is a dimensionless coefficient
%               ht is the time constant


v =-mv-v_rest;  % calculate departure from rest
                % avoid division by 0
    if v==-35; v=-35.01;end;
	bh=1/(exp((v+35)/10)+1);
	ah=.07*exp(v/20);
	% Coefficients
	% ------------
	h=ah/(ah+bh);
    ht=1/(ah+bh);