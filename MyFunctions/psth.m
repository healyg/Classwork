function psth(Spikes, binsize, tlength)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function creates a peristimulus time histogram of spike input data.
% 
% INPUTS: 
% 
% Spikes - the spike data. Can either be in cell or binary form.
%
% color - allows user control over plot color (black is the most
% professional, but orange is pretty damn cool. 
% 
% binsize - the size of the bins the data will be fit into (in seconds). 
%
% tlength - the length of time the data exists over
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

hold all

numbins = (tlength/binsize);
histarray = [];
quickarray = [];

if iscell(Spikes)
    s = size(Spikes);
    for ii = 1:s(2)
        trial = Spikes{:,ii};
        for jj = 1:length(trial)
            spkx = trial(jj);
            quickarray(jj) = spkx;
        end
        histarray = [histarray,quickarray];
        quickarray = [];
    end    
    h = histogram(histarray,numbins)
else 
    disp('wtf')
end 

xlabel('Time(s)');ylabel('Counts');
title('Peri-stimulus Time Histogram');
end