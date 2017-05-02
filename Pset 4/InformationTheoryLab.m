%% PSET #4 
% The file mtNeuron.mat contains the responses of a single directionally 
% tuned MT neuron to random dot stimuli moving coherently in directions 
% varying from -90 to 90 degrees relative to the previously measured 
% preferred direction of the neuron. The thirteen stimuli were each 
% presented 184 times. Each stimulus began at time 0 and continued for 256 
% ms. Recordings continued until 512 ms after the beginning of the stimulus. 
% The array has the dimensions 256×13×184. The first dimension is time in 2 ms 
% bins, the second dimension motion direction, the third dimension the 
% repeated presentations.
% 
%% Part 1 
% 
% Create a raster of all the spike trains in one plot, sorted by stimulus 
% direction. Rasters corresponding to different directions should be 
% represented in different colors.
% 

close all; 

load('mtNeuron.mat')

data = mtNeuron.data;
time = mtNeuron.time;
dirs = mtNeuron.dirs; 
s=size(data);

figure 
hold on 
xlabel('Time(ms)');ylabel('Direction');
title('Neural Spikes by Stimulus Direction');
for i = 1:s(1)
    for j= 1:s(2)
        for k = 1:s(3)
            if data(i,j,k) == 1
                if j == 1
                    vcolor = [1 0 0]; %red
                elseif j == 2
                    vcolor = [1 0.5 0]; %orange
                elseif j == 3
                    vcolor = [1 1 0]; %yellow
                elseif j == 4
                    vcolor = [0 1 0]; %green
                elseif j == 5
                    vcolor = [0 0 1]; %blue
                elseif j == 6
                    vcolor = [0.5 0 1]; %indigo
                elseif j == 7
                    vcolor = [1 0 1]; %magenta
                elseif j == 8
                    vcolor = [0 1 1]; %teal
                elseif j == 9
                    vcolor = [0.5 0 0]; %dim red
                elseif j == 10
                    vcolor = [0.5 0.25 0]; %dim orange
                elseif j == 11
                    vcolor = [0.5 0.5 0]; %dim yellow
                elseif j == 12
                    vcolor = [0 0.5 0]; %dim green
                elseif j == 13
                    vcolor = [0 0 0.5]; %dim blue
                else
                end
                spkx = [i*2, i*2];
                spky = [dirs(j) - 7.4,dirs(j) + 7.4];
                line(spkx,spky,'color',vcolor,'LineWidth',1);
            else 
            end
        end
    end
end
hold off

%% Part 2 
% 
% Compute and plot the mutual information between cumulative spike count 
% and motion direction as a function of time. If you’d like to check your 
% work, this neuron was published in Osborne et al 2004, Figure 3C. 
% 


imax =-1*13*(1/13)*log2(1/13);
vals = zeros(s(1),s(2)); 
trials = zeros(s(1),s(3));
spikes = 0:23;
d = cumsum(data);
tots = zeros(1,s(2));
newinf = zeros(s(1),1);
stdvals = zeros(s(1),1);

for q = 1:s(1)
    counts = zeros(length(spikes),s(2));
    unc = zeros(1,s(2));
    for i = 1:s(2)
        for j = 1:s(3)
           x = d(q,i,j); 
           counts(x+1,i) = counts(x+1,i) + 1;
        end
    end
    for i = 1:s(2)
        for j = 1:length(spikes)
            tot2 = sum(counts(j,:));
            if tot2 ~=0
                prob = counts(j,i)/ 184;
                mprob = mean(counts(j,:))/184;
                if prob~=0
                    unc(i) = unc(i) - prob*log2(prob/mprob);
                end
            end
        end
    end
    newinf(q) = - mean(unc);
    stdvals(q) = std(unc);
end

figure
hold on
%plot(time*1000,newinf)
errorbar(time*1000,newinf,stdvals)
xlabel('Time(ms)');ylabel('Information(bits)');
title('Information Contained in Spike Counts with Regard to Stimulus Direction');
hold off


%% Part 3 
% 
% Determine the latency of this neuron in a principled way. What proportion 
% of the mutual information is available within the first 50 ms of the 
% neural response? Within the first 100 ms?
% 

r = zeros(s(1),s(2));
for i = 1:s(1)
    for j = 1:s(2)
        for k = 1:s(3)
            d = data(i,j,k);
            if d == 1
                r(i,j) = r(i,j) + 1;
            end
        end
    end
end

r = r/184;

m = max(r);
lat = zeros(s(2),1); %Latency values by direction

for i = 1:s(2)
    v=0;
    for j = 1:s(1)
        if r(j,i)>=(m(i)*0.15) && v==0
            lat(i) = time(j)*1000;
            v=1;
        end
    end
end

meanlat = mean(lat);

disp('The latency is about ' + string(meanlat) + ' ms.')

first50 = (newinf(65)/max(newinf))*100;
first100 = (newinf(90)/max(newinf))*100;

disp('About ' + string(first50) + '% of the information is available in the' + ...
    ' first 50 seconds of neural responds, and about ' + string(first100) + ...
    '% of the information is available in the first 100 seconds.');


