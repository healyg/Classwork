%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script is used to answer problem 5: correlation on the midterm 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all 
close all
load('MEA_Data.mat');
c1 = MEA_Data(1,:);
c2 = MEA_Data(2,:);
s = size(MEA_Data);
i_prev = -20;
threshold=0.6*max(c1);
n = 1;

for i=2:length(c1)-1
    if (c1(i)>threshold)                                           
      if((c1(i)>=c1(i-1))&&(c1(i)>=c1(i+1)))      
          if (i-i_prev > 40)                                               
              spikes1(n)=i;
              i_prev=i;
              n=n+1;
          end
      end
    end
end

i_prev = -20;
threshold=0.4*max(c2);
n = 1;

for i=2:length(c2)-1
    if (c2(i)>threshold)                                           
      if((c2(i)>=c2(i-1))&&(c2(i)>=c2(i+1)))       
          if (i-i_prev > 40)                                               
              spikes2(n)=i;
              i_prev=i;
              n=n+1;
          end
      end
    end
end

figure
hold on
plot(c1);
for i = 1:length(spikes1)
    spkx = [spikes1(i), spikes1(i)];
    spky = [96,104];
    line(spkx,spky,'color','k','LineWidth',0.5);
end
xlabel('Time');ylabel('Voltage muV (and non-unit raster)');title('CH1 Electrode Recording and Raster Spike Train');
hold off

figure
hold on
plot(c2);
for i = 1:length(spikes2)
    spkx = [spikes2(i), spikes2(i)];
    spky = [96,104];
    line(spkx,spky,'color','k','LineWidth',0.5);
end
xlabel('Time');ylabel('Voltage muV (and non-unit raster)');title('CH2 Electrode Recording and Raster Spike Train');
hold off

% correlations 

[a1,lags1]=xcorr(spikes1,spikes1);
[a2,lags2]=xcorr(spikes2,spikes2);
[cross,xlags]=xcorr(spikes1,spikes2);

figure 
plot(lags1,a1)
title('Autocorrelation CH1');xlabel('Lag');ylabel('Correlation');

figure
plot(lags2,a2)
title('Autocorrelation CH2');xlabel('Lag');ylabel('Correlation');

figure 
plot(xlags,cross)
title('Cross Correlation CH1&2');xlabel('Lag');ylabel('Correlation');

% spike triggered average

ws = 400;
spikes2_sta = spikes2(spikes2>=ws);
n = length(spikes2_sta);
avg = zeros(1,ws);
savg = zeros(1,ws);

for i = 1:n
    w = spikes2_sta(i)-ws:spikes2_sta(i)-1;
    loc = find(spikes1>w(1)&spikes1<w(ws));
    avg = avg + c1(w);
    if loc ~=0
        for j = 1:length(loc)
            ind = loc(j);
            time = spikes1(ind);
            t = spikes2_sta(i)-time;
            savg(t) = savg(t) + 1;
        end
    end
end

avg = avg./n;
x = -400:-1;
figure 
hold on 
plot(x,avg);
title('Spike Triggered Average, CH1 Stimulus, CH2 Response');
xlabel('Time Before Spike');ylabel('Average Trigger Voltage (muV)');
hold off

figure 
hold on 
plot(x,savg);
title('Spike Triggered Average, CH1 Stimulus, CH2 Response');
xlabel('Time Before Spike');ylabel('Average # of Spikes');
hold off

%% Problem 6 

clear all
close all

t = -0.5:0.1:5;
y = (9*exp(-2*t))-((17/2)*exp(-3*t))+((1/2)*exp(-t));

figure
plot(t,y)
xlabel('Time');ylabel('Signal Transform');title('Laplace Transform of ODE');