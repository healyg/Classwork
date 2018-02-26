

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
xlabel('Time');ylabel('Voltage muV (and non-unit raster)');title('Electrode Recording and Raster Spike Train');
hold off

figure
hold on
plot(c2);
for i = 1:length(spikes2)
    spkx = [spikes2(i), spikes2(i)];
    spky = [96,104];
    line(spkx,spky,'color','k','LineWidth',0.5);
end
xlabel('Time');ylabel('Voltage muV (and non-unit raster)');title('Electrode Recording and Raster Spike Train');
hold off

% x-correlations 

[a1,lags1]=xcorr(spikes1,spikes1);
[a2,lags2]=xcorr(spikes2,spikes2);
[cross,xlags]=xcorr(spikes1,spikes2);

figure 
plot(a1)

figure
plot(a2)

figure 
plot(cross)

% spike triggered average

ws = 400;
spikes2_sta = spikes2(spikes2>=ws);
n = length(spikes2_sta);
avg = zeros(1,ws);

for i = 1:n
    w = spikes2_sta(i)-ws:spikes2_sta(i)-1;
    avg = avg + c1(w);
end

avg = avg./n;

figure 
hold on 
plot(avg);
hold off






