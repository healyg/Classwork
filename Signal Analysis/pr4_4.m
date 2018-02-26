% pr4_4.m
% detection of QRS complexes in neonatal ECG
% this program requires the ascii file subecg

% initialization and vars
clear;
sr=1000;                                                % sample rate
fN=sr/2;                                                % Nyquist
level=0.8;                                              % level for the detection threshold
i_prev=-50;
n=1;

% 1. load data
load('MEA_Data.mat');

% 2. pre-process the data 
[c,d]=butter(2,[15/fN 45/fN],'bandpass');               % 2nd order 15-45 Hz

meaFF=filtfilt(c,d,MEA_Data-mean(MEA_Data));             % use filtfilt to prevent phase-shift
                            % Artifact rejection could be included here
                            
% 3. detect peaks above threshold in D
                            % In this routine we only look into the nearest neighbours (i.e. three subsequent points)
                            % adding additional points will make the algorithm more robust
threshold=level*max(meaFF);                          % detection threshold 
for i=2:length(meaFF)-1
    if (meaFF(i)>threshold)                                           % first check if the level is worthwhile analyzing
      if((meaFF(i)>=meaFF(i-1))&&(meaFF(i)>=meaFF(i+1)))       % is the point a relative maximum (Note the >= sign)?
          if (i-i_prev > 50)                                                % if yes, is it not a subsequent maximum in the same heartbeat
              D(n)=i;
              i_prev=i;
              n=n+1;
          end;
      end;
    end;
end;

% 4. plot the result
H=ones(length(D),1)*max(meaFF);
figure;
hold;
plot(MEA_Data-mean(MEA_Data));
plot(D,H,'r.');
title('ECG and detected QRS complexes (USE ZOOM TO INSPECT THE DETECTION)');
xlabel('Time (ms)');
ylabel('Amplitude (AU)');
