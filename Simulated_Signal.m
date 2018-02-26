function [S1 S2 S3] = Simulated_Signal(SNR, K2, K3);
% Output signals S1, S2, S3 
% the signal-to-noise for S1 is determined by SNR.
% In the example in the text we use SNR = 3. 
% The coupling from S1 to S2 and S3 are determined by K2 and K3
% In the text we employ 0.1 and 0.001 for K2 and K3 resp.
% Delay between S2 and S1 is dly2 and between S3 and S1 is dly3

%parameters
sample_rate=400;                        % sr in Hz
freq=30;                                % f in Hz
tim=40;                                 % time in seconds
         
dly2=5;              % # points delay and coupling strength signal 1 --> 2
dly3=10;             % # points delay and coupling strength signal 1 --> 3

% derived parms
f_Nyq=sample_rate/2;                    % Nyquist
dt=1/sample_rate;                       % time step
N=tim*sample_rate;                      % no of points
t=0:dt:tim;                             % time axis
noise=randn(1,length(t));               % noise component
null=zeros(1,length(t));                % basline of zeros

% Source Signal
S1=sin(2*pi*freq*t).*SNR*std(noise);
% Create the derived signals S2 and S3 by adding delay and noise
S2=null; S3=null;
for k=dly2+1:N
    S2(k)=S1(k-dly2).*K2+randn;         % coupling
end;
for k=dly3+1:N
    S3(k)=S1(k-dly3).*K3+randn;         % coupling
end;
% Noise added to S1
S1=S1+noise;    
% Normalize all signals
S1=(S1-mean(S1))/std(S1);               
S2=(S2-mean(S2))/std(S2);
S3=(S3-mean(S3))/std(S3);

