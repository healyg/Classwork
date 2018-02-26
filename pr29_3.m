% pr29_3.m
% Integrate-and-Fire Neuron
% cascade of a leaky-lowpass filter (~integrate) 
% and a threshold (fire)
% Linear component is a leaky low-pass RC circuit
% we use R=10k and C=3.3uF
% similar to the filter in CH 15 and CH 16

clear;
close all;
% Parameters
R=10e3;
C=3.3e-6;
RC=R*C;
sample_rate=1000; 
dt=1/sample_rate;
time=1;
A=RC/dt;
threshold=20;
AP=100;
y_previous=0;                    % initial value

% Input
Amp=20.1;                        % Amplitude of the input
N=floor(time*sample_rate);       % Size of the input
x=ones(1,N)*Amp;                 % the input is a step; initial part is 0 
x(1:round(N/5))=zeros(1,round(N/5));   
% OR ALTERNATIVELY uncomment the following line and input is random noise 
%Amp=150; x=randn(1,N)*Amp;      
 
    for n=1:N;
        y(n)=(A*y_previous+x(n))/(A+1); % the linear component; see CH 16
        if (y(n) < threshold);          % the nonlinear part
            y_previous=y(n);
            zz(n)=y(n);                 % subthreshold behavior
        else;
            y_previous=0;
            zz(n)=AP;                   % action potential
        end;
    end;
zz=zz-65;                         % Correct Membrane potential 
                                  % for rest potential at -65 mV

figure;
subplot(211),plot(x);
title('input')
ylabel('Amplitude')
subplot(212);plot(zz,'r');
title('output')
xlabel('Time (ms)'),ylabel('Amplitude')
