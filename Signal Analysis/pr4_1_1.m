% pr4_1
% averaging
clear all
close all

sz=256;
NOISE_TRIALS=randn(sz); % a [sz x sz] matrix filled with noise

SZ=1:sz;                        % Create signal with sinus 
SZ=SZ/(sz/2);                   % Divide the array SZ by sz/2
S=sin(2*pi*SZ);

for i=1:sz                    % create a noisy signal 
    NOISE_TRIALS(i,:) = NOISE_TRIALS(i,:) + S;
end

% signal is created, manipulation from here-on out

new = zeros(sz);
for i = 1:sz
    r = ceil(rand*(sz-1))+1;
    new(i,:) = [NOISE_TRIALS(i,r:sz) NOISE_TRIALS(i,1:(r-1))];
end

% with sinus noise

average=sum(NOISE_TRIALS)/sz;   % create the average
odd_average=sum(NOISE_TRIALS(1:2:sz,:))/(sz/2);
even_average=sum(NOISE_TRIALS(2:2:sz,:))/(sz/2);
noise_estimate=(odd_average-even_average)/2;

figure
hold on
plot(NOISE_TRIALS(1,:),'g')
plot(noise_estimate,'k')
plot(average,'r')
plot(S)
title('Average RED, Noise estimate BLACK; Single trial GREEN, Signal BLUE')
hold off

% with sinus noise removed by random interval variance

average=sum(new)/sz;   % create the average
odd_average=sum(new(1:2:sz,:))/(sz/2);
even_average=sum(new(2:2:sz,:))/(sz/2);
noise_estimate=(odd_average-even_average)/2;

figure
hold on
plot(NOISE_TRIALS(1,:),'g')
plot(new(1,:),'m')
plot(noise_estimate,'k')
plot(average,'r')
plot(S)
title('Average RED, Noise estimate BLACK; Single trial GREEN, Single offset trial PURPLE, Signal BLUE')
hold off

% part c, plotting the ms of the noise versus n
ms = zeros(1,sz);
for i = 1:sz
    ms(i) = sum(noise_estimate(1,1:i).^2)/i;
end

x = 1:sz;
figure 
hold on 
plot(x,ms)
xlabel('Sample Size (n)');ylabel('Mean Squared Value');title('MS Changing with Sample Size');
hold off
