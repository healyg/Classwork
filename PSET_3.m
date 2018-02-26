%
%
%
clear all 
close all

pr1_1

figure 
hold all
sq = EEG.*EEG/length(EEG);
ps = (fft(EEG).*conj(fft(EEG)))/length(EEG);
as = (2/length(EEG))*sqrt(fft(EEG).*conj(fft(EEG)));
R = real(fft(EEG));
I = imag(fft(EEG));
phs = atan(I./R);
plot(ps)
hold off
figure 
plot(as)
figure 
plot(phs)

%%

clear all
close all

T=1;
t=0:.1:T;
N=length(t);
y=cos(2*pi*1*t);
figure;plot(t,y);
title('T=1 Time Series');
ps =(fft(y).*conj(fft(y)))/N;
figure
plot(ps)
title('T=1 Power Spec');

T=1.3;
t=0:.1:T;
N=length(t);
y=cos(2*pi*1*t);
figure
ps =(fft(y).*conj(fft(y)))/N;
plot(ps)
title('T=1.3 Power Spec');

T=1;
t=0:.1:T;
N=length(t);
y=cos(2*pi*1*t);
w = hann(N);figure;plot(w);
x = conv(y,w);
figure
ps =(fft(x).*conj(fft(x)))/N;
plot(ps)
title('Hann')

%%
clear all
close all

x = zeros(1,10000);
y1 = zeros(1,10020);
y2 = zeros(1,10020);
x(1) = randn;
x(2) = randn;
for n = 3:10000
    x(n) = randn;
    y1(n) = x(n) + x(n-1);
    y2(n) = 0.5*x(n) + x(n-1) + 3*x(n-2);
end

mx = mean(x);
my1 = mean(y1);
my2 = mean(y2);

[c,l] = xcov(y1,y2);

plot(l,c);axis([-6 6 0 max(c)]);

c2 = zeros(1,13);
for m = -6:6
    for n = 9:10000
        y2(n) = 0.5*x(n+m) + x(n-1+m) + 3*x(n-2+m);
        c2(m+7) = mean(y1.*y2);
    end
end

plot(l,c2);axis([-6 6 0 max(c)]);