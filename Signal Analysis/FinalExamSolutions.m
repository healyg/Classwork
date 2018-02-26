%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Final Exam, Signal Analysis Spring 2017
% 
% Garrett Healy, 21MAY2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Question 2: Filter 
clear all 
close all 

fr = [0.3 0.5 1 2 5 10 20 50 100 200 500 1000];
amprat = [0.05 0.11 0.18 0.35 0.69 0.92 1 1 1 1 1 1];
atten = [-25.38 -19.36 -14.67 -9.21 -3.19 -0.70 0 0 0 0 0 0];

figure 
hold on
semilogx(fr,atten,'-o');

R = 10^4; %kAmps
C = 3.3e-6; %muF 
f = 0:1:1000;
y = 1-exp(-f*R*C); %solution in the time domain for ease of plotting

semilogx(f,y);legend('Measured Data','Theoretical Response');xlabel('Frequency');ylabel('Attenuation');

% the figures look very similar, both grow logarithmically and reach an
% asymptote at 1 after reaching high values. When plotted together,
% however, it can be seen that the theoretical representation does not
% start as low for the low frequencies as does the measured data. 
%% Question 3
clear all
close all

load handel; 
%p = audioplayer(y,Fs);
%play(p,[1 (get(p,'SampleRate')*8)]);
t=1/Fs:1/Fs:length(y)/Fs;

% filters and Bode plots 
[b1,a1] = butter(2,[50 100]/4096);
f1 = filtfilt(b1,a1,y);
figure 
freqz(b1,a1);title('50-100 Hz');
[b2,a2] = butter(2,[100 1000]/4096);
f2 = filtfilt(b2,a2,y);
figure 
freqz(b2,a2);title('100-1000 Hz');
[b3,a3] = butter(2,[1000 2000]/4096);
f3 = filtfilt(b3,a3,y);
figure 
freqz(b3,a3);title('1000-2000 Hz');
[b4,a4] = butter(2,[2000 4000]/4096);
f4 = filtfilt(b4,a4,y);
figure 
freqz(b4,a4);title('2000-4000 Hz');

% signal and bandpass filters on one plot 

figure 
hold on
plot(t,y)
plot(t,f1)
plot(t,f2)
plot(t,f3)
plot(t,f4)
title('Signal and Bandpass Filters');
legend('Original Signal','50-100 Hz','100-1000 Hz','1000-2000 Hz','2000-400 Hz');

%sounds 

p = audioplayer(f1,Fs);
play(p,[1 (get(p,'SampleRate')*8)]);
pause(4);
p = audioplayer(f2,Fs);
play(p,[1 (get(p,'SampleRate')*8)]);
pause(4);
p = audioplayer(f3,Fs);
play(p,[1 (get(p,'SampleRate')*8)]);
pause(4);
p = audioplayer(f4,Fs);
play(p,[1 (get(p,'SampleRate')*8)]);

% f1 was just deep, oscillating sounds.
% f2 sounded like a muffled version of the original 
% f3 was similar to f2, just more high pitched. Sounded like it was coming
% through a radio 
% f4 was similar to f3, just softer and higher pitched. 

% envelope 

r1 = abs(f1);r2 = abs(f2);r3 = abs(f3);r4 = abs(f4);

[be1,ae1] = butter(1,10/4096);
e1 = filtfilt(be1,ae1,r1);
[be2,ae2] = butter(1,10/4096);
e2 = filtfilt(be2,ae2,r2);
[be3,ae3] = butter(1,10/4096);
e3 = filtfilt(be3,ae3,r3);
[be4,ae4] = butter(1,10/4096);
e4 = filtfilt(be4,ae4,r4);

% pulse train

a = 1; %1V bipolar pulse
dur = [0 0.0005 0.0005 0.001]; %1ms duration 

pulse = zeros(length(f1),1);
c = 1;
while c<length(f1)
    pulse(c:c+7) = [1;1;1;1;-1;-1;-1;-1];
    c = c+8;
end

bp1 = pulse.*e1;
bp2 = pulse.*e2;
bp3 = pulse.*e3;
bp4 = pulse.*e4;
figure
subplot(2,2,1);plot(t,y);title('Original Signal');
subplot(2,2,2);plot(t,f1,'r');title('Bandpass Filter');
subplot(2,2,3);plot(t,e1,'k');title('Envelope');
subplot(2,2,4);plot(t,bp1,'g');title('Biphasic Pulse');
suptitle('First Filter');
figure
subplot(2,2,1);plot(t,y);title('Original Signal');
subplot(2,2,2);plot(t,f2,'r');title('Bandpass Filter');
subplot(2,2,3);plot(t,e2,'k');title('Envelope');
subplot(2,2,4);plot(t,bp2,'g');title('Biphasic Pulse');
suptitle('Second Filter');
figure
subplot(2,2,1);plot(t,y);title('Original Signal');
subplot(2,2,2);plot(t,f3,'r');title('Bandpass Filter');
subplot(2,2,3);plot(t,e3,'k');title('Envelope');
subplot(2,2,4);plot(t,bp3,'g');title('Biphasic Pulse');
suptitle('Third Filter');
figure
subplot(2,2,1);plot(t,y);title('Original Signal');
subplot(2,2,2);plot(t,f4,'r');title('Bandpass Filter');
subplot(2,2,3);plot(t,e4,'k');title('Envelope');
subplot(2,2,4);plot(t,bp4,'g');title('Biphasic Pulse');
suptitle('Fourth Filter');

%% Problem 4: Wavelet Analysis 

clear;

alpha(1)=(1+sqrt(3))/(4*sqrt(2));
alpha(2)=(3+sqrt(3))/(4*sqrt(2));
alpha(3)=(3-sqrt(3))/(4*sqrt(2));
alpha(4)=(1-sqrt(3))/(4*sqrt(2));
beta(1)=alpha(4);
beta(2)=-alpha(3);
beta(3)=alpha(2);
beta(4)=-alpha(1);

% # of points
N=1024;                     

% input signal 
for n=1:N;m=(n-1)/N;
    g(n)=20*m^2*(1-m)^4*cos(12*pi*m);
end

syms x y
q(x) = 20*y^2*(1-y)^4*cos(12*pi*y);
zt = ztrans(q(x));
%%

sr = 1000; %sample/s
dur = 1; %epoch length
dt = dur/sr; 
f = 10; %frequency in Hz freqz(zt)

t = 0:dt:dur;
inp = cos(f*t*2*pi)+rand(1,sr/dur+1)-0.5; %10Hz cosine w/ gaussian noise, centered at 0 +-.5

figure
plot(t,inp);

%% Problem 5: Wavelets 

clear all 
close all 

t = -50:1:50;
sigma = 10;
psi = (2/(pi^(1/4))*sqrt(3*sigma))*((t.^2)/(sigma^2) - 1).*exp(-(t.^2)/(2*sigma^2));

fm = fft(psi); %numeric solution

sigma = 0.2;
w = -50:1:50;
fa = (-sqrt(8)*(sigma^(5/2))*(pi^(1/4))/sqrt(3))*(w.^2).*exp((-sigma^2)*(w.^2)/2); %analytic solution

figure 
plot(t,psi);xlabel('t');ylabel('\psi (t)');title('Time Domain');
figure
plot(w,fa);xlabel('w');ylabel('\psi (w)');title('Frequency Domain');
figure
freqz(fm);title('Numeric Data');
figure
freqz(fa);title('Analytic Function');

% The data-driven phase plot has 3 phase plateaus, where the phase remains
% relatively constant for a group of frequencies. The analytic function
% does not have this, and instead decreases linearly until reaching 0.5
% normalized frequency. After that the plot ends. The amplitude ratio for
% the data appear to increase in an oddly oscillating manner, while that of
% the analytic solution decreases to an asymptote and is noisy thereafter. 

%% Problem 6: LTI...

clear all
close all

t = 0:0.01:100;
wo = pi/4;
inp = cos(wo*t);
out = (1/4)*((cos(-wo))*exp(1i*wo) + (cos(wo))*exp(1i*-wo) + cos(wo*t)/pi);
figure
hold on
plot(t,inp);
plot(t,out);

wo = 3*pi/4;
inp = cos(wo*t);
out = (1/4)*((cos(-wo))*exp(1i*wo) + (cos(wo))*exp(1i*-wo) + cos(wo*t)/pi);
figure
hold on
plot(t,inp);
plot(t,out);
%% Problem 7: Nonlinearity and Chaos

clear all 
close all

% change in initial condition

a = 4; 
t = 1:100;
x = zeros(1,100);
x(1) = 0.2;
for i = 2:100
    x(i)=a*x(i-1)*(1-(x(i-1)));
end
figure
subplot(3,1,1);plot(t,x);title('IC = 0.255');
x(1) = 0.25;
for i = 2:100
    x(i)=a*x(i-1)*(1-(x(i-1)));
end
subplot(3,1,2);plot(t,x);title('IC = 0.25');ylabel('Output');
x(1)=0.225;
for i = 2:100
    x(i)=a*x(i-1)*(1-(x(i-1)));
end
subplot(3,1,3);plot(t,x);title('IC = 0.225');xlabel('Time');

% small changes yield drastically different dynamics when the system is
% chaotic. 

% change in the a-value

a = 3; 
x = zeros(1,100);
av = zeros(1,10);
m = zeros(1,3);
x(1) = 0.2;
for j = 1:10
    x(1) = rand;
    for i = 2:100
        x(i)=a*x(i-1)*(1-(x(i-1)));
    end
    av(j) = x(length(x));
end
m(1) = std(av);
a = 3.2; 
for j = 1:10
    x(1) = rand;
    for i = 2:100
        x(i)=a*x(i-1)*(1-(x(i-1)));
    end
    av(j) = x(length(x));
end
m(2) = std(av);
a = 4; 
for j = 1:10
    x(1) = rand;
    for i = 2:100
        x(i)=a*x(i-1)*(1-(x(i-1)));
    end
    av(j) = x(length(x));
end
m(3) = std(av);
avals = [3 3.2 4];
figure 
plot(avals,m,'-o');xlabel('a');ylabel('Standard Deviation of Final State');
title('a value change in STDev for varying initial value');

% autocorrelation failure & return plot success

a = 4; 
t = 1:1000;
x = zeros(1,1000); %logistic (non-linear)
xb = zeros(1,1000);
x(1) = 0.4;
for i = 2:1000
    x(i)=a*x(i-1)*(1-(x(i-1)));
    xb(i-1) = x(i);
end
l = zeros(1,1000); %linear 
a = 1;
lb = zeros(1,1000);
l(1) = 0.6;
for i = 2:1000
    l(i) = a*(l(i-1))+(rand*2)-1;
    lb(i-1) = l(i);
end
r = rand(1,1000)*2 - 1; %random process
rb = zeros(1,1000);
for i = 2:1000
    rb(i-1) = r(i);
end

[autox,lagsx] = xcorr(x);
[autol,lagsl] = xcorr(l);
[autor,lagsr] = xcorr(r);

figure 
subplot(2,3,1);plot(lagsl,autol);title('Linear Function Autocorr');xlabel('Lags');ylabel('Correlation Coeff');
subplot(2,3,4);scatter(l,lb,2,'k');title('Linear Function Return Plot');xlabel('l(i)');ylabel('l(i-1)');
subplot(2,3,2);plot(lagsx,autox);title('Non-Linear Function Autocorr');xlabel('Lags');ylabel('Correlation Coeff');
subplot(2,3,5);scatter(x,xb,2,'k');title('Non-Linear Function Return Plot');xlabel('x(i)');ylabel('x(i-1)');
subplot(2,3,3);plot(lagsr,autor);title('Random Function Autocorr');xlabel('Lags');ylabel('Correlation Coeff');
subplot(2,3,6);scatter(r,rb,2,'k');title('Random Function Return Plot');xlabel('r(i)');ylabel('r(i-1)');

% There is obviously a relationship between previous values and the
% following values for the logistic function, but it is missed by the
% autocorrelation since it is non-linear. The linear correlation is found,
% while the random one is not. 
% The return plot, on the other hand, is highly successful in finding a
% relationship. We can see the parabolic curve when plotting one lag away. 
% These examples demonstrate the difficulty of dealing with non-linearity
% and chaos, as small changes can have drastic effects, and relationships
% can be difficult to find with traditional techniques. However, return
% plots appear to be a strong technique to determine relationships. 

%% Problem 8: Modeling Neural Activity 

clear all; close all;
% hodgkin-huxley 

cm = 1; %membrane capacitance (uF/cm^2)
vk = -102; %potassium reversal potential, mv
vna = 25; %sodium reversal potential, mv
vl = -79.4; %leaky reversal potential, mv
gl = 0.3; %leaky capacitance
vrest = 90;
gnamax = 100;
gkmax = 40; 

t = 0;simt = 1000; %ms
dt = 1/25; 
stim = zeros(1,(simt/dt)+1);
sa(1)=20.1;     % Input current (uA/cm^2)
sd(1)=300;      % Duration of the pulse in ms 
sa(2)=0;        % Amplitude (in uA/cm^2) of a second stimulus 
sd(2)=0;        % Duration of the second pulse in ms
dly(1)=50;      % Delay between first stimulus and simulation onset (ms)
dly(2)=0;       % Delay (in ms) between the stimuli #2 and #1 

MV = zeros(1,(simt/dt)+1);
mv = -vrest;

% initialization 

[n,nt]=init_n(mv,vrest);
[m,mt]=init_m(mv,vrest);
[h,ht]=init_h(mv,vrest);
  
c=1;
while t<simt 
    MV(c) = mv;
    
    n=update_n(mv,vrest,dt,n);
    m=update_m(mv,vrest,dt,m);
    h=update_h(mv,vrest,dt,h);
    gk=gkmax*n^4;
	gna=gnamax*m^3*h;
    
    ik=gk*(vk-mv);
	ina=gna*(vna-mv);
	il=gl*(vl-mv);	
	it=ik+ina+il;
    
    if (t <= sd(1)+dly(1))
        if (t >= dly(1))
            it=it+sa(1);
            stim(c)=sa(1);
        end
    end
    if (t <= (dly(1)+dly(2)+sd(2))) 
        if (t >= dly(1)+dly(2))
            it=it+sa(2);
            stim(c)=sa(2); 
        end
    end
    
    dv=it*dt/cm;
    mv=mv+dv;
    
    T(c)=t;t=t+dt;c = c+1;
end

figure;hold;
plot(T,MV);
plot(T,stim,'r')
xlabel('sample# (1/25 ms)');
ylabel(' Membrane Potential (mV)');
ttt=['Single Neuron; Injected Current ' num2str(sa(1)) ' (red)'];
title(ttt);

% integrate-and-fire

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
x=zeros(1,N);                 
x(50:350)=ones(1,301)*Amp;    
%Amp=150; x=randn(1,N)*Amp;      
 
    for n=1:N
        y(n)=(A*y_previous+x(n))/(A+1); % the linear component; see CH 16
        if (y(n) < threshold)          % the nonlinear part
            y_previous=y(n);
            zz(n)=y(n);                 % subthreshold behavior
        else
            y_previous=0;
            zz(n)=AP;                   % action potential
        end
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

% normalization

MV = (MV+abs(min(MV)));
MV = MV/max(MV);
stim = stim/max(stim);
x = x/max(x);
zz = (zz+abs(min(zz)));
zz = zz/max(zz);

figure
subplot(1,2,1);hold on;plot(T,MV);plot(zz,'r');
subplot(1,2,2);hold on;plot(T,stim);plot(x);







