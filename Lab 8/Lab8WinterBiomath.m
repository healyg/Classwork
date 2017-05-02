%% Lab 8: Gillespie simulations of continuous-time Markov processes
%
% Name: Garrett Healy
%
% This algorithm is used to simulate the behavior of a CTMC where there are
% different species of individuals and random transitions can be described
% by *reactions* that modify the numbers of individuals in at least one
% species. Before implementing a Gillespie simulation, write down all
% possible reactions and how they affect each species, and assign each one
% an index.
%
% Gillespie algorithm
%
%% Part 1: Birth and death processes
% 1.1 *Pure birth process*
%
% A pure birth process is a population where each individual produces one
% offspring with a given rate, independent of the others. In the language
% of the master equation, there is only one reaction possible: increasing
% the size of the population by 1, with effective rate (propensity) N.
% Implement the Gillespie algorithm for a pure birth process, initially
% starting with population of 1, with birth rate of 0.1 per day, and run it
% for 100 days. (Hint: since there is only one reaction, you don't need to
% decide which reaction takes place, only the time of the next birth.)
% The mean of this process is governed by the ODE:
%    dN/dt = 0.1N, with N(0)=1
% Plot the solution of this ODE on the same plot as the 10 stochastic
% realizations of the pure birth process.

clear all
close all

figure 
hold on
for i = 1:10
    M = 1;
    N = 1; 

    c = 0.1; 
    X = 1; 

    t = 0; %sets initial time 
    tplot = zeros(1); 
    tmax = 100;
    k = 2; 
    xplot = ones(1); 
    tau=1; 

    while (t<tmax) && tau~=0
        r1 = rand;
        r2 = rand; 
        a = c*X;  
        A = sum(a); 
        tau = -log(r1)/A;
        t = t+tau; 
        X = X+1; 
        tplot(k) = t; 
        xplot(k) = X; 
        k = k+1; 
    end
    plot(tplot,xplot)
end

[t2,y] = ode45(@(t,y) 0.1*y, [0 100], 1);
plot(t2,y,'-o')
hold off
legend('1','2','3','4','5','6','7','8','9','10','ODE');

% 1.2 *Pure death process*
%
% Modify your code to model a pure death process with a stochastic death
% rate for each individual, starting with a population of 1000 individuals,
% with a death rate of 0.05 per day, and run 10 simulations for 100 days.
% Compute the solution for the mean number of individuals from the ODE:
%
%   dN/dt = -0.05N, with N(0)=1000
%
% Plot the solution of this ODE on the same plot as the 10 stochastic
% realizations of the pure birth process.

figure 
hold on
for i = 1:10 
    M = 1;
    N = 1; 

    c = 0.05; 
    X = 1000; 

    t = 0; %sets initial time 
    tplot = zeros(1); 
    tmax = 100;
    k = 2; 
    xplot = 1000*ones(1); 
    tau=1; 

    while (t<tmax) && tau~=0
        r1 = rand;
        r2 = rand; 
        a = c*X;  
        A = sum(a); 
        tau = -log(r1)/A;
        t = t+tau; 
        X = X-1; 
        tplot(k) = t; 
        xplot(k) = X; 
        k = k+1; 
    end
    plot(tplot,xplot)
end

[t2,y] = ode45(@(t,y) -0.05*y, [0 100], 1000);
plot(t2,y,'-o')
hold off
legend('1','2','3','4','5','6','7','8','9','10','ODE');


% 1.3 *Extinction time statistics*
% Run 100 simulations until all the individuals die and compute the
% mean time to extinction and the standard deviation. How well is it
% predicted by the deterministic mean solution (exponential decay) dipping
% below 1 individual?  

extt = zeros(1,100);

for i = 1:100 
    M = 1;
    N = 1; 

    c = 0.05; 
    X = 1000; 

    t = 0; %sets initial time 
    tplot = zeros(1); 
    tmax = 200;
    k = 2; 
    xplot = 1000*ones(1); 
    tau=1; 

    while (t<tmax) && tau~=0 && X>1
        r1 = rand;
        r2 = rand; 
        a = c*X;  
        A = sum(a); 
        tau = -log(r1)/A;
        t = t+tau; 
        X = X-1; 
        tplot(k) = t; 
        xplot(k) = X; 
        k = k+1; 
    end
    extt(i) = t; 
end

meanext = mean(extt);
stdext = std(extt);

[t2,y] = ode45(@(t,y) -0.05*y, [0 200], 1000);
for i = 1:length(y) 
    if y(i) > 1 
        extt2 = t2(i);
    end
end

% The extinction time is fairly well modeled, but after running it a few
% times it appear to predict a slightly quicker extinction that we would
% expect simply by looking at the ODE. 


%% Part 2. Predator-prey stochastic model
% Now let us model a system with more than one reaction, namely a
% predator-prey system. Let us call the number of prey X and the number of
% predators Y . The three different processes are: 
%
% # Birth of one prey with rate c1 = 20 and propensity a1 = c1*X
% # Death of one predator with rate c2 = 10 and propensity a2 = c2*Y
% # Predation (death of one prey and birth of one predator) with rate c3
% =0.01 and propensity  a3 = c3*X*Y 
%
%
% 2.1 *ODE model*
%
% Write down the two-variable differential equation (Lotka-Volterra) based
% on the model above. Find the nullclines and the fixed points and predict
% the qualitative behavior of the system. This is the model for the mean
% number of predators and prey.

clear all;clc; 
close all

c1 = 20;
c2 = 10; 
c3 = 0.01;

options = odeset('Refine',10,'MaxStep',0.1);

tlimit=1000;
tspan=[0 tlimit]; 
time=zeros(tlimit,1);

figure; 
hold all;
xmin= 999; 
xmax= 1001;
ymin= 1999;
ymax= 2001;
dx =0.05; 
dy =0.05; 
[X,Y]=meshgrid(xmin:dx:xmax, ymin:dy:ymax); 

dy = Y*(-c2 + c3*X);
dx = X*(c1 - c3*Y);

%   dS/dt = S*(a-S-b*R)
%   dR/dt = R*(c-R-d*S)

quiver(X,Y,dx,dy);
xlim([xmin xmax]); 
ylim([ymin ymax]);

% Fixed point, calculated by hand, is (1000,2000). The nullcline is a
% circle, therefore the populations will alternatively increase and
% decrease in a stable fashion (at least in theory).

%% 2.2 *Stochastic model*
%
% Implement the Gillespie algorithm for the three-reaction predator-prey
% model. Using the parameters above and initial conditions X0 = 500, Y0 =
% 500, propagate the system for a large number of steps (e.g. 100000)
% and plot the trajectory in phase space (plot(X,Y)). Describe whether it
% agrees with the qualitative behavior of the continuous ODE model, and how
% the stochasticity manifests itself.


M = 3;
N = 2; 

c = [20 10 0.01]; 
X = [500 500]; 

t = 0; %sets initial time 
tplot = zeros(1); 
tmax = 100000;
k = 2; 
xplot = 500*ones(1); 
yplot = 500*ones(1);
tau=1; 

while (t<tmax) && tau~=0
    r1 = rand;
    r2 = rand; 
    a = [c(1)*X(1); c(2)*X(2);c(3)*X(1)*X(2)]; 
    A = sum(a); 
    tau = -log(r1)/A;
    t = t+tau; 
    check = r2*A;
    if check <= a(1)
        X(1) = X(1) + 1;
    elseif a(1) < check && check <= a(1) + a(2)
        X(2) = X(2)-1;
    elseif a(1) + a(2) < check && check <= A
        X(1) = X(1) - 1; 
        X(2) = X(2) + 1;
    end 
    tplot(k) = t; 
    xplot(k) = X(1); 
    yplot(k) = X(2);
    k = k+1; 
end
plot(xplot,yplot)

% The plot shows significant similarities to the proposed nullcline/ phase
% plane dynamics. The overall shape is the same, and the fixed point can be
% seen in the center of the system. The stochasticity manifests itself in
% different speeds in unraveling. 

% 2.3 *Extinction of species*
%
% Change the initial populations to be 10 each (of predatory and prey) so
% that you see extinction of one of the species in some simulations. Report
% the fraction of 200 simulations which go extinct in 1000 steps. Change
% the initial populations to 20 each, run 200 simulations, and again report
% the fraction that goes extinct in 1000 steps.

ct = 0;

for i = 1:200
    M = 3;
    N = 2; 

    c = [20 10 0.01]; 
    X = [10 10]; 

    t = 0; %sets initial time 
    tplot = zeros(1); 
    tmax = 1000;
    k = 2; 
    xplot = 10*ones(1); 
    yplot = 10*ones(1);
    tau=1; 

    while (t<tmax) && tau~=0 && k < tmax
        r1 = rand;
        r2 = rand; 
        a = [c(1)*X(1); c(2)*X(2);c(3)*X(1)*X(2)]; 
        A = sum(a); 
        tau = -log(r1)/A;
        t = t+tau; 
        check = r2*A;
        if check <= a(1)
            X(1) = X(1) + 1;
        elseif a(1) < check && check <= a(1) + a(2)
            X(2) = X(2)-1;
        elseif a(1) + a(2) < check && check <= A
            X(1) = X(1) - 1; 
            X(2) = X(2) + 1;
        end 
        tplot(k) = t; 
        xplot(k) = X(1); 
        yplot(k) = X(2);
        k = k+1; 
    end
    if X(1) == 0 || X(2) == 0
        ct=ct+1;
    end
end 

frac = ct/200;

ct2 = 0;

for i = 1:200
    M = 3;
    N = 2; 

    c = [20 10 0.01]; 
    X = [20 20]; 

    t = 0; %sets initial time 
    tplot = zeros(1); 
    tmax = 1000;
    k = 2; 
    xplot = 20*ones(1); 
    yplot = 20*ones(1);
    tau=1; 

    while (t<tmax) && tau~=0 && k < tmax
        r1 = rand;
        r2 = rand; 
        a = [c(1)*X(1); c(2)*X(2);c(3)*X(1)*X(2)]; 
        A = sum(a); 
        tau = -log(r1)/A;
        t = t+tau; 
        check = r2*A;
        if check <= a(1)
            X(1) = X(1) + 1;
        elseif a(1) < check && check <= a(1) + a(2)
            X(2) = X(2)-1;
        elseif a(1) + a(2) < check && check <= A
            X(1) = X(1) - 1; 
            X(2) = X(2) + 1;
        end 
        tplot(k) = t; 
        xplot(k) = X(1); 
        yplot(k) = X(2);
        k = k+1; 
    end
    if X(1) == 0 || X(2) == 0
        ct2=ct2+1;
    end
end 

frac2 = ct2/200;

% The fraction is generally significantly lower when initial population is
% increased. 

% 2.4 *Fluctuations around the equilibrium*
%
% Set your initial populations to the non-zero fixed point that you found in 2.1
% and re-run the simulation for 10000 steps. Because this is a stochastic
% process, we do not expect that the populations will remain at the fixed
% points as the ODE predicts. Using the results of your simulation, calculate the
% variance in the Euclidean distance from the fixed point across 10
% simulations as a function of time. Plot variance vs. t. What is the
% relationship? For a pure diffusion process, the variance will increase
% linearly as a function of time. Is that the case here? 

dist = zeros(10,1);
for i = 1:10
    M = 3;
    N = 2; 

    c = [20 10 0.01]; 
    X = [1000 2000]; 

    t = 0; %sets initial time 
    tplot = zeros(1); 
    tmax = 10000;
    k = 2; 
    xplot = 500*ones(1); 
    yplot = 500*ones(1);
    tau=1; 

    while (t<tmax) && tau~=0 && k<tmax
        r1 = rand;
        r2 = rand; 
        a = [c(1)*X(1); c(2)*X(2);c(3)*X(1)*X(2)]; 
        A = sum(a); 
        tau = -log(r1)/A;
        t = t+tau; 
        check = r2*A;
        if check <= a(1)
            X(1) = X(1) + 1;
        elseif a(1) < check && check <= a(1) + a(2)
            X(2) = X(2)-1;
        elseif a(1) + a(2) < check && check <= A
            X(1) = X(1) - 1; 
            X(2) = X(2) + 1;
        end 
        tplot(k) = t; 
        xplot(k) = X(1); 
        yplot(k) = X(2);
        dist(i,k) = sqrt((X(1) - 1000)^2 + (X(2) - 2000)^2);
        k = k+1; 
    end
end
plot(tplot,dist)

% No, that is generally not the case with the functions I see. They
% occasionally look close to linear, given that they are stochastic and
% will never be perfectly linear, but most of the time they are far from
% it, oscillating up and down often. 

% 2.5 *Comparison of amplitude and frequencies of oscillations*
%
% This exercise will require that you run the simulation twice: Once with
% predator and prey initial populations equal to 500, and once with initial
% populations set at the coexistence fixed point, and run the simulation
% for 100000 steps. For each simulation, plot the predator population and
% the prey population as a function of time on the same plot. Using your
% plots, compare the amplitudes and frequencies of population oscillation
% when the initial populations are equal to 500 and when the initial
% populations are equal to the fixed points. How are they similar or
% different? Are amplitude and frequency independent of one another? (this
% is a hallmark of linear oscillations)

close all

c = [20 10 0.01]; 
X = [500 500]; 

t = 0; %sets initial time 
tplot = zeros(1); 
tmax = 100000;
k = 2; 
xplot = 500*ones(1); 
yplot = 500*ones(1);
tau=1; 

while (t<tmax) && tau~=0 && k<tmax
    r1 = rand;
    r2 = rand; 
    a = [c(1)*X(1); c(2)*X(2);c(3)*X(1)*X(2)]; 
    A = sum(a); 
    tau = -log(r1)/A;
    t = t+tau; 
    check = r2*A;
    if check <= a(1)
        X(1) = X(1) + 1;
    elseif a(1) < check && check <= a(1) + a(2)
        X(2) = X(2)-1;
    elseif a(1) + a(2) < check && check <= A
        X(1) = X(1) - 1; 
        X(2) = X(2) + 1;
    end 
    tplot(k) = t; 
    xplot(k) = X(1); 
    yplot(k) = X(2);
    k = k+1; 
end
figure
hold on
plot(tplot,xplot)
plot(tplot,yplot)
legend('X','Y')
hold off

c = [20 10 0.01]; 
X = [1000 2000]; 

t = 0; %sets initial time 
tplot = zeros(1); 
tmax = 100000;
k = 2; 
xplot = 1000*ones(1); 
yplot = 2000*ones(1);
tau=1; 

while (t<tmax) && tau~=0 && k<tmax
    r1 = rand;
    r2 = rand; 
    a = [c(1)*X(1); c(2)*X(2);c(3)*X(1)*X(2)]; 
    A = sum(a); 
    tau = -log(r1)/A;
    t = t+tau; 
    check = r2*A;
    if check <= a(1)
        X(1) = X(1) + 1;
    elseif a(1) < check && check <= a(1) + a(2)
        X(2) = X(2)-1;
    elseif a(1) + a(2) < check && check <= A
        X(1) = X(1) - 1; 
        X(2) = X(2) + 1;
    end 
    tplot(k) = t; 
    xplot(k) = X(1); 
    yplot(k) = X(2);
    k = k+1; 
end
figure 
hold on
plot(tplot,xplot)
plot(tplot,yplot)
legend('X','Y')
hold off

% In both graphs the prey values never reached as high as the predator
% values, but seemed to oscillate as the same frequency in the 500/500
% state. For the steady-state system, the overall frequency appeared to be
% the same, but it also had many more errant frequencies that are difficult
% to track. It would be difficult to tell which is which however, so I
% assume it also has a shared frequency. And while the 500/500 state had
% distinctly larger amplitudes for predators, that was difficult to observe
% in the steady-state: the values for predator population were always
% larger, but their amplitudes appeared similar. Amplitude and frequency do
% generally appear to be independent of one another. 










