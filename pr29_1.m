% pr29_1.m
% Wim van Drongelen July 2006

%   REQUIRED SUBROUTINES
%    n=init_n(mv,ERest)
%    m=init_m(mv,ERest)
%    h=init_h(mv,ERest)
%    n=update_n(mv,ERest,dt,n);
%    m=update_m(mv,ERest,dt,m);
%    h=update_h(mv,ERest,dt,h);

clear;


% parameters
simt=501;       % Simulation time in ms 
sa(1)=input(' Amplitude (in uA/cm^2) of the first stimulus (e.g. 15)  : ');
sd(1)=100;      % Duration of the pulse in ms 
sa(2)=0;        % Amplitude (in uA/cm^2) of a second stimulus 
sd(2)=0;        % Duration of the second pulse in ms
dly(1)=50;      % Delay between first stimulus and simulation onset (ms)
dly(2)=0;       % Delay (in ms) between the stimuli #2 and #1 

% Initialization
% --------------
t=0;                        % Start Time
dt=1/25;                    % Simulation Interval
sz=simt/dt;                 % Size of the Arrays
MV=zeros(1,sz+1);           % Array for the Membrane Potential
stim=zeros(1,sz+1);         % Array for the Membrane Current

ERest=90;                   % Abs value of Resting Potential relative to the outside in mV
EK=-102;                    % Equilibrium Potentials in mV re to the Outside
ENa=25;
EL=-79.4;
mv=-ERest;                  % Current Potential relative to outside in mV (initial value - ERest)

Cm=1;                       % Membrane Capacitance (uF/cm^2)
gnaBAR=100;                 % Max Conductivities in mmho/cm^2 (mmho=mS)
gkBAR=40;
gL=0.3;

% Steady State variables (Using Subroutines [x,xt]=init_x)
% --------------------------------------------------------
[n,nt]=init_n(mv,ERest);
[m,mt]=init_m(mv,ERest);
[h,ht]=init_h(mv,ERest);
  
% ------------------------------
%  SIMULATION LOOP Current Clamp
% ------------------------------
k=1;
while (t<simt);

 	MV(k)=mv;   % Store last Membrane Potential in the Array
     
	% Update the conductances (Using Subroutines w. Euler's Method)
	% -----------------------
    n=update_n(mv,ERest,dt,n);
    m=update_m(mv,ERest,dt,m);
    h=update_h(mv,ERest,dt,h);
    gk=gkBAR*n^4;
	gna=gnaBAR*m^3*h;

	% Update the currents
	% -------------------
	ik=gk*(EK-mv);
	ina=gna*(ENa-mv);
	il=gL*(EL-mv);	
	it=ik+ina+il;

    % Is there a stimulus to add to the currents?
	% -------------------------------------------
		if (t <= sd(1)+dly(1)); 
            if (t >= dly(1));
                it=it+sa(1);
                stim(k)=sa(1);
            end;
        end;
		if (t <= (dly(1)+dly(2)+sd(2))); 
            if (t >= dly(1)+dly(2)); 
                it=it+sa(2);
                stim(k)=sa(2); 
            end; 
        end;

	% The new membrane potential
	% --------------------------
  	dv=it*dt/Cm;
    mv=mv+dv;
    
	% Update t,k
	% --------
	T(k)=t;t=t+dt;k=k+1;
end;

% Plot Result
% -----------
figure;hold;
plot(T,MV);
plot(T,stim,'r')
xlabel('sample# (1/25 ms)');
ylabel(' Membrane Potential (mV)');
ttt=['Single Neurone; Injected Current ' num2str(sa(1)) ' (red)'];
title(ttt);
