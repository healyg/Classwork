%% Lab 5: Simulated annealing Monte Carlo for optimization
%% Part 1: MC for TSP
%
% * Input: a list of lattitudes (x-coordinates), a list of longitudes (y-coordinates), intial path
%
% * Do: find the shortest possible path through all N vertices
%
% * Output: ordered list of the N indices of the verticies that produces the shortest
% path.
%
%
% 1.1: Write a function that takes in the list of coordinates and the
% list of the indices and outputs the total length of the path. This will
% be the objective function. Use the following list of 5 "city" coordinates to
% test your objective function for different order of verticies:

clear all; 
close all; 

xcoord=[0 2 -1 2 -2];
ycoord=[0 1 1 3 1];

load('US_cities.mat');

tlength = objfunk(xcoord,ycoord,[2 3 4 5 1]);

%%
% 1.2: Write a Monte Carlo simulation that uses the objective function
% from task 1.1, and uses the following proposal distribution: take the
% current ordered list of vertices and swap two randomly selected ones
% (this can be done in one line). The proposal candidate will be accepted
% or rejected according to the Metropolis criterion with temperater
% parameter T. Use the test coordinates provided above, and then
% run your simulation on a larger data set that you load from the file
% US_cities.mat (which contains coordinates for 22 US cities.) Run the
% simulation for different values of T and oberve the effect on convergence. 

index = 1:22;
l = length(index);

tlength = objfunk(lat,long,index);
startlength = tlength;
c = 0;
T = 1;
while c<10000
    swap1 = round(1 + (l-1)*rand);
    swap2 = round(1 + (l-1)*rand);
    ind1 = index(swap1);
    ind2 = index(swap2);
    index(swap1) = ind2;
    index(swap2) = ind1;
    test = objfunk(lat,long,index);
    if test<tlength
    else
        p = exp((tlength - test)/T);
        maxp = exp(tlength/T);
        p = p/maxp;
        r = rand;
        if r<p
        else 
            index(swap1) = ind1;
            index(swap2) = ind2;
        end
    end
    tlength = objfunk(lat,long,index);
    c=c+1;
end

endlength = objfunk(lat,long,index);
diff = startlength - endlength;

newlocs = cell(char(0),length(index),1);
for i = 1:length(index)
    ind = index(i);
    newlocs(i) = names(ind);
end

disp(endlength);
disp(newlocs);

%%
% 1.3 Implement a simulated annealing cooling schedule to change
% temperature T from an initial temperature (T_init) down to final
% (T_final) by decreasing it every interation by a constant step (dT).
% Change those three parameters and see which results in improved
% convergence, first for the 5-city test data set, and then for the data
% set of 22 US cities.

index = 1:22;
l = length(index);

tlength = objfunk(lat,long,index);
startlength = tlength;
c = 1;
T_int = 1;
T_final = 0.001;
dT = 0.00001;
while c< (T_int-T_final)/dT
    T= T_int - dT*c;
    swap1 = round(1 + (l-1)*rand);
    swap2 = round(1 + (l-1)*rand);
    ind1 = index(swap1);
    ind2 = index(swap2);
    index(swap1) = ind2;
    index(swap2) = ind1;
    test = objfunk(lat,long,index);
    if test<tlength
    else
        p = exp((tlength - test)/T);
        maxp = exp(tlength/T);
        p = p/maxp;
        r = rand;
        if r<p
        else 
            index(swap1) = ind1;
            index(swap2) = ind2;
        end
    end
    tlength = objfunk(lat,long,index);
    c=c+1;
end

diff = startlength - tlength;

newlocs = cell(char(0),length(index),1);
for i = 1:length(index)
    ind = index(i);
    newlocs(i) = names(ind);
end

disp(newlocs);


%% Part 2: Monte Carlo simulation of protein folding on a lattice
% In this part you will use the HP lattice model to simulate the process of
% protein folding. For this model, assume the following:
%
% 1. The "protein" is an unbranched chain of amino acid residues, each of
% which can reside on a single point on a square lattice in the plane,
% which can be indexed by integer coordinates.
%
% 2. The chain can be arranged in any conformation on the lattice, as long as:
% a) consecutive amino acids reside on adjacent lattice points, and 
% b) the chain does not run into itself (each lattice point can only be
% occupied by one residue).
%
% 3. there are two kinds of residues: H (hydrophobic) and P (polar). The
% energy of a protein conformation is given by the negative sum of all
% adjacent hydrophobic pairs (which means those off by 1 either in
% horizontal or vertical coordinate, diagonal neighbors don't count), all
% other pairs have energy of zero.
%
% Your task it to write a script that uses Monte Carlo simulated annealing to
% find the lowest energy conformation for a given "protein sequence."
%
% * Inputs: Sequence of H and P of length N, objective function for measuring
% total energy of a conformation, and the parameters
% determining the cooling schedule for simulated annealing and
% the maximum number of iterations.
%
% * Do: find the arrangement of the N verticies that meets the conditions
% specified above are met, such that the objective function is minimized.
%
% * Return: The coordinates which minimize the objective function, and the
% value of the objective function.
%
% 2.1 Write a function to calculate the energy of a given conformation
% (negative of the total number of HH pairs on adjacent lattice points).
% Take care not to double count! The inputs should be the sequence and the
% set of coordinates for each "residue". Test your code on the sequence
% 'HHPPHH' by inputting several different conformations (sets of coordinates).

clear all;clc;
close all; 

load('test_seqs.mat')

tseq = 'HHPPHH';
index = [3 4;1 2;5 6];
[a,b]=lattice(tseq,index);


%index = [1:8;9:16;17:24;25:32;33:40;41:48];

% 2.2. Write a function to randomly generate an initial conformation, that is an
% array of length N of two-dimensional coordinates of the "residues". Here
% is a possible outline: start residue 1 at the origin (coordinate (0,0))
% and randomly choose one of the four lattice neighbors for residue 2, and
% then for each subsequent residue randomly chooses only the unoccupied
% neighbors to place the next residue. Note that this can lead into a dead
% end, in which none of the four neighbors are available, in which case
% your code should go back one step and make a different placement choice,
% and if there is only one, take one more step back, until there is a
% different placement that can be chosen. Hint: use the function randi(N)
% to generate random integers from 1 to N, inclusive.

[coords] = latgen(seq1);

%%
% 2.3 Write a function to randomly generate a candidate conformation, given a
% current conformation. This can be combined with the first function to
% generate an initial conformation: for example, you can send the current
% conformation to this function, then choose at random at which residue to
% start rebuilding (anywhere from 1 - total rebuild, to N - no rebuilding)
% and then proceed with the same algorithm you used to generate the random
% initial conformation.

[coords] = latgen(seq1);
[newcoords] = reconstruct(coords,40);

% 2.4. Write a script with Monte Carlo code that calls the above functions
% to generate a new candidate conformation, then evaluates its energy
% compared to the current conformation, and decides to accept or reject it,
% according to the Metropolis criterion with a set temperature parameter T.
% This code should have a criterion for when to stop, e.g. after Nmax total
% iterations, or after Nstall iterations of not finding any lower energy
% conformation. Keep track of the lowest energy conformation, and after the
% simulation is finished, plot the initial and best conformations, and
% write out the number of iterations to finding the best conformation, plus
% the final and initial energies. Test your code on the sequence
% 'HHPPHH'.  First, predict its optimal folded conformation by inspection
% (identify the correct folded configuration, although it has some
% redundancy). Once you have tested the code and are convinced your code
% is bug-free, run it 10 times with different T parameters and report a)
% the lowest-energy conformation and its energy; b) the number of Monte Carlo
% steps it takes to find it. Comment on which T value appears optimal.


seq = 'HHPPHH';
[coords] = latgen(seq);
[placement,energy] = lattice(seq,coords);
startcoords = coords;
finalcoords = coords;
startplacement = placement;
finalplacement = placement;
c = 0;
T = 1;
while c<10000
    r = randi([1,length(seq)]);
    [coords] = reconstruct(coords,r);
    [newplacement, newenergy] = lattice(seq,coords);
    if newenergy<energy
        energy = newenergy;
        finalcoords = coords;
        finalplacement = newplacement;
    else
        p = exp((energy - newenergy)/T);
        maxp = exp(energy/T);
        p = p/maxp;
        r = rand;
        if r<p
        else 
            energy = newenergy;
            finalcoords = coords;
            finalplacement = newplacement;
        end
    end
    c=c+1;
end

x  = zeros(1,length(finalcoords));
y = zeros(1,length(finalcoords));
for i = 1:length(finalcoords)
    cur = finalcoords{i};
    x(i) =cur(1);
    y(i) = cur(2);
end

figure 
plot(x,y);




%% 2.5 Please find three sequences (which were used to test optimization
% schemes in the 1990s) in the file test_seqs.mat, saved as character
% strings of Hs and Ps in variables seq1, seq2, and seq3. For each
% sequence, run batches of 10 different MC runs with the optimal value of T
% you found above. Report the lowest energy conformations that you found
% (show plots of the folded configurations), their energies and the
% statistics for the number of steps until reaching them. How often do
% these optimization runs fail to reach the minimum energy conformation?
% (make sure you set the max iterations large enough to explore the state
% space). 

seq = seq1;
[coords] = latgen(seq);
[placement,energy1] = lattice(seq,coords);
startcoords1 = coords;
finalcoords1 = coords;
startplacement1 = placement;
finalplacement1 = placement;
c = 0;
T = 1;
while c<10000 
    r = randi([1,length(seq)]);
    [coords] = reconstruct(coords,r);
    [newplacement, newenergy] = lattice(seq,coords);
    if newenergy<energy1
        energy1 = newenergy;
        finalcoords1 = coords;
        finalplacement1 = newplacement;
    else
        p = exp((energy1 - newenergy)/T);
        maxp = exp(energy1/T);
        p = p/maxp;
        r = rand;
        if r<p
        else 
            energy1 = newenergy;
            finalcoords1 = coords;
            finalplacement1 = newplacement;
        end
    end
    c=c+1
end

x  = zeros(1,length(finalcoords1));
y = zeros(1,length(finalcoords1));
for i = 1:length(finalcoords1)
    cur = finalcoords1{i};
    x(i) =cur(1);
    y(i) = cur(2);
end

figure 
plot(x,y);

seq = seq2;
[coords] = latgen(seq);
[placement,energy2] = lattice(seq,coords);
startcoords2 = coords;
finalcoords2 = coords;
startplacement2 = placement;
finalplacement2 = placement;
c = 0;
T = 1;
while c<10000 && t < 10
    r = randi([1,length(seq)]);
    [coords] = reconstruct(coords,r);
    [newplacement, newenergy] = lattice(seq,coords);
    if newenergy<energy
        energy2 = newenergy;
        finalcoords2 = coords;
        finalplacement2 = newplacement;
    else
        p = exp((energy2 - newenergy)/T);
        maxp = exp(energy2/T);
        p = p/maxp;
        r = rand;
        if r<p
        else 
            energy2 = newenergy;
            finalcoords2 = coords;
            finalplacement2 = newplacement;
        end
    end
    c=c+1
end

x  = zeros(1,length(finalcoords2));
y = zeros(1,length(finalcoords2));
for i = 1:length(finalcoords2)
    cur = finalcoords2{i};
    x(i) =cur(1);
    y(i) = cur(2);
end

figure 
plot(x,y);

seq = seq3;
[coords] = latgen(seq);
[placement,energy] = lattice(seq,coords);
startcoords = coords;
finalcoords = coords;
startplacement = placement;
finalplacement = placement;
c = 0;
T = 1;
while c<10000 && t < 10
    r = randi([1,length(seq)]);
    [coords] = reconstruct(coords,r);
    [newplacement, newenergy] = lattice(seq,coords);
    if newenergy<energy
        energy = newenergy;
        finalcoords = coords;
        finalplacement = newplacement;
    else
        p = exp((energy - newenergy)/T);
        maxp = exp(energy/T);
        p = p/maxp;
        r = rand;
        if r<p
        else 
            energy = newenergy;
            finalcoords = coords;
            finalplacement = newplacement;
        end
    end
    c=c+1;
end

x  = zeros(1,length(finalcoords));
y = zeros(1,length(finalcoords));
for i = 1:length(finalcoords)
    cur = finalcoords{i};
    x(i) =cur(1);
    y(i) = cur(2);
end

figure 
plot(x,y);
