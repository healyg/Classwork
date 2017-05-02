function [info,imax,idt] = PF(n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is a function to use the formula put forward in Obsorne et al 2004
% to determine the information contained in a population of neurons.
%
% INPUT
% 
% n - the # of neurons you would like to check.
%
% OUTPUT
% 
% info - the information contained, found using Shannon's formula
% 
% imax - the maximum information that could be contained in the system
% 
% idt - the uncertainty in the system 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


data = dataload(n);
f = fieldnames(data);

% The data needs to be in 8ms time bins, currently it is in 2
% Since we are looking at it binarily, any point must either be 0 or 1. Any
% value 2 or higher is removed.
% It can only create words for the number of trials in the data set with
% the smallest number of trials (otherwise you lose that letter of the
% word)

s = size(data.first);

setn = zeros(s(1)/4, s(2));

datan = cell(1,n);
sizes = zeros(1,n);

for g = 1:n
    nm = f{g};
    set = data.(nm);
    s = size(set);
    for i = 1:s(1)/4
        for j = 1:s(2)
            for k = 1:s(3)
                setn(i,j,k) = set((4*i),j,k) + set(((4*i)-1),j,k) + set(((4*i)-2),j,k);
                if setn(i,j,k)>=2
                    setn(i,j,k)=1;
                end
            end
        end
    end
    datan{g} = setn;
    sizes(g) = s(3);
end

s = size(datan{1});
m = min(sizes);

words = strings(s(1),s(2),m);
for g = 1:n
    for i = 1:s(1)
        for j = 1:s(2)
            for k = 1:m
                words(i,j,k) = words(i,j,k) + string(datan{g}(i,j,k));
            end
        end
    end
end

% make an array with all possible words, in order to determine the
% probabilities of getting each word

nwords = 2^n;
posswords = strings(nwords,1);
i = 0;
while i < nwords
    i = i+1;
    posswords(i) = '';
    for j = 1:n
        posswords(i) = posswords(i) + string(0+round(rand));
    end
    j = 1;
    while j < nwords
        if j~=i && isequal(posswords(i),posswords(j))
            i = i-1;
            j = nwords;
        end
        j = j+1;
    end
end 
posswords(nwords+1:length(posswords),1) = '';
pwords = zeros(nwords,1);
disp('Words done')
disp(' ')

s = size(words);
for t = 1:nwords
    c = 0;
    for i = 1:s(1)
        for j = 1:s(2)
            for k = 1:s(3)
                if isequal(posswords(t),words(i,j,k))
                    c = c+1;
                end
            end
        end
    end
    pwords(t) = c / (s(1)*s(2)*s(3));
end


%% Compute the probability of finding a certain word in a time and direction

pworddt = zeros(s(1),s(2)); %probability matrix
pwordc = cell(nwords,1);

for t = 1:nwords
    for i = 1:s(1)
        for j = 1:s(2)
            c = 0;
            for k = 1:s(3)
                if isequal(posswords(t),words(i,j,k))
                    c = c+1;
                end
            end
            cur = c/s(3);
            if cur == 0
                cur =1;
            end
            pworddt(i,j) = -cur*log2(cur);
        end
    end
    pwordc{t} = pworddt;
end

tot = zeros(s(1),s(2));
for i = 1:nwords 
    tot = tot + pwordc{i};
end

check = sum(sum(tot));
check = check/(14*s(1));

%%

% Calculate the maximum information 

imax = 0;
for i =1:nwords
    lol = pwords(i);
    if lol == 0
        lol = 1;
    end
    imax = imax + (lol)*log2(lol);
end

imax = -1*imax;
idt = check;
% Calculate the information 
info = imax - idt;


end
