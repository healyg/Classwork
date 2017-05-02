function [newcoords] = reconstruct(coords,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function reconstructs the lattice (coords) from point n
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    error('You need a previous coordinate cell vector as well as a starting' ...
        + ' point for the reconstruction.')
elseif n>length(coords)
    error('n is larger than the length of your coordinate cell')
end 

c = n;

if c ==1 || c==2
    coords{1} = [0 0];
    last = coords{1};
    r = rand;
    if r<0.25
        coords{2} = [last(1) last(2)+1]; %move right
    elseif r>=0.25 && r<0.5
        coords{2} = [last(1)+1 last(2)]; %move up
    elseif r>=0.5 && r<0.75
        coords{2} = [last(1) last(2)-1]; %move left
    else 
        coords{2} = [last(1)-1 last(2)]; %move down
    end
    c = 3;
end

while c<=length(coords)
    last = coords{c-1};
    coords{c} = coords{c-2};
    
    %makes sure the snake doesn't double back
    while isequal(coords{c},coords{c-2})
        r = rand;
        if r<0.25
            coords{c} = [last(1) last(2)+1]; %move right
        elseif r>=0.25 && r<0.5
            coords{c} = [last(1)+1 last(2)]; %move up
        elseif r>=0.5 && r<0.75
            coords{c} = [last(1) last(2)-1]; %move left
        else 
            coords{c} = [last(1)-1 last(2)]; %move down
        end
    end
    cur = coords{c};
    t=0;
    
    %makes sure the snake hasn't hit a dead end (if c is not the last
    %number)
    if c == length(coords)
    else
        for i =1:c
            if isequal(coords{i},[cur(1) cur(2)+1]) || isequal(coords{i},[cur(1) cur(2)-1]) ...
                    || isequal(coords{i},[cur(1)+1 cur(2)]) || isequal(coords{i},[cur(1)-1 cur(2)])
                t = t+1;
            end
        end
    end
    
    %corrects the snake if it has hit a dead end
    q = c;
    while t==4
        p = 0;
        last = coords{q-1};
        for i = 1:q-1
            if isequal(coords{i},[last(1) last(2)+1]) || isequal(coords{i},[last(1) last(2)-1]) ...
                    || isequal(coords{i},[last(1)+1 last(2)]) || isequal(coords{i},[last(1)-1 last(2)])
                p = p+1;
            end
        end
        if p <2
            t=0;
        else 
            q = q-1;
        end
    end
    %makes sure the snake is in an unoccupied spot
    if q == c
        d = 1;
        while d<c
            if isequal(coords{d},coords{c})
                r = rand;
                if r<0.25
                    coords{c} = [last(1) last(2)+1]; %move right
                elseif r>=0.25 && r<0.5
                    coords{c} = [last(1)+1 last(2)]; %move up
                elseif r>=0.5 && r<0.75
                    coords{c} = [last(1) last(2)-1]; %move left
                else 
                    coords{c} = [last(1)-1 last(2)]; %move down
                end
                d=0;
            end
            d = d+1;
        end
        c = q+1;
    else 
        c=q;
    end
end

newcoords = coords;

end