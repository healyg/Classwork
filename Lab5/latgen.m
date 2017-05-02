function [coords] = latgen(seq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Generates a random lattice structure for the protein with boundaries of m
% and n. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coords = cell(1,length(seq));

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
while c<=length(seq)
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
    if c == length(seq)
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
        last = coords{c-1};
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

x  = zeros(1,length(coords));
y = zeros(1,length(coords));
for i = 1:length(coords)
    cur = coords{i};
    x(i) =cur(1);
    y(i) = cur(2);
end

figure 
plot(x,y);


end