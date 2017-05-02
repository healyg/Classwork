clear coords 

coords = cell(1,length(tseq));

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
while c<=length(tseq)
    last = coords{c-1};
    coords{c} = coords{c-2};
    
    %makes sure the snake doesn't double back
    while isequal(coords{c},coords{c-2})
        r = rand;
        if r>0.25
            coords{c} = [last(1) last(2)+1]; %move right
        elseif r>=0.25 && r<0.5
            coords{c} = [last(1)+1 last(2)]; %move up
        elseif r>=0.5 && r<0.75
            coords{c} = [last(1) last(2)-1]; %move left
        else 
            coords{c} = [last(1)-1 last(2)]; %move down
        end
    end
    
    %makes sure the snake is in an unoccupied spot
    d = 1;
    while d<c
        if isequal(coords{d},coords{c})
            r = rand;
            if r>0.25
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
    
    cur = coords{c};
    t=0;
    
    %makes sure the snake hasn't hit a dead end (if c is not the last
    %number)
    if c == length(tseq)
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
    while t==3
        q = q-1;
        p = 0;
        for i = 1:q
            if isequal(coords{i},[last(1) last(2)+1]) || isequal(coords{i},[last(1) last(2)-1]) ...
                    || isequal(coords{i},[last(1)+1 last(2)]) || isequal(coords{i},[last(1)-1 last(2)])
                p = p+1;
            end
        end
        if p <2
            t=0;
            q=c;
        else 
            q = q-1;
        end
    end
    c = q+1;
end

%%
%checking the max 

mx1 = 0;
mx2 = 0;
mn1 = 0;
mn2 = 0;
for i = 2:length(coords)
    cur = coords{i};
    if cur(1)>mx1 
        mx1 = cur(1);
    end 
    if cur(2)>mx2
        mx2 = cur(2);
    end
    if cur(1)<mn1
        mn1 = cur(1);
    end 
    if cur(2)<mn2
        mn2 = cur(2);
    end
end

borders = [mx1 mx2; mn1 mn2];
m = borders(1,1)- borders(2,1);
n = borders(1,2) - borders(2,2);



