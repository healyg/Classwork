function [lat,objvalue] = lattice(seq,index)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Solve the lattice structure problem associated with Lab 5
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2
    error('Two inputs are needed, see function code for more info')
end

if iscell(index)
    mx1 = 0;
    mx2 = 0;
    mn1 = 0;
    mn2 = 0;
    for i = 2:length(index)
        cur = index{i};
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
    m = (borders(1,1)- borders(2,1))+1;
    n = (borders(1,2) - borders(2,2))+1;
    temp = zeros(m,n);
    for i= 1:length(index)
        new = index{i};
        new(1) = new(1) + abs(mn1)+1;
        new(2) = new(2) + abs(mn2)+1;
        temp(new(1),new(2)) = i;
    end
end

index = temp;

s = size(index);

m = s(1); 
n = s(2);

lat = char(m,n);

for i = 1:m
    for j = 1:n
        ind = index(i,j);
        if ind == 0
        else
            lat(i,j) = seq(ind);
        end
    end
end

objvalue=0;
tlat = lat;

for i = 1:m
    for j = 1:n
        if lat(i,j)=='H'
            if i~=m && lat(i+1,j)=='H'
                objvalue = objvalue - 1;
            end
            if j~= n && lat(i,j+1)=='H'
                objvalue = objvalue - 1;
            end
        end
    end
end


end