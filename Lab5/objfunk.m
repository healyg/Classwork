function [tlength] = objfunk(lat,long,path)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This function takes in coordinates and path indices and outputs the total
% length of the path. 
% 
% INPUT
% 
% lat/long - the latitude and longitudinal coordinates of the 'cities'
% 
% path - the indices for the path you want to travel
% 
% OUTPUT 
% 
% tlength - the lenght, in units, of the chosen path. 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3
    error('You need three vectors, two for the coordinates and one for the indices')
elseif length(lat) ~= length(long) || length(lat)~=length(path)
    error('Your coordinates and indices must be the same length')
end

tlength = 0;
newlat = zeros(1,length(path));
newlong = zeros(1,length(path));

for i = 1:length(path)
    ind = path(i);
    newlat(i) = lat(ind);
    newlong(i) = long(ind);
end

for i = 2:length(newlat)
    tlength = tlength + sqrt((newlat(i) - newlat(i-1))^2 + (newlong(i) - ...
        newlong(i-1))^2);
end


end