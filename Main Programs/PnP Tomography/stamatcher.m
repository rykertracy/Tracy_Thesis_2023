function imatch = stamatcher(x, y, tol, XR, YR)
% stafinder - Match stations to reference coordinates based on minimum distance.
% This function finds the index of the closest reference coordinates (XR, YR) 
% for each given station coordinate (x, y).
%
% Inputs:
%   x   - Vector of x-coordinates of the stations
%   y   - Vector of y-coordinates of the stations
%   tol - Tolerance distance (not currently used in the function)
%   XR  - Vector of x-coordinates of reference points
%   YR  - Vector of y-coordinates of reference points
%
% Outputs:
%   imatch - Vector of indices corresponding to the closest reference coordinates for each station

% Loop through each station coordinate
for i = 1:length(x)
    % Calculate the Euclidean distance from the current station (x(i), y(i)) to all reference points (XR, YR)
    d = sqrt((XR - x(i)).^2 + (YR - y(i)).^2);
    
    % Find the index of the closest reference coordinate and the distance to it
    [junk, imatcht] = min(abs(d));
    
    % Store the index of the closest reference coordinate for the current station
    imatch(i) = imatcht;
end
