function [xf, yf] = stafinder(x, y, tol)
% stafinder - Find unique station coordinates based on a distance tolerance.
% This function takes x and y coordinates of stations and a tolerance value.
% It returns the unique station coordinates that are separated by at least 'tol' distance.
%
% Inputs:
%   x   - Vector of x-coordinates of the stations
%   y   - Vector of y-coordinates of the stations
%   tol - Tolerance distance for determining uniqueness of stations
%
% Outputs:
%   xf  - Vector of unique x-coordinates of the stations
%   yf  - Vector of unique y-coordinates of the stations

test = 9;  % Initialize a variable for controlling the while loop
i = 0;     % Initialize counter for unique stations

while test == 9  % Loop until all unique stations are found
    i = i + 1;   % Increment the unique station counter
    xf(i) = x(1); % Assign the first remaining x-coordinate to the output
    yf(i) = y(1); % Assign the first remaining y-coordinate to the output
    
    % Compute distances from the newly found unique station to all other stations
    d = sqrt((xf(i) - x).^2 + (yf(i) - y).^2);
    
    % Find indices of stations that are farther away than the tolerance
    others = find(d > tol);
    
    % Update the lists of x and y coordinates to only include those that are farther away
    x = x(others);
    y = y(others);
    
    % If no more stations are left, exit the loop
    if isempty(y) == 1
        test = 3;  % Change test variable to exit condition
    end
end
