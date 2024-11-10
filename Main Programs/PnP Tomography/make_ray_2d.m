function [X, Y] = make_ray_2d(x, y, dl)
% make_ray_2d - Generate a ray with points along its path in 2D space.
% This function takes two points (x, y) and creates evenly spaced points 
% along the line segment connecting them based on the specified distance 
% interval (dl).
%
% Inputs:
%   x  - A vector containing two x-coordinates [x1, x2]
%   y  - A vector containing two y-coordinates [y1, y2]
%   dl - Distance interval for spacing the points along the ray
%
% Outputs:
%   X  - Vector of x-coordinates of points along the ray
%   Y  - Vector of y-coordinates of points along the ray

% Calculate the length of the ray segment between points (x1, y1) and (x2, y2)
L = sqrt((x(2) - x(1))^2 + (y(2) - y(1))^2);
x1 = x(1);  % First x-coordinate
x2 = x(2);  % Second x-coordinate

% Adjust x2 slightly to avoid division by zero if x1 equals x2
if x1 == x2
    x2 = x2 + 0.01 * dl;  % Shift x2 by a small amount
end

% Calculate the slope (m) of the line segment
m = (y(2) - y(1)) / (x2 - x1);

% Determine the number of points (nl) along the ray based on the specified distance interval
nl = floor(L / dl);
dx = (x2 - x1) / nl;  % Calculate the increment in x for each point

% Initialize the arrays to hold the coordinates of the ray
% X(1) = x(1); 
% Y(1) = y(1);
for i = 1:nl-1
    % Calculate the x and y coordinates for each point along the ray
    X(i) = x(1) + i * dx;         % X-coordinate of the current point
    Y(i) = y(1) + i * dx * m;     % Y-coordinate of the current point based on the slope
end

% Plot the points along the ray in red
plot(X, Y, 'r*');   % Plot the ray points
hold on;
plot(x, y, 'r*');   % Plot the original points (start and end)
plot(X, Y, 'r*');    % Plot the ray points again (redundant line)

end
