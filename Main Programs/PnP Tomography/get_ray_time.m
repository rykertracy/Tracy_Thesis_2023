function [T, NROW] = get_ray_time(X, Y)
% get_ray_time - Calculates travel times along rays from given coordinates.
%
% Inputs:
%   X - x-coordinates of the ray
%   Y - y-coordinates of the ray
%
% Outputs:
%   T - Total travel time for the ray
%   NROW - Array of accumulated contributions to the ray times for grid points

global x y XV YV VV  % Declare global variables

% Get the size of the meshgrid
[a, b] = size(XV);
% Reshape the meshgrid coordinates and velocity values into 1D arrays
xv = reshape(XV, 1, a * b);  % x-coordinates in a 1D array
yv = reshape(YV, 1, a * b);  % y-coordinates in a 1D array
vv = reshape(VV, 1, a * b);  % velocity values in a 1D array
nv = 1 ./ vv;  % Calculate inverse velocities (slowness)
NROW = zeros(1, a * b);  % Initialize NROW to store accumulated contributions

% Extend the X and Y arrays to include the starting points
X = [x(1) X x(2)];
Y = [y(1) Y y(2)];

% Calculate midpoints between points on the ray
XM = (X(1:end-1) + X(2:end)) / 2;  % Midpoint x-coordinates
YM = (Y(1:end-1) + Y(2:end)) / 2;  % Midpoint y-coordinates

% Calculate distances between consecutive points on the ray
LM = sqrt((X(1:end-1) - X(2:end)).^2 + (Y(1:end-1) - Y(2:end)).^2);

% Loop through each segment of the ray
for i = 1:length(LM)
    xm = XM(i);  % Current midpoint x-coordinate
    ym = YM(i);  % Current midpoint y-coordinate
    
    % Calculate the distance from the midpoint to all grid points
    td = sqrt((xm - xv).^2 + (ym - yv).^2);
    [ddist, id] = sort(td);  % Sort distances and get indices of sorted array
    
    % Get the inverse velocities of the three closest grid points
    n1 = nv(id(1)); n2 = nv(id(2)); 
    x1 = xv(id(1)) - xv(id(1));    y1 = yv(id(1)) - yv(id(1));  % Reference point
    x2 = xv(id(2)) - xv(id(1));    y2 = yv(id(2)) - yv(id(1));  % First neighbor
    x3 = xv(id(3)) - xv(id(1));    y3 = yv(id(3)) - yv(id(1));  % Second neighbor
    n3 = nv(id(3));  % Inverse velocity for the second neighbor
    xt = xm - xv(id(1));  % x-distance to the reference point
    yt = ym - yv(id(1));  % y-distance to the reference point
    kk = 3;  % Counter for the third neighbor

    % Find a valid third neighbor if the first three points are collinear
    while round(x1 + x2 + x3) == 0 || round(y1 + y2 + y3) == 0
        kk = kk + 1;
        x3 = xv(id(kk)) - xv(id(1));    % Update x3
        y3 = yv(id(kk)) - yv(id(1));    % Update y3
        n3 = nv(id(kk));  % Update inverse velocity for the third neighbor
    end
    
    % Create the system of equations to calculate the travel time
    G = [x1 y1 1; x2 y2 1; x3 y3 1];  % Geometry matrix
    N = [n1 n2 n3]';  % Inverse velocity vector
    W = inv(G' * G) * G';  % Weight matrix using the least squares solution
    M = W * N;  % Solve for the slowness model
    
    % Calculate travel time for the segment
    nm = [xt yt 1] * M;  % Interpolated slowness at the midpoint
    t(i) = LM(i) * nm;  % Travel time for this segment
    
    % Accumulate contributions to NROW for each grid point
    nr1 = W(1, 1) * xt + W(2, 1) * yt + W(3, 1);  % Contribution from first neighbor
    nr2 = W(1, 2) * xt + W(2, 2) * yt + W(3, 2);  % Contribution from second neighbor
    nr3 = W(1, 3) * xt + W(2, 3) * yt + W(3, 3);  % Contribution from third neighbor
    NROW(id(1)) = NROW(id(1)) + LM(i) * nr1;  % Update NROW for first neighbor
    NROW(id(2)) = NROW(id(2)) + LM(i) * nr2;  % Update NROW for second neighbor
    NROW(id(3)) = NROW(id(3)) + LM(i) * nr3;  % Update NROW for third neighbor
end

T = sum(t);  % Total travel time is the sum of all segment travel times
