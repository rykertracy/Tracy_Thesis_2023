function runallTTsizestcorrZZ(file, xcell, Vavg, LagS, dLag, LagR, tol, elim)
% elim tosses cells from map
% tol is .05 to toss data with 5% velocity higher than path with Vavg
% LagR is Lagrange multiplier for model over corrections
% dLag is multiplier from one Lagrange multiplier to the next
% LagS is the first Lagrange multiplier divided by dLag
% runallTTsizestcorrZZ('PN_MRG_RF',50,8.2,30,2.0,3000,.1,0.003)

global x y XV YV VV  % Declare global variables for use within the function

load(file)  % Load the input data from the specified file

% Extract x-coordinates (longitudes) of events and stations
xt = [cevent_x cstation_x]; 

% Extract y-coordinates (latitudes) of events and stations
yt = [cevent_y cstation_y]; 

% Calculate the distance between event and station for each pair
d = sqrt((xt(:,1) - xt(:,2)).^2 + (yt(:,1) - yt(:,2)).^2); 

% Find indices where the event-station distance is greater than 100 km
gt40 = find(d > 100);

% Extract the event and station coordinates for the filtered pairs
xs = xt(gt40, :);  
ys = yt(gt40, :);

% Extract the corresponding Moho travel times for these pairs
Mtt = Moho_tt_positiveONLY(gt40); 

% Shift the coordinate system (possibly to align with a reference point)
Xshift = -991.41; 
Yshift = -645.7;

figure(201)  % Create or switch to figure window 201 for potential plotting

% Create a 2D model grid with x and y dimensions defined by 1600 xcell and 1300 xcell
% Initialize velocity model VV with average velocity Vavg
[XV, YV, VV] = make_mod_2d(1600, xcell, 1300, xcell, Vavg);

% Get the size of the velocity model grid VV
[yc, xc] = size(VV);  

% Get the size (number of rows) of the xs array (event-station pairs)
[a, b] = size(xs);  

% Identify unique station locations within 35 units from the x and y coordinates of the station positions
[XST, YST] = stafinder(xs(:,2), ys(:,2), 35);

% Match the stations based on proximity within 35 units, returning the index of stations
IST = stamatcher(xs(:,2), ys(:,2), 35, XST, YST);

% Identify unique event locations within 35 units from the x and y coordinates of the event positions
[XEV, YEV] = stafinder(xs(:,1), ys(:,1), 35);

% Match the events based on proximity within 35 units, returning the index of events
IEV = stamatcher(xs(:,1), ys(:,1), 35, XEV, YEV);

% Determine the maximum station and event indices
NST = max(IST);  % Total number of unique stations
NEV = max(IEV);  % Total number of unique events

% Initialize the counter for valid ray paths
j = 0;  

% Loop over each event-station pair
for i = 1:a
    
    % Extract the x and y coordinates for the current event-station pair
    x = xs(i, :);
    y = ys(i, :);
    
    % Generate a 2D ray path between the event and station
    [XT, YT] = make_ray_2d(x, y, xcell / 4);
    
    % Calculate the travel time along the ray path and get the ray path values (NROWT)
    [TT, NROWT] = get_ray_time(XT, YT);
    
    % Calculate the travel-time residual (difference between observed and calculated travel times)
    DTJ = Mtt(i) - TT;
    
    % Check if the ray path contains NaN values (i.e., invalid paths)
    allnan = isnan(NROWT);  
    gotnan = find(allnan == 1);  % Find indices where NaN values exist
    
    % Only proceed if there are no NaN values and the residual is within the tolerance
    if isempty(gotnan) == 1 && abs(DTJ / TT) < tol
        
        % Initialize zero arrays for station and event identifiers
        STROW = zeros(1, NST);  
        STROW(IST(i)) = 1;  % Mark the corresponding station
        
        EVROW = zeros(1, NEV);  
        EVROW(IEV(i)) = 1;  % Mark the corresponding event
        
        % Increment the valid ray path counter
        j = j + 1;  
        
        % Store the ray path, station/event associations, and travel times
        G(j, :) = NROWT;  % Ray path values
        STM(j, :) = STROW;  % Station association matrix
        EVM(j, :) = EVROW;  % Event association matrix
        T(j) = Mtt(i);  % Observed travel time
        DT(j) = DTJ;  % Travel-time residual
    end
    % Uncomment the following lines for debugging or inspection
    % [sum(NROWT * 1 / 8) T(j) DT(j)]
    % pause
end

% Output the final values of the loop counters
i = i;  % Last index of the loop
j = j;  % Number of valid ray paths

[yc, xc] = size(YV);  % Get the size of the YV matrix (probably a grid of y-coordinates)

% Initialize the identity matrix D1 with dimensions based on the number of grid points (xc * yc)
D1 = eye(xc * yc);

% Modify D1 by shifting the last column to the first position, creating a finite difference operator in the x-direction
D1 = D1 - [D1(:, end) D1(:, 1:end-1)];

% Make D1 symmetric by adding its transpose, completing the finite difference in the x-direction
D1 = D1 + D1';

% Initialize another identity matrix D with the same dimensions (xc * yc)
D = eye(yc * xc);

% Modify D to create a finite difference operator in the y-direction by shifting columns by 'yc' (grid size in y-direction)
D = D - [D(:, (yc+1):end) D(:, 1:yc)];

% Further modify D to create diagonal finite differences by shifting the columns again
D = D - [D(:, end-(yc-1):end) D(:, 1:end-yc)];

% Make D symmetric by adding its transpose, combining the finite differences in x and y directions
D = D + D1;

% Get the size of the velocity model grid (YV)
[vyc, vxc] = size(YV);

% Compute the sum of ray path values (G) over all rows, storing the sum for each grid cell
SG = sum(G);

% Find grid cells where the sum is below a threshold defined by 'elim * max(SG)' (cells with insufficient ray coverage)
SG0 = find(SG < elim * max(SG));

% Set all ray path sums to 1 (this is likely initializing or resetting a mask)
SG = ones(size(SG));

% Set the elements in SG corresponding to the low-coverage cells (SG0) to NaN, marking them as invalid
SG(SG0) = 0 / 0;  % NaN assignment to low-coverage cells

% Reshape SG into the size of the velocity model grid (vyc x vxc) for further processing
SG = reshape(SG, vyc, vxc);

% Concatenate the ray path matrix (G), station matrix (STM), and event matrix (EVM) into a single matrix BM
BM = [G STM EVM];

% Get the size of YV again (the velocity model grid)
[yc, xc] = size(YV);

% Compute the normal equation matrix BMTBM = BM' * BM (used in least squares fitting)
BMTBM = BM' * BM;

% Initialize SZ0 as an identity matrix with dimensions matching BMTBM
SZ0 = eye(size(BMTBM));

% Scale the top left block of SZ0 (corresponding to the velocity model) by the Lagrange multiplier LagR
SZ0(1:vyc*vxc, 1:vyc*vxc) = LagR * eye(size(D));

% Assign SZ0 to SZ (this may be the matrix used for regularization in the inversion)
SZ = SZ0;

% Set the Lagrange multiplier to LagS (could be the starting value for regularization)
Lag = LagS;

% Apply coordinate shifts to the x and y velocity grids (XV and YV), likely for alignment purposes
XV = XV + Xshift;
YV = YV + Yshift;

for i = 1:10  % Loop to perform 10 iterations of inversion adjustment and visualization
    figure(i)  % Create a new figure for each iteration
    clf        % Clear the figure
    
    subplot(2,1,1)  % Create a subplot for displaying velocity model updates
    
    % Update the Lagrange multiplier by scaling with dLag
    Lag = Lag * dLag;
    
    % Solve the regularized least squares problem:
    % N is the solution vector containing updates to the velocity model and potentially station/event corrections
    N = inv(BMTBM + Lag * SZ) * BM' * DT';
    
    % Extract the velocity model updates (NV) from the solution vector N
    NV = N(1:vyc*vxc);
    
    % Convert NV to new velocity values using the starting average velocity (Vavg)
    NV = 1 ./ ((1/Vavg) + NV);
    
    % Compute the mean of the updated velocity model (MNV)
    MNV = mean(NV);
    
    % Normalize NV so that it represents percentage changes from the mean velocity
    NV = 100 * (1/MNV) * NV - 100;
    
    % Compute the mean and standard deviation of the normalized velocity values
    MNV2 = mean(NV);
    STDNV = std(NV);
    
    % Limit the standard deviation to 5 if it exceeds this value, possibly to avoid large outliers
    if STDNV > 5, STDNV = 5; end
    
    % Reshape NV into a 2D grid (matching the size of the velocity model)
    NP = reshape(NV, vyc, vxc);
    
    % Plot the updated velocity model with the ray coverage (SG) as a mask
    pcolor(XV, YV, SG .* NP)
    shading interp  % Interpolate colors between grid cells for smooth visualization
    colorbar        % Display a color scale bar
    grid on         % Show grid lines on the plot
    
    % Title showing the current Lagrange multiplier and average velocity
    title(['lagrange ' num2str(Lag) ' with avg vel of ' num2str(MNV)])
    
    % Set color axis limits based on the mean and standard deviation of the updated velocity model
    caxis([MNV2 - 1.5 * STDNV, MNV2 + 1.5 * STDNV])
    
    % Extract the station corrections from N
    NNST = N((vyc * vxc):(vyc * vxc + NST - 1))';
    
    % Interpolate the station corrections onto the velocity grid for visualization
    MNST = griddata(XST, YST, NNST, XV - Xshift, YV - Yshift);
    
    % Plot the station corrections on the second subplot
    subplot(2,1,2)
    pcolor(XV, YV, MNST)
    
    % Compute the mean and standard deviation of the station corrections
    MNNST = mean(NNST);
    STDNST = std(NNST);
    
    grid on           % Show grid lines
    shading interp    % Interpolate colors between grid cells
    caxis([MNNST - 1.5 * STDNST, MNNST + 1.5 * STDNST])  % Set color axis limits for station corrections
    colorbar          % Display a color scale bar for station corrections
    
    % Pause is commented out, but it would allow stepping through each iteration
    % pause
end

% Save important variables for future use or review
save('G', 'G', 'DT', 'T', 'N', 'xc', 'yc', 'XV', 'YV')
