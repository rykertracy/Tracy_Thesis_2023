function runallTTsmthstcorrZZ(file,xcell,Vavg,LagS,dLag,LagR,tol,elim)
% elim - cells with low ray coverage below a threshold are excluded from the map
% tol - sets the threshold for discarding data with more than 5% velocity higher than the average velocity (Vavg)
% LagR - Lagrange multiplier controlling the smoothness of the model corrections
% dLag - factor used to update the Lagrange multiplier in each iteration
% LagS - initial Lagrange multiplier value
% Example call: runallTTsmthstcorrZZ('PN_MRG_RF',50,8.2,30,2.0,3000,.1,0.003)

global x y XV YV VV  % Declare global variables used throughout the program

load(file)  % Load the dataset from a file

% Extract coordinates of events and stations from loaded data
xt = [cevent_x cstation_x];
yt = [cevent_y cstation_y];

% Compute the distances between events and stations
d = sqrt((xt(:,1)-xt(:,2)).^2 + (yt(:,1)-yt(:,2)).^2);

% Select ray paths with distances greater than 100 units
gt40 = find(d > 100);

% Extract event and station coordinates for these paths
xs = xt(gt40,:);
ys = yt(gt40,:);
Mtt = Moho_tt_positiveONLY(gt40);  % Travel times corresponding to the selected paths

% Define grid shifts to align the velocity model with the data
Xshift = -1003.9;
Yshift = -658.2;

% Create a figure to visualize the data
figure(201)

% Create a 2D grid for the velocity model with cell size defined by xcell and initial average velocity (Vavg)
[XV, YV, VV] = make_mod_2d(1600, xcell, 1300, xcell, Vavg);

[yc, xc] = size(VV);  % Get the size of the velocity model grid
[a, b] = size(xs);    % Get the number of rays

% Identify and match stations and events with adjusted tolerance of 5
[XST, YST] = stafinder(xs(:,2), ys(:,2), 5);  % Find unique stations
IST = stamatcher(xs(:,2), ys(:,2), 35, XST, YST);  % Match stations to rays

[XEV, YEV] = stafinder(xs(:,1), ys(:,1), 5);  % Find unique events
IEV = stamatcher(xs(:,1), ys(:,1), 35, XEV, YEV);  % Match events to rays

% Maximum station and event indices
NST = max(IST);
NEV = max(IEV);

j = 0;  % Initialize a counter for valid rays
for i = 1:a  % Loop over all rays
    
    x = xs(i,:);  % Extract event and station x-coordinates
    y = ys(i,:);  % Extract event and station y-coordinates
    [XT, YT] = make_ray_2d(x, y, xcell / 4);  % Compute ray path coordinates on the grid
    
    [TT, NROWT] = get_ray_time(XT, YT);  % Compute travel times and ray path weights
    DTJ = Mtt(i) - TT;  % Compute the residual travel time (observed - calculated)
    
    % Check for NaN values in the ray path weights
    allnan = isnan(NROWT);
    gotnan = find(allnan == 1);
    
    % Include only valid rays (no NaN values) and whose residuals are within tolerance
    if isempty(gotnan) && abs(DTJ / TT) < tol
        STROW = zeros(1, NST);  % Initialize station row for this ray
        STROW(IST(i)) = 1;      % Mark the station corresponding to this ray
        EVROW = zeros(1, NEV);  % Initialize event row for this ray
        EVROW(IEV(i)) = 1;      % Mark the event corresponding to this ray
        j = j + 1;              % Increment the valid ray counter
        G(j,:) = NROWT;         % Store the ray path weights
        STM(j,:) = STROW;       % Store the station row
        EVM(j,:) = EVROW;       % Store the event row
        T(j) = Mtt(i);          % Store the observed travel time
        DT(j) = DTJ;            % Store the residual travel time
    end
end

[yc, xc] = size(YV);  % Get the size of the velocity grid
D1 = eye(xc*yc);  % Create identity matrix for regularization
D1 = D1 - [D1(:,end) D1(:,1:end-1)];  % Apply finite difference operator along x-direction
D1 = D1 + D1';
D = eye(yc*xc);  % Another identity matrix for regularization in the y-direction
D = D - [D(:,(yc+1):end) D(:,1:yc)];
D = D - [D(:,end-(yc-1):end) D(:,1:end-yc)];
D = D + D1;  % Combine the regularization operators

[vyc, vxc] = size(YV);  % Get the size of the velocity model grid

% Create a ray coverage mask and filter out cells with insufficient ray coverage
SG = sum(G);
SG0 = find(SG < elim * max(SG));
SG = ones(size(SG));
SG(SG0) = NaN;
SG = reshape(SG, vyc, vxc);  % Reshape coverage mask to grid size

% Assemble the full system matrix (G for ray paths, STM for stations, EVM for events)
BM = [G STM EVM];
[yc, xc] = size(YV);  % Get size of the velocity model
BMTBM = BM' * BM;     % Compute the data-fitting term in the least-squares system

% Set up the regularization term
SZ0 = eye(size(BMTBM));
SZ0(1:vyc*vxc, 1:vyc*vxc) = LagR * D;  % Apply Lagrange multiplier with grid regularization
SZ = SZ0;
Lag = LagS;  % Initialize Lagrange multiplier
XV = XV + Xshift;  % Apply grid shift in x-direction
YV = YV + Yshift;  % Apply grid shift in y-direction

for i = 1:10  % Loop through inversion iterations
    figure(i)  % Create a new figure for each iteration
    clf        % Clear the figure
    
    subplot(2,1,1)  % First subplot for the updated velocity model
    
    % Update Lagrange multiplier
    Lag = Lag * dLag;
    
    % Solve the least-squares inversion problem
    N = inv(BMTBM + Lag * SZ) * BM' * DT';
    
    % Extract and update the velocity values
    NV = N(1:vyc*vxc);
    NV = 1 ./ ((1 / Vavg) + NV);
    MNV = mean(NV);  % Compute mean velocity
    NV = 100 * (1 / MNV) * NV - 100;  % Normalize velocities as percent differences from mean
    MNV2 = mean(NV);  % Compute mean of normalized velocities
    STDNV = std(NV);  % Compute standard deviation of velocities
    if STDNV > 5, STDNV = 5; end  % Cap the standard deviation at 5%
    
    % Reshape normalized velocities into a grid and plot
    NP = reshape(NV, vyc, vxc);
    pcolor(XV, YV, SG .* NP)  % Plot the velocity model
    shading interp  % Use interpolated shading
    colorbar  % Add a colorbar
    grid on   % Show grid
    % axis equal  % Uncomment if needed to make axis equal
    
    title(['Lagrange multiplier: ' num2str(Lag) ', Average velocity: ' num2str(MNV)])
    caxis([MNV2 - 1.5 * STDNV, MNV2 + 1.5 * STDNV])  % Set color axis limits based on velocity distribution
    
    % Second subplot for station corrections
    NNST = N((vyc*vxc):(vyc*vxc+NST-1))';  % Extract station correction terms
    MNST = griddata(XST, YST, NNST, XV - Xshift, YV - Yshift);  % Interpolate onto the grid
    subplot(2,1,2)  % Second subplot
    pcolor(XV, YV, MNST)  % Plot the station corrections
    MNNST = mean(NNST);  % Compute mean station correction
    STDNST = std(NNST);  % Compute standard deviation of station corrections
    grid on  % Show grid
    shading interp  % Use interpolated shading
    caxis([MNNST - 1.5 * STDNST, MNNST + 1.5 * STDNST])  % Set color axis limits
    colorbar  % Add a colorbar
end

% Save relevant data to a .mat file
save('G', 'G', 'DT', 'T', 'N', 'xc', 'yc', 'XV', 'YV')
