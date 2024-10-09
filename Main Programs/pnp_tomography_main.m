function pnp_tomography_main(file, xcell, Vavg, LagS, dLag, LagR, tol, elim)
    % Display the menu options
    disp('Select an inversion program to run:');
    disp('1. runallTTsizestcorr.m');
    disp('2. runallTTsizestcorrZZ.m');
    disp('3. runallTTsmthstcorr.m');
    disp('4. runallTTsmthstcorrZZ.m');
    disp('5. Exit.');

    % Get user input
    choice = input('Enter your choice (1-5): ');

    % Process the user's choice
    switch choice
        case 1
            disp('Running runallTTsizestcorr...');
            runallTTsizestcorr(file, xcell, Vavg, LagS, dLag, LagR, tol, elim);  % Call your first inversion program function
        case 2
            disp('Running Inversion Program 2...');
            runallTTsizestcorrZZ(file, xcell, Vavg, LagS, dLag, LagR, tol, elim);  % Call your second inversion program function
        case 3
            disp('Running runallTTsmthstcorr...');
            runallTTsmthstcorr(file, xcell, Vavg, LagS, dLag, LagR, tol, elim);  % Call your third inversion program function
        case 4
            disp('Running runallTTsmthstcorrZZ.m...');
            runallTTsmthstcorrZZ(file, xcell, Vavg, LagS, dLag, LagR, tol, elim);
        case 5
            disp('Exiting...');
            return;  % Exit the main function
        otherwise
            disp('Invalid choice. Please select a valid option.');
            main();  % Restart the menu if the choice is invalid
    end
end

%% Inversion Function 1 (Size)
function runallTTsizestcorr(file, xcell, Vavg, LagS, dLag, LagR, tol, elim)
% Main function to run 2D velocity correction analysis using travel times
% Input parameters:
%   - file: the filename containing data
%   - xcell: cell size in x direction (in km)
%   - Vavg: average velocity (km/s)
%   - LagS: initial Lagrange multiplier for smoothing
%   - dLag: factor to adjust Lagrange multiplier in each iteration
%   - LagR: regularization term for the model
%   - tol: tolerance for excluding data with a large velocity discrepancy
%   - elim: threshold for removing low-sensitivity cells

global x y XV YV VV

% Load the data from the specified file
load(file)

% Adjust the x and y coordinates by shifting them so the minimum is centered
Xshift = min([cevent_x; cstation_x]) + xcell / 2;
xt = [cevent_x cstation_x] - Xshift;
Yshift = min([cevent_y; cstation_y]) + xcell / 2;
yt = [cevent_y cstation_y] - Yshift;

% Calculate distances between event-station pairs and find pairs greater than 50 km apart
d = sqrt((xt(:, 1) - xt(:, 2)).^2 + (yt(:, 1) - yt(:, 2)).^2);
gt40 = find(d > 50);  % Filter out pairs with distances under 50 km

% Extract the relevant x, y coordinates and Moho travel times for the selected pairs
xs = xt(gt40, :);
ys = yt(gt40, :);
Mtt = Moho_tt_positiveONLY(gt40);  % Positive Moho travel times only
IEV = eventnum(gt40);  % Event numbers for the selected pairs
NEV = max(IEV);  % Maximum event number

% Plotting
figure(201)

% Generate the 2D velocity model grid with the specified dimensions and average velocity
[XV, YV, VV] = make_mod_2d(1800, xcell, 1300, xcell, Vavg);

% Get grid sizes
[yc, xc] = size(VV);
[a, b] = size(xs);

% Find stations based on their x, y positions and match station IDs
[XST, YST] = stafinder(xs(:, 2), ys(:, 2), 35);
IST = stamatcher(xs(:, 2), ys(:, 2), 35, XST, YST);
NST = max(IST);  % Number of unique stations

j = 0;  % Initialize counter for valid ray paths

% Loop over all the selected ray paths
for i = 1:a
    x = xs(i, :);
    y = ys(i, :);

    % Create ray segments between stations and events
    [XT, YT] = make_ray_2d(x, y, xcell / 4);

    % Calculate travel time along the ray using the velocity model
    [TT, NROWT] = get_ray_time(XT, YT);
    DTJ = Mtt(i) - TT;  % Calculate the residual between observed and calculated travel time

    % Check if any NaNs exist in NROWT and if the residual is within tolerance
    allnan = isnan(NROWT);
    gotnan = find(allnan == 1);
    if isempty(gotnan) && abs(DTJ / TT) < tol
        % Build the system of equations to solve
        STROW = zeros(1, NST);
        STROW(IST(i)) = 1;  % Associate ray with station
        EVROW = zeros(1, NEV);
        EVROW(IEV(i)) = 1;  % Associate ray with event
        j = j + 1;
        G(j, :) = NROWT;  % Store sensitivity for each grid cell
        STM(j, :) = STROW;  % Store station information
        EVM(j, :) = EVROW;  % Store event information
        T(j) = Mtt(i);  % Store travel time
        DT(j) = DTJ;  % Store residual travel time
    end
end

i = i;
j = j;

% Compute the difference operator matrix for smoothing (used in regularization)
D1 = eye(xc * yc);
D1 = D1 - [D1(:, end) D1(:, 1:end-1)];
D1 = D1 + D1';
D = eye(yc * xc);
D = D - [D(:, (yc+1):end) D(:, 1:yc)];
D = D - [D(:, end-(yc-1):end) D(:, 1:end-yc)];
D = D + D1;

% Initialize smoothing mask based on grid sensitivity
SG = sum(G);
SG0 = find(SG < elim * max(SG));  % Cells with low sensitivity
SG = ones(size(SG));
SG(SG0) = NaN;  % Set low-sensitivity cells to NaN
SG = reshape(SG, yc, xc);

% Build the matrix for the inversion
BM = [G STM EVM];
BMTBM = BM' * BM;
SZ0 = eye(size(BMTBM));  % Regularization term
SZ0(1:yc*xc, 1:yc*xc) = LagR * eye(size(D));
SZ = SZ0;
Lag = LagS;

% Adjust X and Y to account for previous shifts
XV = XV + Xshift;
YV = YV + Yshift;

% Iterative inversion process
for i = 1:10
    figure(i)
    clf
    Lag = Lag * dLag;  % Update the Lagrange multiplier
    N = inv(BMTBM + Lag * SZ) * BM' * DT';  % Solve the system of equations
    NV = N(1:yc * xc);
    NV = 1 ./ ((1 / Vavg) + NV);  % Convert velocity perturbation to absolute velocity
    MNV = mean(NV);
    NV = 100 * (1 / MNV) * NV - 100;  % Normalize perturbation as percentage
    MNV2 = mean(NV);
    STDNV = std(NV);
    if STDNV > 5, STDNV = 5; end  % Cap the standard deviation
    NP = reshape(NV, yc, xc);  % Reshape the velocity perturbation to a 2D grid
    Q = SG .* NP;  % Apply the smoothing mask

    % Convert local coordinates to latitude and longitude
    [lat, lon] = local2latlon((XV * 1000), (YV * 1000), 0, [31, -101, 0]);
    pcolor(lon, lat, Q)  % Plot the velocity perturbation
    shading interp
    texas = shaperead('usastatehi', 'UseGeoCoords', true, 'Selector', {@(name) strcmpi(name, 'Texas'), 'Name'});
    ok = shaperead('usastatehi', 'UseGeoCoords', true, 'Selector', {@(name) strcmpi(name, 'Oklahoma'), 'Name'});
    NewMex = shaperead('usastatehi', 'UseGeoCoords', true, 'Selector', {@(name) strcmpi(name, 'New Mexico'), 'Name'});
    geoshow(texas, 'FaceColor', 'none')
    geoshow(ok, 'FaceColor', 'none')
    geoshow(NewMex, 'FaceColor', 'none')
    colorbar
    axis equal
    title('Vp Variation Along Moho | 100 Km Cells', 'FontSize', 15)
    subtitle(['Lagrange ' num2str(Lag) ' With Average Velocity of ' num2str(MNV), ' | Model B'])
    clim([MNV2 - 1.5 * STDNV, MNV2 + 1.5 * STDNV])
    xlabel('Longitude (Degrees)')
    ylabel('Latitude (Degrees)')
    ylabel(colorbar, 'Percent Variation from Average', 'FontSize', 13, 'Rotation', 270)

    % Pause between iterations for visual inspection
end

% Save the results for future use
save('G', 'Q', 'G', 'DT', 'T', 'N', 'xc', 'yc', 'XV', 'YV')

Mtt=Moho_tt_positiveONLY(gt40);
IEV=eventnum(gt40);
NEV=max(IEV);
figure(201)

[XV,YV,VV]=make_mod_2d(1800,xcell,1300,xcell,Vavg);

[yc,xc]=size(VV)
[a,b]=size(xs)
[XST,YST]=stafinder(xs(:,2),ys(:,2),35);
IST=stamatcher(xs(:,2),ys(:,2),35,XST,YST);

NST=max(IST);

j=0
for i=1:a

    x=xs(i,:);
    y=ys(i,:);
    [XT,YT]=make_ray_2d(x,y,xcell/4);

    [TT,NROWT]=get_ray_time(XT,YT);
    DTJ=Mtt(i)-TT;

    allnan=isnan(NROWT);
    gotnan=find(allnan==1);
    if isempty(gotnan)==1 & abs(DTJ/TT) < tol,
        STROW=zeros(1,NST);
        STROW(IST(i))=1;
        EVROW=zeros(1,NEV);
        EVROW(IEV(i))=1;
        j=j+1;
        G(j,:)=NROWT;
        STM(j,:)=STROW;
        EVM(j,:)=EVROW;
        T(j)=Mtt(i);
        DT(j)=DTJ;
    end
    %     [sum(NROWT*1/8) T(j) DT(j)]=
    %     pause
end
i=i
j=j

[yc,xc]=size(YV);
D1=eye(xc*yc);
D1=D1-[D1(:,end) D1(:,1:end-1)];
D1=D1+D1';
D=eye(yc*xc);
D=D-[D(:,(yc+1):end) D(:,1:yc)];
D=D-[D(:,end-(yc-1):end) D(:,1:end-yc)];
D=D+D1;

[vyc,vxc]=size(YV);

SG=sum(G);
SG0=find(SG<elim*max(SG));
SG=ones(size(SG));
SG(SG0)=0/0;
SG=SG;
SG=reshape(SG,vyc,vxc);

BM=[G STM EVM];
[yc,xc]=size(YV);
BMTBM=BM'*BM;
SZ0=eye(size(BMTBM));
SZ0(1:vyc*vxc,1:vyc*vxc)=LagR*eye(size(D));
SZ=SZ0;
Lag=LagS
XV=XV+Xshift;
YV=YV+Yshift;

for i=1:10
    figure(i)
    clf
    % subplot(2,1,1)
    Lag=Lag*dLag;
    N=inv(BMTBM+Lag*SZ)*BM'*DT';
    NV=N(1:vyc*vxc);
    NV=1./((1/Vavg)+NV);
    MNV=mean(NV);
    NV=100*(1/MNV)*NV-100;
    MNV2=mean(NV);
    STDNV=std(NV); if STDNV > 5, STDNV=5; end
    NP=reshape(NV,vyc,vxc);
    Q = SG.*NP;
    [lat,lon]=local2latlon((XV.*1000),(YV.*1000),0,[31,-101,0]);
    pcolor(lon,lat,Q)
    shading interp
    texas = shaperead('usastatehi', 'UseGeoCoords', true,'Selector',{@(name) strcmpi(name,'Texas'), 'Name'});
    ok = shaperead('usastatehi', 'UseGeoCoords', true,'Selector',{@(name) strcmpi(name,'Oklahoma'), 'Name'});
    NewMex = shaperead('usastatehi', 'UseGeoCoords', true,'Selector',{@(name) strcmpi(name,'New Mexico'), 'Name'});
    geoshow(texas,'FaceColor','none')
    geoshow(ok,'FaceColor','none')
    geoshow(NewMex,'FaceColor','none')
    colorbar
    axis equal
    title('Vp Variation Along Moho | 100 Km Cells','FontSize',15)
    subtitle(['Lagrange ' num2str(Lag) ' With Average Velocity of ' num2str(MNV), ' | Model B'])
    clim([MNV2-1.5*STDNV MNV2+1.5*STDNV])
    xlabel('Longitude (Degrees)')
    ylabel('Latitude (Degrees)')
    ylabel(colorbar,'Percent Variation from Average','FontSize',13,'Rotation',270)

    % NNST=N((vyc*vxc):(vyc*vxc+NST-1))';
    % MNST=griddata(XST,YST,NNST,XV-Xshift,YV-Yshift);
    % % subplot(2,1,2)
    % pcolor(lon,lat,MNST)
    % geoshow(texas,'FaceColor','none')
    % geoshow(ok,'FaceColor','none')
    % geoshow(NewMex,'FaceColor','none')
    % MNNST=mean(NNST);
    % STDNST=std(NNST);
    % title('Vp Variation Along Moho','FontSize',15)
    % subtitle(['Lagrange ' num2str(Lag) ' With Average Velocity of ' num2str(MNV), ' | Model B | Crustal Correction'])
    % % axis equal
    % xlabel('Longitude (Degrees)')
    % ylabel('Latitude (Degrees)')
    % caxis([MNNST-1.5*STDNST MNNST+1.5*STDNST])
    % colorbar


    %     NNST=N((vyc*vxc+NST-1):end)'
    %     pause
end
% Save relevant data to a .mat file
save('G', 'G', 'DT', 'T', 'N', 'xc', 'yc', 'XV', 'YV')
end

%% Inversion Function 2 (Size ZZ)
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
% Save important variables for future use or review
end
% Save relevant data to a .mat file
save('G', 'G', 'DT', 'T', 'N', 'xc', 'yc', 'XV', 'YV')
end

%% Inversion Function 3 (Smooth)
function runallTTsmthstcorr(file,xcell,Vavg,LagS,dLag,LagR,tol,elim)
% elim: threshold to toss cells from the map with poor sensitivity
% tol: threshold to toss data with velocities 5% higher than the average velocity (Vavg)
% LagS: initial Lagrange multiplier for regularization
% dLag: factor to increase Lagrange multiplier at each iteration
% LagR: Lagrange multiplier for smoothness constraints on the model
% file: the name of the input data file
% xcell: grid cell size
% Vavg: average velocity for the initial model
% runallTTsmthstcorr('PN_MRG_1M', 50, 8.2, 0.1, 3.0, 3000, .1, 0.003)

global x y XV YV VV

% Load the data from the file
load(file)

% Shift the X and Y coordinates of the events and stations for grid alignment
Xshift = min([cevent_x; cstation_x]) + xcell/2;
xt = [cevent_x cstation_x] - Xshift;
Yshift = min([cevent_y; cstation_y]) + xcell/2;
yt = [cevent_y cstation_y] - Yshift;

% Calculate the distances between event and station pairs
d = sqrt((xt(:,1)-xt(:,2)).^2 + (yt(:,1)-yt(:,2)).^2);

% Find ray paths where the distance between event and station is greater than 100 km
gt40 = find(d > 100);
xs = xt(gt40, :);
ys = yt(gt40, :);
Mtt = Moho_tt_positiveONLY(gt40); % Positive ray travel times along the Moho
IEV = eventnum(gt40); % Event numbers for the rays
NEV = max(IEV); % Total number of events

figure(201)

% Create a 2D model grid with initial average velocity (Vavg)
[XV, YV, VV] = make_mod_2d(1800, xcell, 1300, xcell, Vavg);
[yc, xc] = size(VV);
[a, b] = size(xs);

% Find station coordinates and match them with indices
[XST, YST] = stafinder(xs(:,2), ys(:,2), 35);
IST = stamatcher(xs(:,2), ys(:,2), 35, XST, YST);
NST = max(IST); % Total number of stations

% Initialize data storage
j = 0;
for i = 1:a
    x = xs(i, :); % Current ray's X coordinates
    y = ys(i, :); % Current ray's Y coordinates

    % Generate ray path for the current ray
    [XT, YT] = make_ray_2d(x, y, xcell/4);

    % Get the calculated ray travel time (TT) and the number of ray segments (NROWT)
    [TT, NROWT] = get_ray_time(XT, YT);

    % Compute the travel time residual (observed minus calculated)
    DTJ = Mtt(i) - TT;

    % Check if the ray path contains any NaN values and if the residual is within tolerance
    allnan = isnan(NROWT);
    gotnan = find(allnan == 1);
    if isempty(gotnan) && abs(DTJ/TT) < tol
        STROW = zeros(1, NST); % Station row for this ray
        STROW(IST(i)) = 1; % Assign current station

        EVROW = zeros(1, NEV); % Event row for this ray
        EVROW(IEV(i)) = 1; % Assign current event

        % Store the relevant ray data
        j = j + 1;
        G(j, :) = NROWT; % Raypath matrix
        STM(j, :) = STROW; % Station matrix
        EVM(j, :) = EVROW; % Event matrix
        T(j) = Mtt(i); % Travel times
        DT(j) = DTJ; % Travel time residuals
    end
end

% Set up smoothing regularization matrices
[yc, xc] = size(YV);
D1 = eye(xc * yc);
D1 = D1 - [D1(:,end) D1(:,1:end-1)];
D1 = D1 + D1';
D = eye(yc * xc);
D = D - [D(:,(yc+1):end) D(:,1:yc)];
D = D - [D(:,end-(yc-1):end) D(:,1:end-yc)];
D = D + D1;

% Compute sensitivity for the cells
[vyc, vxc] = size(YV);
SG = sum(G);
SG0 = find(SG < elim * max(SG)); % Cells with low sensitivity
SG = ones(size(SG));
SG(SG0) = NaN; % Mark low-sensitivity cells as NaN
SG = reshape(SG, vyc, vxc); % Reshape to the model grid

% Build data matrix and apply regularization
BM = [G STM EVM];
BMTBM = BM' * BM; % Matrix product for least-squares inversion
SZ0 = eye(size(BMTBM));
SZ0(1:vyc*vxc, 1:vyc*vxc) = LagR * D; % Apply smoothness regularization
SZ = SZ0;
Lag = LagS; % Set initial Lagrange multiplier

% Adjust model coordinates to match the shift applied earlier
XV = XV + Xshift;
YV = YV + Yshift;

% Plot the model results
amapc = colormap;
colormap(amapc);
for i = 1:10
    figure(i)
    clf
    subplot(2,1,1)
    Lag = Lag * dLag; % Increment Lagrange multiplier
    N = inv(BMTBM + Lag * SZ) * BM' * T'; % Perform inversion to solve for model parameters
    NV = N(1:vyc*vxc);
    NV = 1./(NV); % Compute velocities
    MNV = mean(NV);
    NV = 100 * (1/MNV) * NV - 100; % Calculate percent deviation from average velocity
    MNV2 = mean(NV);
    STDNV = std(NV); if STDNV > 5, STDNV = 5; end % Cap standard deviation at 5
    NP = reshape(NV, vyc, vxc); % Reshape model result to grid

    % Apply sensitivity mask and plot the model
    Q = SG .* NP;
    [lat, lon] = local2latlon((XV.*1000), (YV.*1000), 0, [31,-101,0]);
    pcolor(lon, lat, Q)
    colormap(amapc);
    shading interp
    % Display state boundaries
    texas = shaperead('usastatehi', 'UseGeoCoords', true, 'Selector', {@(name) strcmpi(name,'Texas'), 'Name'});
    ok = shaperead('usastatehi', 'UseGeoCoords', true, 'Selector', {@(name) strcmpi(name,'Oklahoma'), 'Name'});
    NewMex = shaperead('usastatehi', 'UseGeoCoords', true, 'Selector', {@(name) strcmpi(name,'New Mexico'), 'Name'});
    geoshow(texas,'FaceColor','none')
    geoshow(ok,'FaceColor','none')
    geoshow(NewMex,'FaceColor','none')
    colorbar
    title('Vp Variation Along Moho | 100 Km Cell - Smoothed','FontSize',16)
    subtitle(['Lagrange ' num2str(Lag) ' With Average Velocity of ' num2str(MNV), ' | Model B | Smoothed'])
    clim([MNV2-1.5*STDNV MNV2+1.5*STDNV])
    xlabel('Longitude (Degrees)','FontSize',13)
    ylabel('Latitude (Degrees)','FontSize',13)
    ylabel(colorbar,'Percent Variation from Average','FontSize',13,'Rotation',270)

    % Plot station variations in the model
    NNST = N((vyc*vxc):(vyc*vxc+NST-1))';
    MNST = griddata(XST, YST, NNST, XV-Xshift, YV-Yshift);
    subplot(2,1,2)
    pcolor(lon, lat, MNST)
    geoshow(texas,'FaceColor','none')
    geoshow(ok,'FaceColor','none')
    geoshow(NewMex,'FaceColor','none')
    MNNST = mean(NNST);
    STDNST = std(NNST);
    title('Vp Variation Along Moho','FontSize',15)
    subtitle(['Lagrange ' num2str(Lag) ' With Average Velocity of ' num2str(MNV), ' | Model B | Crustal Correction'])
    clim([MNV2-1.5*STDNV MNV2+1.5*STDNV])
    xlabel('Longitude (Degrees)')
    ylabel('Latitude (Degrees)')
    ylabel(colorbar,'Percent Variation from Average','FontSize',13,'Rotation',270)
    caxis([MNNST-2*STDNST MNNST+2*STDNST])
    colorbar

    pause
end
% Save relevant data to a .mat file
save('G', 'G', 'DT', 'T', 'N', 'xc', 'yc', 'XV', 'YV')
end

%% Inversion Functino 4 (Smooth ZZ)
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
end

%% get_ray_time
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
end

%% make_mod_2D
function [XV, YV, VV] = make_mod_2d(xmax, dx, ymax, dy, vbar)
% make_mod_2d - Create a 2D model grid for visualization and manipulation.
%
% Inputs:
%   xmax - Maximum x-coordinate for the grid
%   dx   - Grid spacing in the x-direction
%   ymax - Maximum y-coordinate for the grid
%   dy   - Grid spacing in the y-direction
%   vbar - Base velocity to be added to the model values
%
% Outputs:
%   XV   - Meshgrid of x-coordinates
%   YV   - Meshgrid of y-coordinates
%   VV   - 2D array representing model values

% Calculate the number of points in the y-direction (ny) and x-direction (nx)
ny = round(ymax/dy) + 1;  % Number of y points
nx = round(xmax/dx) + 1;  % Number of x points

% Initialize the model grid with zeros
VV = 0 * ones(ny, nx);  % Create a ny by nx matrix filled with zeros

% Create vectors for x and y coordinates based on the specified grid spacing
xv = [0:dx:(nx-1)*dx];  % x-coordinates from 0 to xmax
yv = [0:dy:(ny-1)*dy];  % y-coordinates from 0 to ymax

% Generate a meshgrid from the x and y coordinates for plotting
[XV, YV] = meshgrid(xv, yv);

% Display the current workspace variables (for debugging)
whos 

% Plot the initial model grid with the base velocity added
pcolor(XV, YV, VV + vbar);
colorbar;  % Add a color bar to the plot

% Define a smoothing kernel for model smoothing
smoother = 1/9 * ones(3, 3);  % Averaging filter of size 3x3
nsmooth = 10;  % Number of smoothing iterations

% Initialize control variables
do_what = 1;  % Control variable for user input
x = 1;  % Dummy variable for controlling input
ipick = 0;  % Counter for picked anomalies

% Main loop to handle user interaction
while do_what < 2.5
    do_what = 3;  % Reset to a default value for menu options

    % Uncomment to display a menu (commented out for demo)
    % do_what = menu('do what','add anomalies to model','smooth model','exit');

    % Case for smoothing the model
    if do_what == 2
        for ismooth = 1:nsmooth
            VV = conv2(VV, smoother, 'same');  % Apply smoothing
            for i = 1:ipick
                VV(iz(i), ix(i)) = amp(i);  % Re-add anomalies after smoothing
            end
        end
        pcolor(xv, yv, VV + vbar);  % Plot the smoothed model
        colorbar;
    end

    % Case for adding anomalies to the model
    if do_what == 1
        RVN = round([-50:50]) / 100 + 0;  % Define perturbation values
        RVA = hgnum2char(RVN, 2);  % Convert to character array for menu display
        IRVT = menu('Vp perturbation', RVA);  % User selects perturbation
        
        % Loop for picking locations on the plot
        while x > 0
            [x, z] = ginput(1);  % Get user input for coordinates
            if x > 0  % Ensure valid input
                ipick = ipick + 1;  % Increment picked counter
                ix(ipick) = round(x/dx);  % x index for picked anomaly
                iz(ipick) = round(z/dy);  % y index for picked anomaly
                amp(ipick) = RVN(IRVT);  % Assign perturbation value
                VV(iz(ipick), ix(ipick)) = amp(ipick);  % Update the model with the anomaly
                whos  % Display current workspace variables (for debugging)
                pcolor(xv, yv, VV + vbar);  % Plot the updated model
                colorbar;  % Add color bar
            end
        end
    end
    x = 1;  % Reset x for the next iteration
end

% Add the base velocity to the final model values
VV = VV + vbar;
hold on;  % Keep the current plot for further modifications
end

%% make_ray_2D
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

%% Stamatcher
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
end

%% Stafinder
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
end
