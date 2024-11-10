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

% Save the results
save('G','G','DT','T','N','xc','yc','XV','YV')
