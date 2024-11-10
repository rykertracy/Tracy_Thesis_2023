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

save('G','Q','G','DT','T','N','xc','yc','XV','YV')
