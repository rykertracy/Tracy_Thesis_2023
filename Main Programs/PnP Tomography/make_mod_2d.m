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
