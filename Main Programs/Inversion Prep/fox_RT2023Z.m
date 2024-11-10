function [cstation_x,cstation_y,cevent_x,cevent_y]=fox_RT2023Z(file,origin)
%% Introduction
% This program will convert latitudes and longitudes for surface
% event/station and moho-piercing points into local X-Y coordinates that are
% used in the Pn-inverstion program (i.e. runallTTsmthstcorr.m).
%
%
% The file is going to be some variation of the 'stuff_needed.mat' file. The
% origin is the center of the study area.
%
% The 'stuff_needed.mat' file needs to contain at least the following
% variables:
% Moho piercing points: cStlat, cStlon, cEQlat, cEQlon (in meters, or alter code below)
% surface points: c_stlat, c_stlon, c_evlat, c_evlon (in meters, or alter
% code below)
%
% The output of fox_RT2023Z will be cstation_x,cstaion_y,cevent_x,cevent_y.
% These are simple, arbitrary positive and negative X-Y coordinates with
% (0,0) at the center.
%
% Example execution for New Mexico, Oklahoma, and Texas study area (Tracy
% Thesis): fox_RT2023Z('stuff_needed2.mat',[31,-101,0]);
%

%% Program
% Load variables
load(file)

% Convert from lat/lon to "local" X-Y coordinate system
for n = 1:length(cStlat)
    [cstation_x(n), cstation_y(n)]=latlon2local(cStlat(n),cStlon(n),0,origin);
    [cevent_x(n), cevent_y(n)]=latlon2local(cEQlat(n),cEQlon(n),0,origin);


    % [station_x(n), station_y(n)]=latlon2local(c_stlat(n),c_stlon(n),0,origin);
    % [event_x(n), event_y(n)]=latlon2local(c_evlat(n),c_evlon(n),0,origin);
end

% Convert from meters to kilometers
cstation_x = cstation_x ./ 1000;
cstation_y = cstation_y ./ 1000;
cevent_x = cevent_x ./ 1000;
cevent_y = cevent_y ./ 1000;

%Make the data a single column for the inversion program
cevent_x = cevent_x';
cevent_y = cevent_y';
cstation_y = cstation_y';
cstation_x = cstation_x';

% Plot Events:
 % plot(event_x,event_y,'*','MarkerEdgeColor','r','MarkerFaceColor','r')
% Plot Stations:
 % plot(station_x,station_y,'*','MarkerEdgeColor','r','MarkerFaceColor','r')

% Save the new data
save('moho_xy.mat','cevent_x','cevent_y','cstation_x','cstation_y')


end