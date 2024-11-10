% function D_names=_data_map_RT_2023Z
%% Instructions and Details
%The following program, used for Ryker Tracy Thesis 2023, opens all
%directories that begin with the string 'Pn_DATA' (standing for Pn
%tomography data), gathers variable information, and plots all event
%locations on a lat lon grid, plots all stations on a lat lon grid, and
%plots the estimated plan-view ray paths on a lat lon grid.

%The idea of this program is to get a feel for data cluster locations and
%determine where more data may need to be acquired.

%NOTE: This program needs to be run after running the 'SACtoMATB_2023Z'
%function.

%% Loop through directories and store variables for plotting
%This section accesses all directores that begin with the string 'Pn_DATA'
%and finds the filenames in the further MAT folder
% files=[];
% Ddirs=dir('Pn_DATA*');
% k=0;
% for n=14:length(Ddirs)
%     files_t=dir([Ddirs(n).name '\MAT']);
%     for m=3:length(files_t);
%       k=k+1;
%       D_names(k).name=[pwd '\' Ddirs(n).name '\MAT\' files_t(m).name ];
%     end
% end

%Program can take as long as 140 minutes to run.
c = {};
directory = dir('sorted_by_events');
directory = directory(3:end);
load('Clean Signals From Events.mat')
figure(1)
hold on
for i = 1:length(directory)
    files = dir([directory(i).folder '\' directory(i).name]);
    files = files(3:end);
    for n=1:length(files)
        % file=D_names(n).name;
        c{n} = 'k';
        file = [files(n).folder '\' files(n).name];
        load(file);
        STlat(n)=data.station_lat;
        STlon(n)=data.station_lon;
        EVlat(n)=data.event_lat;
        EVlon(n)=data.event_lon;
        if STlat(n) < 20 || STlat(n) > 40 || STlon(n) < -115 || STlon(n) > -90 || EVlat(n) < 20 || EVlat(n) > 40 || EVlon(n) < -115 || EVlon(n) > -90
            delete(file)
            continue;
        end
        if ismember(file,clean_signal_all)
            c{n} = 'r';
        end
        plot([STlon(n) EVlon(n)],[STlat(n) EVlat(n)])
        hold on
    end

    title('Map Ray Paths for Events');
end
whos

%% Further Plots
%Station location plot
figure(2)
plot(STlon,STlat,'.')
title('Map of Station Locations');

disp(c);

%Event location plot with magnitude-dependent point-size
figure(3)
plot(EVlon,EVlat,[c, '.'])
title('Map of Event locations');

%Station and event location comparison
figure(4)
plot(STlon,STlat,'r.')
hold on
plot(EVlon,EVlat,'k.')

title('Map of Stations and Events');
% end

