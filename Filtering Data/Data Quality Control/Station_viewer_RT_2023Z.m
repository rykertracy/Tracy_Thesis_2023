%% Instructions and Details
%The following program, written for the Ryker Tracy Thesis 2023, will
%access directores that have been sorted into a master directory called
%'sorted_by_stations' (see stationsorter_RT_2023Z). From that master
%directory, the program will access station directories and display 20
%waveforms on the screen for rapid determination of whether or not the
%whole station is bad quality. 

%Once 20 waveforms have been displayed, the user should press enter if
%everything looks good, press 1 if there are bad waveforms, but the entire
%station is not bad, and press 2 if the entire station is bad. Pressing 1
%sill store the station name in an array called 'someBad', and pressing 2
%will store the station name in an array called 'badQ' (bad quality).


%% Access directories and initialize variables.
% Define the master directory
master_directory = 'D:\Seismic Data\sorted_by_station';

% Get a list of all the subdirectories in the master directory
subdirs = dir(master_directory);
subdirs = subdirs([subdirs.isdir]);  % Only keep the subdirectories

badQ = {};  % Initialize badQ to an empty cell array
someBad = {};

%% Loop through all directores and plot up to 20 available station waveforms.
% Loop through each subdirectory
for i = 1:length(subdirs)
    subdir = subdirs(i).name;
    if strcmp(subdir, '.') || strcmp(subdir, '..')
        continue  % Skip the current and parent directory listings
    end
    
    % Set the initial figure position and size
    fig_pos = [-400 700 400 300];
    k=1;
    
    % Access the files in the current subdirectory
    sub_dir_path = fullfile(master_directory, subdir);
    file_list = dir(sub_dir_path);  % Get all files in the subdirectory
    for j = 1:length(file_list)
        file_name = file_list(j).name;
        if ~endsWith(file_name, '.mat')
            continue  % Skip non-.mat files
        end
        file_path = fullfile(sub_dir_path, file_name);
        % Load the data from the file and plot it here

        % Set the figure position and size
        fig_pos(1) = fig_pos(1) + 400;
        fig_pos(2) = fig_pos(2);
        if fig_pos(1) > 1600
            fig_pos(1) = 0;
            fig_pos(2) = fig_pos(2)-300;
        end
        figure('Position', fig_pos)
        load(file_path);
        figure(j-2)
        plot(data.data)
        title(subdir)
        k=k+1;
        if k > 20
            break
        end
    end

%% Ask the user if they would like to store the station.    
    % Ask the user if the current station is bad
    userpermission = input('Is this station bad? Press ENTER to proceed, type 1 for further review, type 2 to store bad station: ','s');
    if ~isempty(userpermission) && strcmp(userpermission, '2')
        badQ{end+1} = subdir;
    elseif ~isempty(userpermission) && strcmp(userpermission, '1')
        someBad{end+1} = subdir;
    end
    close all
end
unique_badQ = unique(badQ);