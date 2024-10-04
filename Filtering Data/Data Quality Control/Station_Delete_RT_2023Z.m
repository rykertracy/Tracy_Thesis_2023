%% Instructions and Details
%The following program, written for Ryker Tracy Thesis 2023, will delete
%stations that have been considered unusable from the 'Station_viewer.m'
%program. 'Station_viewer.m' stores an array of station names called 'Bad
%Stations'.

%% Load the file name strings, and delete the entire station directory
% Load the list of bad stations from the "BadQ.mat" file
master_dir = 'D:\Seismic Data\sorted_by_station\';

load('Bad Stations.mat');

% Loop through each bad station and delete its directory
for i = 1:length(badQ)
    bad_station = badQ{i};
    bad_station = char(bad_station);
    directory_path = char(fullfile(master_dir, bad_station));
    if isfolder(directory_path)
        rmdir(directory_path, 's');
        sprintf('Deleted directory for station %s.\n', bad_station);
    else
        sprintf('Directory for station %s does not exist.\n', bad_station);
    end
end
