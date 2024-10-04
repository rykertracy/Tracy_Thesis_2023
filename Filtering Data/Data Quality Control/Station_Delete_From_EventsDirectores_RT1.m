%% Instructions and Details
%The following program, used for Ryker Tracy Thesis 2023, will use the .mat
%file called 'Bad Stations.mat' from the Station_viewer.m program. This
%program, as the name suggests, will delete bad stations from the events
%directory, whereas the other 'Station_Delete.m' program will delete the
%stations from the directory containing the files sorted by station.


%% Define the Master Directory and sort through each event folder and delete the files
% specify the path to the master directory

% load the list of bad files
load('files_to_delete.mat');
badFilenames = files_to_delete(:,2);

master_dir = 'D:\Seismic Data\sorted_by_events\';

% loop over each event folder in the master directory
event_dirs = dir(master_dir);
for ii = 1:length(event_dirs)
    if event_dirs(ii).isdir && ~strcmp(event_dirs(ii).name,'.') && ~strcmp(event_dirs(ii).name,'..')
        event_dir = fullfile(master_dir,event_dirs(ii).name);
        
        % loop over each file in the event folder
        event_files = dir(fullfile(event_dir,'*.mat'));
        
        for jj = 1:length(event_files)
            filename = event_files(jj).name;
            file_path = fullfile(event_dir, filename);
            
            % check if the file name contains any of the bad strings
            is_bad = ismember(filename,badFilenames);
            if is_bad
                % delete the bad file
                delete(file_path);
                fprintf('Deleted file: %s\n',file_path);
            end
        end
    end
end
