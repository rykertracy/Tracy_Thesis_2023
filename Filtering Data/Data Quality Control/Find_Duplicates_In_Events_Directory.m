%% Introduction
%This program, used in Ryker Tracy 2023 Thesis, searches the
%'sorted_by_events' subdirectories (i.e. the events) for duplicate files in
%case some files were placed in multiple directories. It creates a growing
%list of filenames as it searches the next directories and determines if
%they are identical strings to any in the growing list of file names.

%% Program
% Set the path to the directory containing the subdirectories
path_to_sorted_by_events = 'D:\Seismic Data\sorted_by_events\';

% Get a list of all subdirectories in the main directory
subdirs = dir(path_to_sorted_by_events);
subdirs = subdirs([subdirs(:).isdir] & ~ismember({subdirs(:).name},{'.','..'}));

% Initialize an empty cell array to store filenames
filenames = {};

% Loop over each subdirectory
for ii = 1:length(subdirs)
    
    % Get a list of all files in the subdirectory
    subdir_path = fullfile(path_to_sorted_by_events, subdirs(ii).name);
    files = dir(fullfile(subdir_path, '*'));
    files = files(~[files(:).isdir]);
    
    % Loop over each file in the subdirectory
    for jj = 1:length(files)
        
        % Get the filename and check if it's already in the list
        filename = files(jj).name;
        if any(strcmp(filename, filenames))
            fprintf('Duplicate file: %s\n', fullfile(subdir_path, filename));
        else
            % Add the filename to the list
            filenames = [filenames, filename];
        end
        
    end
    
end
