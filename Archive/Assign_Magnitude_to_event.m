%function assign_magnitude(events_path)
% Define the path to the directory containing the subdirectories
directory_path = events_path;

% Get a list of all subdirectories in the directory
subdirs = dir(directory_path);
subdirs = subdirs([subdirs(:).isdir]); % only keep directories
subdirs = subdirs(~ismember({subdirs(:).name},{'.','..'})); % remove . and ..

% Loop through the subdirectories
for i = 1:length(subdirs)
    subdir_path = fullfile(directory_path, subdirs(i).name);
    
    % Load the first file in the subdirectory
    files = dir(fullfile(subdir_path, '*.mat'));
    data = load(fullfile(subdir_path, files(1).name));
    
    % Extract the magnitude
    magnitude = data.data.magnitude;
    
    % Create a new subdirectory with the magnitude in the name
    new_subdir_name = sprintf('%g_%s', magnitude, subdirs(i).name);
    new_subdir_path = fullfile(directory_path, new_subdir_name);
    mkdir(new_subdir_path);
    
    % Move all files in the subdirectory to the new subdirectory
    for j = 1:length(files)
        file_path = fullfile(subdir_path, files(j).name);
        movefile(file_path, new_subdir_path);
    end
end
end