% Step 1
directory_path = 'D:\Seismic Data\sorted_by_station\';

% Step 2
subdirectories = dir(directory_path);
subdirectories = subdirectories([subdirectories.isdir]);
subdirectories = subdirectories(~ismember({subdirectories.name},{'.','..'}));

% Step 3
bad_files = {}; % initialize empty cell array to hold names of bad files
std_threshold = 0.01; % set a threshold for standard deviation
for i = 1:numel(subdirectories)
    subdirectory_path = fullfile(directory_path, subdirectories(i).name);
    
    % Step 4
    mat_files = dir(fullfile(subdirectory_path, '*.mat'));
    
    % Step 5
    for j = 1:numel(mat_files)
        mat_file_path = fullfile(subdirectory_path, mat_files(j).name);
        
        % Step 6
        load(mat_file_path);
        
        % Step 7
        if std(data.data) < std_threshold
            % Step 8
            bad_files{end+1} = [mat_files(j).folder '\' mat_files(j).name]; % save name of bad file to cell array
        end
    end
end

% Display list of bad files
disp('The following files have majority linear or constant values:');
for i = 1:numel(bad_files)
    disp(bad_files{i});
end
