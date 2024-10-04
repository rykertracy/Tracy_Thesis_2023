% Get a list of all the subdirectories in the 'sorted_by_events' directory
subdirs = dir('sorted_by_events');
subdirs = subdirs([subdirs.isdir]); % Remove non-directories from the list
subdirs = subdirs(3:end); % Remove the '.' and '..' directories from the list
counter = 0;

% Loop through each subdirectory
for i = 1:length(subdirs)
    subdir = subdirs(i).name;
    subdir_path = fullfile('sorted_by_events', subdir);
    
    % Get a list of all the files in the current subdirectory
    files = dir(subdir_path);
    files = files(~[files.isdir]); % Remove directories from the list
    
    % Loop through each file in the current subdirectory
    for j = 1:length(files)
        file = files(j).name;
        
        % Check if the current file matches any of the filenames in the second column of 'files_to_delete'
        if any(strcmp(file, more_bad_files))
            % If it does, delete the file
            file_path = fullfile(subdir_path, file);
            delete(file_path);
            counter = counter + 1
        end
    end
end
