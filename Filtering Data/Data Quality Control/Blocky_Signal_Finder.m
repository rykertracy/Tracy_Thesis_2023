% Step 1
directory_path = 'D:\Seismic Data\sorted_by_station\';

% Step 2
subdirectories = dir(directory_path);
subdirectories = subdirectories([subdirectories.isdir]);
subdirectories = subdirectories(~ismember({subdirectories.name},{'.','..'}));

% Step 3
blocky_files = {}; % initialize empty cell array to hold names of bad files
threshold = 0.5; % set a threshold for standard deviation
window_size = 50;
pinger = 0;
for i = 500:numel(subdirectories)
    subdirectory_path = fullfile(directory_path, subdirectories(i).name);

    % Step 4
    mat_files = dir(fullfile(subdirectory_path, '*.mat'));

    % Step 5
    for j = 1:numel(mat_files)
        mat_file_path = fullfile(subdirectory_path, mat_files(j).name);

        % Step 6
        load(mat_file_path);
        data = data.data - mean(data.data);
        unique_amps = unique(data);
        unique_amps = unique_amps(unique_amps ~= 0);
        num_unique_amps = length(unique_amps);
        num_nonzero_amps = length(find(data ~= 0));
        num_nonunique_amps = num_nonzero_amps - num_unique_amps;

        % Check if the waveform is potentially blocky
        if num_nonunique_amps / num_nonzero_amps > threshold
            pinger = pinger + 1

            % Analyze the waveform using windows
            num_windows = floor(length(data) / window_size);
            repeat_counts = zeros(num_windows, 1);

            for l = 1:num_windows
                window = data((l-1)*window_size + 1:l*window_size);
                repeat_counts(l) = length(find(diff(window) == 0)) + 1;
            end


            % Check if any windows contain a significant number of repeats
            if max(repeat_counts) / window_size > threshold
                % Waveform is likely blocky and should be removed
                disp([mat_file_path, 'Blocky'])
                blocky_files{end+1} = [mat_files(j).folder '\' mat_files(j).name]; % save name of bad file to cell array
            end
        end
    end
end


