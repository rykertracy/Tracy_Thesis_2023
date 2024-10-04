%% Introduction
%This program will delete files if they are not the correct sample lenth to
%reach the intended wave length. This program was used to ensure there were
%no errors in the sample number using dt for a 60 second wave. I have
%allowed a tolerance of 1 second in case there are discrepencies that don't
%impact P-arrival time very much.

%% Program
directory_path = 'D:/Seismic Data/sorted_by_events/';

% Get a list of all subdirectories in the "sorted_by_station" directory
subdirectories = dir(directory_path);
subdirectories = subdirectories([subdirectories(:).isdir]);
subdirectories = subdirectories(3:end);
deletefiles = {};
counter = 0;

for i = 1:length(subdirectories)
    files = dir(fullfile('sorted_by_events', subdirectories(i).name, '*.mat'));
    for j = 1:length(files)
        filename = fullfile('sorted_by_events', subdirectories(i).name, files(j).name);
        load(filename);
        dt = data.dt;
        num_samples = numel(data.data);
        duration = num_samples * dt;
        if abs(duration) < 10
            counter = counter+1
            files_to_delete{end+1,1} = ['D:\Seismic Data\' filename];
            delete(filename);
        end
    end
end

