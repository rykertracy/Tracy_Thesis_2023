function qc_main(master_path)
%% Documentation
%{ This program is developed for Ryker Tracy Thesis 2023, and aims to run
%the downloaded IRIS data through a Quality Control Process. This program
%is to be run after downloaded data using PyWEED. The program will convert
%the data from .sac to .mat files using the SACtoMATB_2023Z function. Then
%it will sort the files both by seismic events and by recording stations
%using the sorter_RT_2023Z function.

%To run, execute main(master_path) where master_path is the full file path 
% to the location of the data. It is acceptable for the data to be stored
% in subdirectories. When initial run for the thesis, I had my
% subdirectories categorized by event magnitude, which was a parameter
% available in PyWEED upon downloading the data. The program will search
% subdirectories (only one subdirectory layer deep) for .sac files before
% converting.

%       Example: main('C:\Users\ryker\Tracy_Thesis_2023')

% Importantly, in order for the sac to mat program to run, the rdsac.m file
% must be added to the path.

%}
cd(master_path) % Enter master_path directory
subdirs = dir(master_path); % List contents of master_path directory
subdirs = subdirs([subdirs(:).isdir]); % Remove everything non-directory from the variable.
subdirs = subdirs(~ismember({subdirs(:).name},{'.','..'})); % Remove . and ..
SACtoMATB_2023Z; % Run the .sac to .mat program

%Cycle through subdirectories
for i=1:length(subdirs)
    currentDir = fullfile(master_path,subdirs(i).name);
    cd(currentDir);
    sacFiles = dir('*.sac');
    
    %Apply the .sac to .mat conversion if applicable and copy the MAT
    %folder contents to the one in the master_path.
    if ~isempty(sacFiles)
        SACtoMATB_2023Z;
        matPath = fullfile(pwd, 'MAT');
        matFiles = dir(fullfile(matPath, '*.*'));
        matFiles = matFiles(~ismember({matFiles.name}, {'.', '..'})); % Remove '.' and '..'
        for k = 1:length(matFiles)
            sourceFile = fullfile(matPath, matFiles(k).name);
            destinationFile = fullfile(master_path, 'MAT');
            
            copyfile(sourceFile, destinationFile);
        end
    
    %Print that no .sac files are found so that the user can determine if
    %that's expected or not.
    else
        fprintf('No .sac files found in directory %s\n', currentDir);
    end
    
end
cd(master_path) % Return to master_path
sorter_RT_2023Z(master_path) %Perform sorting.
end

function sorter_RT_2023Z(master_path) %FUNCTION REFACTOR IN PROGRESS
%% Documentation
%{ This program was developed to sort all of the PyWEED-downloaded files
%into a folder for seismic events and a folder for recording events (for
%the purposes of ray-tracing). The time of the event and the station name
%are both listed in the title of the file, so this function performs string
%operations to extract that information. The function will create new
%directories for events and stations within the a 'sorted_by_events' and
%'sorted_by_station' subdirectory respectively.

%The function will run with the execution of the main(master_path)
%function, but for stand alone running:

%       Example: sorter_RT_2023Z('C:\Users\ryker\Tracy_Thesis_2023')

%}
%% Information Extraction
% Defining directory path strings
cd(master_path); mkdir  sorted_by_events; mkdir sorted_by_station; %Change path to master directory and create a new folder for sorting by events and storing by station.
events_path = [master_path, '\sorted_by_events']; % Define the 'sorted_by_events' path string.
stations_path = [master_path, '\sorted_by_station']; % Define the 'sorted_by_station' path string.

% Extract information from the files
files=[]; % Initialize the files array
files_t=dir(['MAT']); % List the contents of the MAT subfolder (when run with the main() function, it will run the MAT in the master_path).

for n=3:length(files_t) % Loop through files for relevant information.
    file=[master_path '\MAT\' files_t(n).name]; % Isolate specific file.
    load(file); % Load file
    files(n).name = files_t(n).name; % Store file name
    files(n).eventnames = data.evdate_string; % Store event name for event sorting
    files(n).eventlat = data.event_lat; % Store event lat for event sorting
    files(n).eventlon = data.event_lon; %Store event lon for event sorting
end

%% Sorting
%{ It was found that the same event was recorded several times across all
%subdirectories. This is likely because of rounding errors between events.
%To prevent this, the following code was produced to round the event time
%and gather all of the information into a unique event.

%The following block will, additionally, isolate information needed to sort
%by station.

%}

filename = {files(3:end).name};
newfilename = extractBefore(filename,'HZ'); % Extract station name (which is stored before the channel 'HZ' in the string

% Loop through the file and shift the indices of the string.
for i = 1:length(newfilename);
    newfilename{i} = newfilename{i}(1:end-2);
end
unique_station = unique(newfilename); %Store unique station names for grouping.

% Define event names
events = {files(3:end).eventnames};
rounded_events = cellfun(@(x) [x(1:end-2)], events, 'UniformOutput', false); % Round events (strings)
lat = {files(3:end).eventlat}; % Store event lat
lon = {files(3:end).eventlon}; % Store event lon
unique_events = unique(rounded_events); % Isolate unique strings

%% Optimize file names, create folders with unique events as names, and copy files into them.
%Replace the colons in the event times with underscores for file creation.
rounded_events = cellfun(@(x) strrep(x, ':', '_'), rounded_events, 'UniformOutput', false);
unique_events = cellfun(@(x) strrep(x, ':', '_'), unique_events, 'UniformOutput', false);

cd(events_path); %Enter sorted by events directory
%Create a bunch of directories with unique event times as names
for ii=1:numel(unique_events)
    mkdir(unique_events{ii}); 
end

%Take the filename that corresponds to the event time and place into new
%directory
for ii=1:length(files(3:end))
    current_rounded_event = rounded_events{ii};
    current_filename = [master_path '\MAT\' filename{ii}];.
    copyfile(current_filename, [events_path '\' current_rounded_event]);
end 

%% Make Directories based on station name
%Create a bunch of directories with unique stations as names
cd(stations_path);
for ii=1:numel(unique_station)
    mkdir(unique_station{ii});
end

%Take the filename that corresponds to the station and place it in the
%directory
for ii=1:length(files(3:end))
    current_station = newfilename{ii};
    current_filename = [master_path '\MAT\' filename{ii}];
    copyfile(current_filename, [stations_path '\' current_station]);
end
cd(master_path) % Return to master_path
end