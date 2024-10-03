%% Instructions and Details
%This program is used for Ryker Tracy Thesis 2023. It sorts files
%downloaded from PyWEED into folders titled the time of the event. Each
%event directory will contain all files that claim it as the associated
%event. Admittedly, it's not the most effecient, where you will need to
%change all lines that contain desired directories every time you want to
%run it.

%I strongly recommend reading all comments in this program before running
%it to understand relevant syntax.

%% Create 'sorted_by_events' directory, then sort through files and gather relevant information
mkdir 'D:\Seismic Data\sorted_by_events' %Change Path to reflect directory
files=[];
files_t=dir(['Pn_DATA_3p5_3p6_2014_2023_East_Apply Shit\MAT']); %Change with for every folder you want to sort

for n=3:length(files_t) %Goes through the files and gathers relevant information
    file=['D:\Seismic Data\Pn_DATA_3p5_3p6_2014_2023_East_Apply Shit\MAT\' files_t(n).name]; %Change directory accordingly.
    load(file);
    files(n).name = files_t(n).name;
    files(n).eventnames = data.evdate_string;
    files(n).eventlat = data.event_lat;
    files(n).eventlon = data.event_lon;
end

%% Define variables and sort them into arrays.
%We found that the same event was recorded several times with a small fraction of a second difference.
%The following bit will round the event time and gather all of the information into a unique event.  
filename = {files(3:end).name};
events = {files(3:end).eventnames};
rounded_events = cellfun(@(x) [x(1:end-2)], events, 'UniformOutput', false);
lat = {files(3:end).eventlat};
lon = {files(3:end).eventlon};
tableofshit = [filename; rounded_events; lat; lon];
stuff = [rounded_events; lat; lon];
unique_events = unique(rounded_events);

%Create a table of event times, latitude, and longitude to determine if the
%difference in .001 seconds between unique events is the
% same lat and lon
%or if they're truly different events. It appears to be a rounding error.
[~, idx] = ismember(unique_events, rounded_events);
stuff1 = stuff(:,idx);

%% Optimize file names, create folders with unique events as names, and copy files into them.
%Replace the colons in the event times with underscores for file creation.
rounded_events = cellfun(@(x) strrep(x, ':', '_'), rounded_events, 'UniformOutput', false);
unique_events = cellfun(@(x) strrep(x, ':', '_'), unique_events, 'UniformOutput', false);


%Create a bunch of directories with unique event times as names
for ii=1:numel(unique_events);
    mkdir(['D:\Seismic Data\sorted_by_events\' unique_events{ii}]);
end

%Take the filename that corresponds to the event time and place into new
%directory
for ii=1:length(files(3:end));
    current_rounded_event = rounded_events{ii};
    current_filename = ['D:\Seismic Data\Pn_DATA_3p5_3p6_2014_2023_East_Apply Shit\MAT\' filename{ii}]; %Change path to reflect copy.
    copyfile(current_filename, ['D:\Seismic Data\sorted_by_events\' current_rounded_event]); %Change Path here too.
end 
