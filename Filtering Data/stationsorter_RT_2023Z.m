%% Instructions and Details
%The following program, written for Ryker Tracy Thesis 2023, will sort all
%files into their unique stations be accessing the first few characters in
%the file name string. Like the event sorter, it lacks efficiency and the
%directores to access needs to be changed each time the program is ran.


%% Make directores and access current directores.
mkdir 'D:\Seismic Data\sorted_by_station' %Change Path to reflect directory
files=[];
files_t=dir(['Pn_DATA_3p5_3p6_2014_2023_East_Apply Shit\MAT']); %Change with for every folder you want to sort

for n=3:length(files_t) %Goes through the files and gathers relevant information
    file=['D:\Seismic Data\Pn_DATA_3p5_3p6_2014_2023_East_Apply Shit\MAT\' files_t(n).name]; %Change directory accordingly.
    load(file);
    files(n).name = files_t(n).name;
end

%% Store station name
%The following bit will round the event time and gather all of the information into a unique event.  
filename = {files(3:end).name};
newfilename = extractBefore(filename,'HZ');
for i = 1:length(newfilename);
    newfilename{i} = newfilename{i}(1:end-2);
end
unique_station = unique(newfilename);

%Create a table of event times, latitude, and longitude to determine if the
%difference in .001 seconds between unique events is the
% same lat and lon
%or if they're truly different events. It appears to be a rounding error.
[~, idx] = ismember(unique_station, newfilename);
listofstations = newfilename(idx);


%% Make Directories based on station name
%Create a bunch of directories with unique event times as names
for ii=1:numel(unique_station);
    mkdir(['D:\Seismic Data\sorted_by_station\' unique_station{ii}]);
end

%Take the filename that corresponds to the event time and place into new
%directory
for ii=1:length(files(3:end));
    current_station = newfilename{ii};
    current_filename = ['D:\Seismic Data\Pn_DATA_3p5_3p6_2014_2023_East_Apply Shit\MAT\' filename{ii}]; %Change path to reflect copy.
    copyfile(current_filename, ['sorted_by_station\' current_station]); %Change Path here too.
end
