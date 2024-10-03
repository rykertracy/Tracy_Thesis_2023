%% Instructions and Details

%For Ryker Tracy thesis 2023, files containing LHZ, DHZ, UHZ, and VHZ were
%chosen to be deleted due to inadequate sampling. The following program
%will delete them from the current directory.

%% Program
% Specify the path to the "sorted_by_station" directory
directory_path = 'D:/Seismic Data/sorted_by_events/';

% Get a list of all subdirectories in the "sorted_by_station" directory
subdirectories = dir(directory_path);
subdirectories = subdirectories([subdirectories(:).isdir]);
subdirectories = subdirectories(3:end)

% Loop through each subdirectory
for i = 1:length(subdirectories)
    % Get the name of the current subdirectory
    subdirectory_name = subdirectories(i).name;
    
    % Skip any subdirectories that start with a dot (i.e., are hidden)
    if subdirectory_name(1) == '.'
        continue;
    end
    
    % Construct the full path to the current subdirectory
    subdirectory_path = fullfile(directory_path, subdirectory_name);
    
    % Change the current working directory to the current subdirectory
    cd(subdirectory_path);
    
    % Execute the code to delete the specified files
    deletefiles1 = dir('*.LHZ_*'); %This variable picks all files including 'LHZ', which we have chosen to throw out
    deletefiles2 = dir('*.DHZ_*'); %This variable picks all files including 'DHZ', which we have chosen to throw out
    deletefiles3 = dir('*.UHZ_*'); %This variable picks all files including 'UHZ', which we have chosen to throw out
    deletefiles4 = dir('*.VHZ_*'); %This variable picks all files including 'VHZ', which we have chosen to throw out
    
    if ~isempty(deletefiles1) || ~isempty(deletefiles2) || ~isempty(deletefiles3) || ~isempty(deletefiles4)
         if ~isempty(deletefiles1)
             for j = 1:length(deletefiles1)
                 delete(deletefiles1(j).name);
                 fprintf('Deleted file: %s\n', deletefiles1(j).name);
             end
         end
         if ~isempty(deletefiles2)
             for j = 1:length(deletefiles2)
                 delete(deletefiles2(j).name);
                 fprintf('Deleted file: %s\n', deletefiles2(j).name);
             end
         end
         if ~isempty(deletefiles3)
             for j = 1:length(deletefiles3)
                 delete(deletefiles3(j).name);
                 fprintf('Deleted file: %s\n', deletefiles3(j).name);
             end
         end
         if ~isempty(deletefiles4)
             for j = 1:length(deletefiles4)
                 delete(deletefiles4(j).name);
                 fprintf('Deleted file: %s\n', deletefiles4(j).name);
             end
         end
    end
    
    % Change the current working directory back to the parent directory
    cd('..');
end
