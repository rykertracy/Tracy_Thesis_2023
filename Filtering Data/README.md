# Filtering Data
The programs in this folder filter the data after downloading from IRIS with PyWEED and being transitioned to MATLAB files.

## Programs 
eventsorter_RT_2023Z.m - sorts files into their respective seismic event. Created March 30, 2023.

stationsorter_RT_2023Z.m - groups files by their respective recording station. Created March 30, 2023.

Delete_Files_Based_On_wavelet_Length.m - removes files from the sorted_by_events directory. Can also be applied to sorted_by_stations if the code is edited. Created April 7, 2023.

Station_viewer_RT_2023Z.m - will allow a user to view and select files from stations that do not visually look correct. **WARNING: NEEDS TO SAVE A Bad Stations.mat FILE IN ORDER FOR OTHER PROGRAMS TO RUN**.

Throw_away_files_with_majority_zero_amplitudes.m - **NEEDS DOCUMENTATION**.

butterfilter_RT.m && firfilter2023_RT.m - Filters for analyzing wavelets. Created in February and March of 2023 respectively. **NEED DOCUMENTATION**.

Find_Duplicates_In_Events_Directory.m - Finds files with identical names across events. Sometimes fringe errors in code would cause this.

Assign_Magnitude_to_event.m - Assigns the events magintude to the directory name. **NEEDS DOCS**.

Blocky_Signal_Finder.m - Sometimes wavelets would have blocky visual errors. This codes finds most of those digitizer error problems. **NEEDS DOCS**

Delete_badfiles_from_events.m - Delets files from the events directories. **NEEDS DOCS AND TAILORING**

Clean_Signal_Identifier_2023Z.m - Finds signals with a clean signal to noise ratio and places it into an array called 'clean_signal.mat'.

make_data_list.m - creates a .mat file that serves as a Master for all the data available. Created on April 14, 2023.



## Notes to self/TO DO 
Waveform_Auditor was important, but I don't quite remember to how to operate it. It's not in gitHub yet. let's make sure we understand how it's used properly.

