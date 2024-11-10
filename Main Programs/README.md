# Main Programs
This folder contains all of the programs needed for PnP inversion. Run **qc_main.m FIRST**

## Folders
Inversion Prep: Contains files and some sample data largely to prepare the ray paths for inversion programs.

PnP Tomography: Contains both "smooth" and "size" inversion models that perform inversion with slightly different parameters.

## Programs
qc_main.m: Runs quality control programs on the data first.

CleanSignalSort_RT_2023Z.m - Displays clean signals and allows the user to click where the P wave arrives. Relies on external programs: firfilter and PPicker

plot_data_map_RT_2023Z.m - Should be run immediately after SACtoMAT. It gives an idea of the clusters of events and stations spatially. Created March 2, 2023.

Manual_P_Pick.m - Allows a user to click on the location of the P wave arrival on a wavelet and store the pick location.