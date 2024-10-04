# Tracy Thesis 2023
This repo contains the MATLAB codes used in the Texas Tech University thesis: Geophysical Investigation of the Moho Beneath New Mexico, Oklahoma, and Texas.

## TO DO / Notes
Current Task: Load all relevant programs into GitHub. This is done in the **documentation** branch.

Next: Create a workflow from beginning to end of program running.

Then: Update documentation on all programs. Commit each change to a program with a "docs" type.

Then: Begin reworking programs. This includes making them more user friendly, combining files where possible, and amplifying effeciency.

## Original Workflow
Data downloaded from PyWeed was originally in SAC format. The first step was to translate from SAC to MAT through the "SACtoMATB_2023Z.m" program. "rdsac.m" and "SACTOOLB_2023Z.m" are depended on by the program, and I did not write them--the first came from Mathworks and the second from a previous Texas Tech student. 

After MAT files created, "eventsorter_RT_2023Z.m" and "stationsorter_RT_2023Z" will copy and group data into seismic events and recording stations. "Delete_Files_Based_On_wavelet_Length.m" will delete files based on sampling length.
## Code Layout

## Dependencies

## Credit