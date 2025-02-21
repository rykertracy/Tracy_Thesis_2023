# Tracy Thesis 2023
### **Disclaimer**

<<<<<<< HEAD
## Programs 
The [Main Programs](https://github.com/rykertracy/Tracy_Thesis_2023/tree/master/Main%20Programs) folder contains all files that need to be run after importing data from PyWeed. Running the **qc_main.m** program will execute the programs in both [Data Importing from Pyweed](https://github.com/rykertracy/Tracy_Thesis_2023/tree/master/Data%20Importing%20from%20Pyweed) and [Filtering Data](https://github.com/rykertracy/Tracy_Thesis_2023/tree/master/Filtering%20Data).

> The [Inversion Prep](https://github.com/rykertracy/Tracy_Thesis_2023/tree/master/Main%20Programs/Inversion%20Prep) subfolder contains all of the programs necessary to prepare the cleaned ray path data (see the Main Programs folder for some programs to help clean the data and choose P arrivals) for inversion.
>
> The [PnP Tomography](https://github.com/rykertracy/Tracy_Thesis_2023/tree/master/Main%20Programs/PnP%20Tomography) subfolder contains all programs necessary to perform inversion. There are 4 different inversion programs with will invert based on slightly different parameters (largest difference being resolution).

**qc_main(master_path).m** - Is a combination of functions that will allow for a more streamline QC process. Elements of the combined files have been removed based on lack of use, like event_information tables that were unused. Do access these, see the archive folder.
||||||| e64db37
## TO DO / Notes
Current Task: Load all relevant programs into GitHub. This is done in the **documentation** branch.

Next: Create a workflow from beginning to end of program running.

Then: Update documentation on all programs. Commit each change to a program with a "docs" type.

Then: Begin refactoring programs. This includes making them more user friendly, combining files where possible, and amplifying effeciency.

## Programs 
qc_main(master_path).m - Is a combination of functions that will allow for a more streamline QC process. Elements of the combined files have been removed based on lack of use, like event_information tables that were unused. Do access these, see the archive folder.

## Original Workflow
Data downloaded from PyWeed was originally in SAC format. The first step was to translate from SAC to MAT through the "SACtoMATB_2023Z.m" program. "rdsac.m" and "SACTOOLB_2023Z.m" are depended on by the program, and I did not write them--the first came from Mathworks and the second from a previous Texas Tech student. "plot_data_map_RT_2023Z.m" can spatially display the events and stations.

After MAT files created, "eventsorter_RT_2023Z.m" and "stationsorter_RT_2023Z" will copy and group data into seismic events and recording stations. "Delete_Files_Based_On_wavelet_Length.m" will delete files based on sampling length.

The user should use their own methods for determining if sampling (dt) is the same across all files. I wrote my own simple scripts, that are not included here, for such things.

## Code Layout

## Dependencies

## Credit
=======
This repository contains the MATLAB codes used in the Texas Tech University [thesis](): Geophysical Investigation of the Moho Beneath New Mexico, Oklahoma, and Texas. The intention of this repository is to store and provide the programs necessary to execute PnP tomography inversion with the guidance of the thesis collaborators. This repository is not necessarily intended to provide an end-to-end pipeline for PnP tomography. Creating this repository did not occur to me until nearly two years after the programs were produced, and some programs have been lost--most of which were files relating to cleaning the massive amount of data. With that admission in mind, my passion for other projects relating to machine learning and AI have dominated my free time such that I am unwilling to create new versions of the lost data-cleaning files at this time.

That disclaimer aside, the files in this repository will complete a model to be displayed and anlyzed in MATLAB if the data has been properly formatted.

## Directory Structure:
1) Archive: Contains outdated programs that may have use for offshoot applications.
2) Data Importing from Pyweed: Programs to import and transform PyWeed from .sac to .mat file format.
3) Filtering Data: Programs for assessing, cleaning, filtering, and QA/QC.
4) Modeling: Programs necessary to produce the thesis model.
5) NonCritical Code: Programs that may be useful to some, but not others. Could also be called "helpers".

## Suggested Workflow
>>>>>>> documentation
