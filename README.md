# Tracy Thesis 2023
This repo contains the MATLAB codes used in the Texas Tech University thesis: Geophysical Investigation of the Moho Beneath New Mexico, Oklahoma, and Texas.

## Programs 
The [Main Programs](https://github.com/rykertracy/Tracy_Thesis_2023/tree/master/Main%20Programs) folder contains all files that need to be run after importing data from PyWeed. Running the **qc_main.m** program will execute the programs in both [Data Importing from Pyweed](https://github.com/rykertracy/Tracy_Thesis_2023/tree/master/Data%20Importing%20from%20Pyweed) and [Filtering Data](https://github.com/rykertracy/Tracy_Thesis_2023/tree/master/Filtering%20Data).

> The [Inversion Prep](https://github.com/rykertracy/Tracy_Thesis_2023/tree/master/Main%20Programs/Inversion%20Prep) subfolder contains all of the programs necessary to prepare the cleaned ray path data (see the Main Programs folder for some programs to help clean the data and choose P arrivals) for inversion.
>
> The [PnP Tomography](https://github.com/rykertracy/Tracy_Thesis_2023/tree/master/Main%20Programs/PnP%20Tomography) subfolder contains all programs necessary to perform inversion. There are 4 different inversion programs with will invert based on slightly different parameters (largest difference being resolution).

**qc_main(master_path).m** - Is a combination of functions that will allow for a more streamline QC process. Elements of the combined files have been removed based on lack of use, like event_information tables that were unused. Do access these, see the archive folder.