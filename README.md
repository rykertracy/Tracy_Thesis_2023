# Tracy Thesis 2023
### **Disclaimer**

This repository contains the MATLAB codes used in the Texas Tech University [thesis](https://github.com/rykertracy/Tracy_Thesis_2023/blob/master/TRACY-THESIS-2023.pdf): Geophysical Investigation of the Moho Beneath New Mexico, Oklahoma, and Texas. The intention of this repository is to store and provide the programs necessary to execute PnP tomography inversion with the guidance of the thesis collaborators. This repository is not necessarily intended to provide an end-to-end pipeline for PnP tomography. Creating this repository did not occur to me until nearly two years after the programs were produced, and some programs have been lost--most of which were files relating to cleaning the massive amount of data. With that admission in mind, my passion for other projects relating to machine learning and AI have dominated my free time such that I am unwilling to create new versions of the lost data-cleaning files at this time.

That disclaimer aside, the files in this repository will complete a model to be displayed and anlyzed in MATLAB if the data has been properly formatted.

## Directory Structure:
1) Archive: Contains outdated programs that may have use for offshoot applications.
2) Data Importing from Pyweed: Programs to import and transform PyWeed from .sac to .mat file format.
3) Filtering Data: Programs for assessing, cleaning, filtering, and QA/QC.
4) Main Programs: Programs necessary to produce the thesis model.
5) NonCritical Code: Programs that may be useful to some, but not others. Could also be called "helpers".

## Suggested Workflow
1) Import data from PyWeed (IRIS).
2) Filter, QA/QC, clean data.
3) Review and run Main Programs.