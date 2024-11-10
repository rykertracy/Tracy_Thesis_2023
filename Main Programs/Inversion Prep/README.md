# Inversion Prep
This folder contains programs necessary to prepare the data for inversion and visualization produced by the PnP Tomography programs in a separate folder.

### Programs
Data_from_picks.m: Stores data from wavelets stored in filepaths into a "master" database (.mat format).

mohopierce_lat_lot_RT_2023Z.m: Calculates and stores the latitude and longitude for each earthquake ray path at the moho depth. It is based on the "ALL_MODS.mat" file that contains a velocity profile beneath each seismic station.

Compute_Traveltime.m: Calculates the travel time from wavelets with a P arrival chosen and appends it to the "master" database.

pick_moho_depth.m: Allows a user to choose the depth of the moho based on vertical velocity data beneath a station stored in "ALL_MODS.mat".

Coverage_map_2023.m: Displays a map of the ray paths to determine, visually, if adequate ray paths are present.

RayPath_Cell_Analyzer_RT_2023Z.m: Allows a user to visually determine if adequate ray paths are present by splitting the study area into "cells" and review the rays that cross that cell for azimuth and length coverage.

fox_RT2023Z.m: convert lat lon values to local x-y coordinates based on data in the "stuff_needed.mat" file, which is described within.

### Sample Data and Files
stuff_needed_Positive_Moho_tt.zip: The latitudes and longitudes of stations, events, and moho piercing points where the moho travel time is only positive.

ALL_MODS_RT2023.mat: Vertical velocity profile variables for the study area.