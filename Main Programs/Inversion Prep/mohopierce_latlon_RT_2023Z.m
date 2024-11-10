function [cStlat,cStlon,cEQlat,cEQlon,Moho_tt] = mohopierce_latlon_RT_2023Z(stationlat,stationlon,eventlat,eventlon,depth,realtime,outelev)
%% Introduction
%this program corrects for initial lat/lon to Moho penetration lat/lon for
%both station and event lat/lons, also produces Travel time correction for
%travel time along moho specifically.

% Program written by Harold Gurrola for Matthew Tave at TTU October 25, 2012
% ORIGINAL PROGRAM NAME: new_latlon_HG

% Program was modified by Ryker Tracy & Harold Gurrola for Ryker Tracy
% Thesis summer 2023. The original program was not commented well
% (typical), so we left some lines in here we didn't end up using or
% thinking was necessary.

% depth = Earthquake event depth
% Stlat = station latitude
% Stlong = station longitude
% EQlat = latitude of Earthquake
% EQlon = longitude of Earthquake
% realtime = traveltime for P wave from event to receiver
% outelev = elevation of station. NOTE: station elevation was not necessary
%   for computations in Tracy's thesis. 

%% Program
% All_MODS_RT2023 contains Velocity Profile data from stations. This data is a concatenation of data processed by Tave, Castille, Harrington, and a USGS catelogue. 
load('ALL_MODS_RT2023_Handpicked.mat')

%Depths every 0.1 km to 120 km. Gurrola/Tracy think its odd that this
%program interpolates to 0.1 km, then back to 1 km.... but whatever.
Z=[0:0.1:120];

for i=1:length(stationlat)
    
    %individual station distance to all lats/lons for Vp modeled stations
    [DST,AZST]=distance(stationlat(i),stationlon(i),LATS,LONS); 

    %Find the minimum distance between select station and modeled stations
    [val,IST]=min(abs(DST)); 

    %Interpolate velocity profile beneath nearest station to 1 km depth ranges.
    VST=interp1(Z,stVp_all(:,IST),[1:1:(MD_selected_rounded(IST)+2)]);   % use moho depth fromt  1st click ************************************

    %i8 is the mohodepth index for the nearest velocity profile, plus 2 kilometers.
    i8=MD_selected_rounded(IST)+2; 
    
    % if isempty(i8)==0, VST(i8)=7.5; end

    %Use the deepest velocity value (plus a little) to represent the moho.
    vmohoST = max(VST)+1;

    %Reset the VST to not include the extra 2 kilometers.
    VST = VST(1:end-2);

    %Find the distances between a selected event and lat/lons with velocity
    %profiles
    [Dev,AZev]=distance(eventlat(i),eventlon(i),LATS,LONS);

    %Find the shortest distance between the event and lat/lon velocity
    %profile
    [val,IEV]=min(abs(Dev));

    %Interpolate the event velocity proile to 1 km, beginning at the event
    %depth.
    VEV=interp1(Z,stVp_all(:,IEV),[depth(i):1:(MD(IEV)+2)]);


    i8=MD(IEV)+2;              %find(VEV>7.5);
    % if isempty(i8)==0, VEV(i8)=7.5; end

    %Use the deepest velocity value (plus a little) to represent the moho.
    vmohoEV = max(VEV)+1;

    %Reset the VEV to not include the extra 2 kilometers.
    VEV = VEV(1:end-2);

    %Ray Tracing beneath the station, using the station's interpolated
    %velocity profile.
    theta_ST=asin(crustle_vp1(IST)/mantle_vp1(IST));                            %theta_ST=asin(VST(end)/vmohoST); %make VST(end) the crust velocity from second ginput, and make vmohoST the third ginput.***********
    thetaST=asin(VST(1:end)*sin(theta_ST)/VST(end));
    XST=sum(tan(thetaST(1:end)));
    XST=km2deg(XST);
    dtST=sum((1./VST(1:end))./cos(thetaST(1:end))); %+((outelev(i)./VST(1))./cos(thetaST(1)));
    
    %Ray Tracing beneath the event, using the interpolated velocity profile
    %beneath the event.
    theta_EV=asin(crustle_vp1(IEV)/mantle_vp1(IEV));                                                                        %theta_EV=asin(VEV(end)/vmohoEV);
    thetaEV=asin(VEV(1:end)*sin(theta_EV)/VEV(end));
    XEV=sum(tan(thetaEV(1:end)));
    XEV=km2deg(XEV);
    dtEV=sum((1./VEV(1:end))./cos(thetaEV(1:end)));
    
    %Travel time correction (time taken to reach/exit to station from the moho 
    tt_c=dtEV+dtST; 
    
    %Travel time along moho = total P wave travel time - time taken to
    %reach the moho.
    Moho_tt=realtime-tt_c;
    
    %Compute arc distances between the events and stations
    [arclength_EQ2St,azEQ]=distance(eventlat,eventlon,stationlat,stationlon);
    [arclength_St2EQ,azSt]=distance(stationlat,stationlon,eventlat,eventlon);
    
    %Find the latitude and longitude of where the ray path pierces the
    %moho.
    [cEQlat,cEQlon]=reckon(eventlat,eventlon,XEV,azEQ);
    [cStlat,cStlon]=reckon(stationlat,stationlon,XST,azSt);
    
end

save('stuff_needed_Handpicked_narrow','cStlat','cStlon','cEQlat','cEQlon','Moho_tt')
%After this code runs, you need to change the lat/lon to UTM and you're
%ready to invert! YAY!

