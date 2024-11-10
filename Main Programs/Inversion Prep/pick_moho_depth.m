%% Introduction
% This program allows a user to pick the moho depth, mantle velocity, and
% crustal velocity (Vp) for Pn Inverstion from the "ALL_MODS.mat" variables. The
% variable of primary interest here is the stVp_all, that contains the
% velocity profile beneath all station 'LATS' and 'LONS'. 
%
% stVp_all should be stored as single columns representing a velocity
% profile beneath a station. The amount of columns in stVp_all should be
% the amount of stations are in your study area, or at least stations that
% have an available Vp profile.
%
% The first click from the ginput command is the moho depth. The second
% click is the mantle velocity, and the third is the crustal velocity.
%
%% Program
% For this program to run correctly, you HAVE TO select elements in the
% following order:
% 1) choose the Moho depth
% 2) choose a representative mantle velocity
% 3) choose a representative low-crustal velocity
%
load("All_MODS_RT2023.mat")
for i = 1:width(stVp_all)
    figure(1)
    hold on
    plot(LONS(i),LATS(i),'k*')
    axis([-109 -93 25 37])

    figure(2)
    plot(stVp_all(:,i))
    
    [depth,velocity] = ginput(3);
    MD_selected(i) = depth(1) / 10;
    moho_vp(i) = velocity(1);
    mantle_vp(i) = velocity(2);
    crustal_vp(i) = velocity(3);
    
    clf(figure(2))
end
% save('All_MODS_RT2023.mat','MD_selected','crustal_vp','mantle_vp','moho_vp',"-append")