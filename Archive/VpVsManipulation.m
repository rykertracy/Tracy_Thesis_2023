%Setting the Vs data from the USGS model

% load('ALL_MODS_RT2023.mat')
% for i=1:length(vs_index)
%     vs_lat(i) = A(vs_index(i),1);
%     vs_long(i) = A(vs_index(i)+1,1);
%     for k = 1:length(latindex)
%         if vs_index(i) == latindex(k)
%             latindex_next = latindex(k+1);
%             y=k;
%         else
%             continue;
%         end
%     end
%     usgs_mods(y).vs = D(vs_index(i):latindex_next-1);
% end

%Including lat/lon into the USGS model array

% for i=1:length(latindex)
%     usgs_mods(i).latlon = [Q(latindex(i)),-Q(latindex(i)+1)];
% end

%Plotting latlons for stations with vs data

% k=0
% for i=1:length(usgs_mods)
%     if ~isempty(usgs_mods(i).vs)
%         k=k+1;
%         latlon(k,1:2) = [usgs_mods(i).latlon];
%         plot(latlon(k,2),latlon(k,1),'.')
%         hold on
%     end
% end

%Computing Vp/Vs

% k=0
% for i = 1:length(usgs_mods)
%     if ~isempty(usgs_mods(i).vs)
%         tv = [usgs_mods(i).vs; 5];
%         tt = [usgs_mods(i).tops; 120];
%         newd = [0:0.1:120];
%         usgs_mods(i).vs_interp = interp1(tt,tv,newd);
%     end
% end
% for i = 1:length(usgs_vpvs)
%     try
%         tv = [usgs_vpvs(i).vs; 5];
%         tt = [usgs_vpvs(i).tops; 120];
%         newd = [0:0.1:120];
%         usgs_vpvs(i).vs_interp = interp1(tt,tv,newd);
%     catch
%         continue;
%     end
%     tv = [usgs_vpvs(i).vs; 5];
%     tt = [usgs_vpvs(i).tops; 120];
%     newd = [0:0.1:120];
%     usgs_vpvs(i).vs_interp = interp1(tt,tv,newd);
% end

% for i=1:length(usgs_vpvs)
%         usgs_vpvs(i).Vp_Vs = usgs_vpvs(i).vp_interp ./ usgs_vpvs(i).vs_interp;
% end
% k=0;
% for i =1:length(usgs_mods)
%     if ~isempty(usgs_mods(i).Vp_Vs)
%         k=k+1;
%         usgs_vpvs(k) = usgs_mods(i);
%     end
% end
% j = []
% for i = 1:length(lats)
%     try
%     j(end+1) = find(LATS == lats(i) & LONS == lons(i),1);
%     catch
%         continue;
%     end
% end
% for i = 1:length(usgs_vpvs)
%     usgs_vpvs(i).Vp_Vs = usgs_vpvs(i).Vp_Vs';
% end
