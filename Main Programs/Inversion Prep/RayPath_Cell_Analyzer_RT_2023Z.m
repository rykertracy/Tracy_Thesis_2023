%% Introduction
% RayPach_Cell_Analyzer_RT_2023Z
% The following program moves across cells within the study area and
% analyzes them for ray paths.
% 'c_' variables stand for variables associated with files that have
% already had a p arrival picked. The c stands for clean because I first
% picked from clean files.

%% Initialization of Variables
fillsquare_filename = {};
fillsquare_pick = [];
fillsquare_az = [];
fillsquare_dist = [];
fillsquare_evlat = [];
fillsquare_evlon = [];
fillsquare_stlat = [];
fillsquare_stlon = [];
useless_cells = [];

%% Check Cell with Current Pick
for lon_middle=-108.75:0.5:-93.25 %Change back to -108.75 %The middle longitude for any cell
    for lat_middle=25.25:0.5:36.75 %Change back to 25.25 %The middle latitude for any cell
        [perpendicular_distance1,dm1,dse1,daz1] = find_dist(lon_middle,lat_middle,c_evlon,c_evlat,c_stlon,c_stlat); %Calculates the perpendicular distance of any line to the center of the cell.
        [perpendicular_distance2,dm2,dse2,daz2] = find_dist(lon_middle,lat_middle,c_stlon,c_stlat,c_evlon,c_evlat); %Calculates the perpendicular distance of any line to the center of the cell.
        dazf=daz1; dsef=dse1; dmf=dm1; perpendicular_distancef=perpendicular_distance1;
        ndaz=find(daz2<daz1);

        dazf(ndaz)=daz2(ndaz);
        dsef(ndaz)=dse2(ndaz);
        dmf(ndaz)=dm2(ndaz);
        perpendicular_distancef(ndaz) = perpendicular_distance2(ndaz);

        %Finding the distance from the middle of the cell to the ray path
        %endpoints
        % [dist1,azimuth1] = distance(lat_middle, lon_middle, c_evlat, c_evlon);
        % [dist2,azimuth2] = distance(lat_middle, lon_middle, c_stlat, c_stlon);
        % delta_azimuth = abs(azimuth2 - azimuth1);

        %Isolate the radii that are within the cell
        in_radius = find(abs(perpendicular_distancef) < 25 & dsef > dmf);

        % figure(1)
        % clf
        % for k = 1:length(in_radius)
        %     rectangle('Position',[lon_middle-.25 lat_middle-.25 0.5 0.5])
        %     title(['Target cell: ', num2str(lat_middle), ' ', num2str(lon_middle)])
        %     hold on
        %     plot([c_stlon(in_radius(k)), c_evlon(in_radius(k))],[c_stlat(in_radius(k)) c_evlat(in_radius(k))])
        % end
        % hi=1
        % pause

        figure(3)
        clf
        %Plot histograms for arc distances and azimuths ALREADY picked.
        subplot(2,1,1)
        histogram(c_dist(in_radius), [0:30:900])
        title('Distance to Events That Cross the Square ')

        subplot(2,1,2)
        histogram(c_az(in_radius), [0:5:360])
        title('Azimuth of Events That Cross the Square')

%% Ask User to look for more ray paths and find ray paths
        %Ask the user what to do with any given cell
        do_what = input('enter a 1 to look for more events a 0 if this looks good');
        if do_what == 1

            %Calculate perpendicular distances of ray paths from total pool
            %of events
            [perpendicular_distance1,dm1,dse1,daz1] = find_dist(lon_middle,lat_middle,EVloc(:,2),EVloc(:,1),STloc(:,2),STloc(:,1)); %Calculates the perpendicular distance of any line to the center of the cell.
            [perpendicular_distance2,dm2,dse2,daz2] = find_dist(lon_middle,lat_middle,STloc(:,2),STloc(:,1),EVloc(:,2),EVloc(:,1)); %Calculates the perpendicular distance of any line to the center of the cell.
            dazf=daz1; dsef=dse1; dmf=dm1; perpendicular_distancef=perpendicular_distance1;
            ndaz=find(daz2<daz1);
            dazf(ndaz)=daz2(ndaz);
            dsef(ndaz)=dse2(ndaz);
            dmf(ndaz)=dm2(ndaz);
            perpendicular_distancef(ndaz) = perpendicular_distance2(ndaz);

            %Find ray paths within the cell radius less than 500 km
            %arcdistance
            in_radius2 = find(perpendicular_distancef < 25 & dsef > dmf);

            %If statement to skip cells with no paths that go through it.
            if isempty(in_radius2)
                disp(['No ray paths for cell with mid-point: ', num2str(lat_middle), ' lat, ', num2str(lon_middle), ' lon. '])
                useless_cells(end+1,1:2) = [lat_middle, lon_middle];
                continue;
            end

            figure(3)
            clf
            %Plot the new potential histograms of distances and azimuths
            subplot(2,1,1)
            histogram(arcd(in_radius2),[0:30:900])
            title('NEW distance to events that cross the square ')

            subplot(2,1,2)
            histogram(azimuth(in_radius2),[0:5:360])
            title('NEW azimuth of events that cross the square')

            [value,ind]=sort(snrMAX(in_radius2),'descend');
            sorted_index = in_radius2(ind);

            %Loop through indices of the master list with the highest
            %signal to noise mean value that allegedly travel through the
            %cell.
            for i = 1:length(sorted_index)
                load(filename_cell{sorted_index(i)},'data')
                D=detrend(data.data);
                FS = firfilter2023_RT(D,data.dt,[.05, 2],'bandpass');
                cut1 = 8/data.dt;
                cut2 = 13/data.dt;
                figure(1)
                clf

                subplot(2,1,1)
                rectangle('Position',[lon_middle-.25 lat_middle-.25 0.5 0.5])
                axis([-110 -92 24 38])
                title(['Target cell: ', num2str(lat_middle), ' ', num2str(lon_middle)])
                hold on
                plot([data.station_lon, data.event_lon],[data.station_lat data.event_lat])

                subplot(2,1,2)
                plot(D*(1/max(abs(D))))
                hold on
                plot(FS*(1/max(abs(FS))))
                subtitle(['Current Event: ', num2str(data.event_lat), ' ', num2str(data.event_lon) ' Current Station: ' num2str(data.station_lat), ' ', num2str(data.station_lon)])
                axis([cut1 cut2 -1 1])

                %Pick the P arrival for the other files. Clicking to the
                %left of the graph does not save the pick. Clicking to the
                %right of the graph indicates satisfaction with catelogue
                %and continues to next cell.
                [pick,~] = ginput(1)
                if pick < cut1
                    clf
                    continue;

                    %Store variables
                elseif pick < cut2 && pick > cut1
                    fillsquare_filename{end+1}= filename_cell{sorted_index(i)};
                    fillsquare_pick(end+1) = pick
                    fillsquare_az(end+1) = azimuth(sorted_index(i));
                    temp_az(i)=azimuth(sorted_index(i));
                    fillsquare_dist(end+1) = arcd(sorted_index(i));
                    temp_dist(i) = arcd(sorted_index(i));
                    fillsquare_evlat(end+1) = EVloc(sorted_index(i),1);
                    fillsquare_evlon(end+1) = EVloc(sorted_index(i),2);
                    fillsquare_stlat(end+1) = STloc(sorted_index(i),1);
                    fillsquare_stlon(end+1) = STloc(sorted_index(i),2);
                    figure(3)
                    clf
                    subplot(2,1,1)
                    histogram([c_dist(in_radius) temp_dist],[0:30:900])
                    title('Fill Distance to Events That Cross the Square ')
                    subplot(2,1,2)
                    histogram([c_az(in_radius) temp_az],[0:5:360])
                    title('Fill Azimuth of Events that Cross the Square')
                elseif pick > cut2
                    break;
                end

            end
        end
    end
    start_lon=lon_middle;
end
