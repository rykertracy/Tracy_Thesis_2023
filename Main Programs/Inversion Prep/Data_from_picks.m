%The following program collects data from the wavelet filepaths and stores
%it in "Master_current_picks_DATE.mat"

path = updated_paths;
for i=1:length(path)
    load(path{i}, 'data')
    c_az(i) = data.azimuth;
    c_dist(i) = data.arc_distance;
    c_evlat(i) = data.event_lat;
    c_evlon(i) = data.event_lon;
    c_stlat(i) = data.station_lat;
    c_stlon(i) = data.station_lon;
    c_p_arrival(i) = updated_picks(i);
    c_dt(i) = data.dt;
    c_st_elevation(i) = data.station_elevation;
    c_event_depth(i) = data.event_depth;
    c_magnitude(i) = data.magnitude;
        
    % These if statements are used it you want to filter by azimuth.
        % if az(i) >= 180
        %     plot([stlon(i), evlon(i)], [stlat(i) evlat(i)], 'r')
        % elseif az(i) < 180
        %     plot([stlon(i) evlon(i)], [stlat(i) evlat(i)], 'k')
        % end
end
save('Master_current_picks_05182023.mat')


