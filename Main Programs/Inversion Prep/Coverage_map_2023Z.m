%The following program produces a map of ray paths over a study area

path = Final_P_arrival.path;
p_arrival = Final_P_arrival.p_arrival;

paths = path(1:length(p_arrival));

figure(2)
hold on
for i=1:length(new_paths)
    load(new_paths{i}, 'data')
    az(i) = data.azimuth;
    dist(i) = data.arc_distance;
    evlat(i) = data.event_lat;
    evlon(i) = data.event_lon;
    stlat(i) = data.station_lat;
    stlon(i) = data.station_lon;


    plot([c_stlon(i), c_evlon(i)], [c_stlat(i) c_evlat(i)],'k')
    for j=1:length(latBounds)-1
        for k=1:length(lonBounds)-1
            if stlat(i)>=latBounds(j) && stlat(i)<latBounds(j+1) && ...
               stlon(i)>=lonBounds(k) && stlon(i)<lonBounds(k+1) && ...
               evlat(i)>=latBounds(j) && evlat(i)<latBounds(j+1) && ...
               evlon(i)>=lonBounds(k) && evlon(i)<lonBounds(k+1)
                rayCounts(j,k) = rayCounts(j,k) + 1;
            end
        end
    end
end

for k = 1:length(fillsquare_filename)
    load(fillsquare_filename{k})
    evlat1(k) = data.event_lat;
    evlon1(k) = data.event_lon;
    stlat1(k) = data.station_lat;
    stlon1(k) = data.station_lon;
    plot([stlon1(k), evlon1(k)], [stlat1(k) evlat1(k)],'b')
end

%Plot the ray path coverage for each grid cell
for latIndex=1:size(rayCounts, 1)
    for lonIndex=1:size(rayCounts, 2)
        latCenter = latBounds(1) + (latIndex-0.5)*cellSize;
        lonCenter = lonBounds(1) + (lonIndex-0.5)*cellSize;
        if rayCounts(latIndex, lonIndex) < 3
            plot(lonCenter, latCenter, 'ro', 'MarkerSize', 10)
        end
    end
end