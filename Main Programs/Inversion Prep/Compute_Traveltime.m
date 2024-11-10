% This program calculates the travel time for seismic ray paths based on geophysical data.

% Define the array 'p_arrival' using the pre-calculated arrival times for P-waves from clean signals (c_p_arrival).
p_arrival = c_p_arrival;

% Loop through each ray path to compute time picks and retrieve event and file times.
for i = 1:length(path)
    % Load seismic data from the current path file, which contains event details.
    load(path{i}, 'data')
    
    % Store the time pick (arrival time) for each path. 
    % 'c_dt' is a sample rate.
    time_pick(i) = p_arrival(i) * c_dt(i);
    
    % Store the event time of the seismic event.
    event_time(i, :) = data.evtime;
    
    % Store the start time of the seismic record (when the recording began).
    file_time(i, :) = data.begintime;
end

% Adjust milliseconds in event times for improved precision in calculations.
% Add fractional seconds (milliseconds) from the 6th column into the seconds column (5th).
event_time(:, 5) = event_time(:, 5) + (event_time(:, 6) ./ 1000);
% Remove the millisecond column after adjustment.
event_time(:, 6) = [];

% Calculate the initial travel time by subtracting the event time from the file time.
traveltime = file_time(:, 1:end) - event_time(:, 1:end);

% Loop through each computed travel time and adjust for cases where crossing into the next minute/hour is needed.
for k = 1:length(traveltime)
    % Adjust for cases where there’s a carry-over in seconds or minutes.
    if traveltime(k, 4) == 1
        % If there's a carryover in minutes (1 minute difference), adjust seconds.
        traveltime(k, 4) = 0;
        traveltime(k, 5) = (60 - event_time(k, 5)) + file_time(k, 5);
    elseif traveltime(k, 4) == 2
        % If there's a two-minute difference, adjust seconds and minutes accordingly.
        traveltime(k, 4) = 0;
        traveltime(k, 5) = (60 - event_time(k, 5)) + 60 + file_time(k, 5);
    elseif traveltime(k, 3) == 1
        % Adjust if there's an hour carry-over (1 hour difference).
        traveltime(k, 3) = 0;
        traveltime(k, 4) = 0;
        traveltime(k, 5) = (60 - event_time(k, 5) + file_time(k, 5));
    end
end

% Add the time pick offsets to the adjusted travel times.
for n = 1:length(traveltime)
    travel_time(n) = traveltime(n, 5) + time_pick(n);
end

% Save the computed travel times to a file, appending to the existing file.
save('Master_current_picks_05182023.mat', 'travel_time', "-append")

% Uncomment the following line to plot arc distance against travel time.
% plot(arcdist, travel_time, '*')
