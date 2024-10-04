%make_data_list
% Initialize cell array to store all clean signals
clean_signal_all = {};

% Get a list of all subdirectories in 'sorted_by_events'
subdirs = dir('sorted_by_events');
subdirs = subdirs([subdirs(:).isdir]); % Remove non-directory entries

% Loop through all subdirectories
n=0;
for subdir_idx = 1:length(subdirs)
    if strcmp(subdirs(subdir_idx).name, '.') || strcmp(subdirs(subdir_idx).name, '..')
        continue % Skip current iteration if current directory is '.' or '..'
    end
    subdir_path = fullfile(subdirs(subdir_idx).folder, subdirs(subdir_idx).name);
    files = dir(fullfile(subdir_path, '*.mat')); % Get list of all .mat files in subdirectory
    
    % Loop through all files in current subdirectory
    for file_idx = 1:length(files)
        file_path = fullfile(files(file_idx).folder, files(file_idx).name);
        load(file_path);
        if data.dt == 0.013
            delete(file_path)
            continue;
        end
        filtered_signal = firfilter2023_RT(data.data,data.dt,1,'low');
        filtered_signal = firfilter2023_RT(filtered_signal,data.dt,.05,'high');
        signal_interval = (8/data.dt):(13/data.dt);
        noise_interval = (2/data.dt):(8/data.dt);
        n=n+1;
        snrMAX(n) = max(abs(filtered_signal(signal_interval)))/max(abs(filtered_signal(noise_interval)));
        snrSTD(n) = std(filtered_signal(signal_interval))/std(filtered_signal(noise_interval));
        snrMEAN(n) = mean(abs(filtered_signal(signal_interval)))/mean(abs(filtered_signal(noise_interval)));
        filename(n).name=file_path;
        arcd(n)=data.arc_distance;
        garc(n)=data.gcarc;
        STloc(n,:)=[data.station_lat, data.station_lon, data.station_elevation];
        EVloc(n,:)=[data.event_lat, data.event_lon, data.event_depth];
        magnitude(n)=data.magnitude;
        azimuth(n)=data.azimuth;
        backaz(n)=data.backazimuth;
        dt(n)=data.dt;

        
        
        %Check if current signal is clean and add to clean_signal_all if it is
        % if snr(1,n) <= 5 && snr(1,file_idx) >= 1.4*snr(2,n) && snr(3,n) > 2.2
        %     clean_signal_all{end+1} = file_path;
        % elseif snr(1,n) > 5 && snr(1,n) >= 1.2*snr(2,n) && snr(3,n) >2.2
        %     clean_signal_all{end+1} = file_path;
        % elseif snr(1,n) > 15
        %     clean_signal_all{end+1} = file_path;
        % end

    end
end
save Master_data.mat
