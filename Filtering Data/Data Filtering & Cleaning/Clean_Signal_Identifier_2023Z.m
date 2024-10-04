%% Introduction
%The following code, used in Ryker Tracy 2023 Thesis, looks through an
%files fed into the function, finds the SNR, std of SNR, mean of nonsignal
%vs signal, and creates a list of file paths that are very clean.

%Files must be loaded as a structured array that contains files.folder and 
%files.name

%% Program
function clean_signal = Clean_Signal_Identifier_2023Z(files)


% Sorts the files into a clean signal array of file names

clean_signal = {};

for i = 1:length(files)
    % Load the current file
    load([files(i).folder '\' files(i).name]);
    
    % Filter the signal
    filtered_signal = firfilter2023_RT(data.data, data.dt, 1, 'low');
    filtered_signal = firfilter2023_RT(filtered_signal, data.dt, .05, 'high');
    
    % Calculate the SNR
    signal_interval = (8/data.dt):(13/data.dt);
    noise_interval = (2/data.dt):(8/data.dt);
    snr(1,i) = max(abs(filtered_signal(signal_interval))) / max(abs(filtered_signal(noise_interval)));
    snr(2,i) = std(filtered_signal(signal_interval)) / std(filtered_signal(noise_interval));
    snr(3,i) = mean(abs(filtered_signal(signal_interval))) / mean(abs(filtered_signal(noise_interval)));
    
    % Check if the file belongs in the clean signal array
    if snr(1,i) <= 5 && snr(1,i) >= 1.4*snr(2,i) && snr(3,i) > 2.2
        clean_signal{1,end+1} = [files(i).folder '\' files(i).name];
    elseif snr(1,i) > 5 && snr(1,i) >= 1.2*snr(2,i) && snr(3,i) > 2.2
        clean_signal{1,end+1} = [files(i).folder '\' files(i).name];
    elseif snr(1,i) > 15
        clean_signal{1,end+1} = [files(i).folder '\' files(i).name];
    end
end

end
