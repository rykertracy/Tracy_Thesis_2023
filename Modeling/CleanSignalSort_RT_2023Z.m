

%% Sort Into Clean Signal Array of File Names
clean_signal = {};
files = files(3:end);
for i = 1:length(files)
    load([files(i).folder '\' files(i).name]);
    filtered_signal = firfilter2023_RT(data.data,data.dt,1,'low');
    filtered_signal = firfilter2023_RT(filtered_signal,data.dt,.05,'high');
    signal_interval = (8/data.dt):(13/data.dt);
    noise_interval = (2/data.dt):(8/data.dt);
    snr(1,i) = max(abs(filtered_signal(signal_interval)))/max(abs(filtered_signal(noise_interval)));
    snr(2,i) = std(filtered_signal(signal_interval))/std(filtered_signal(noise_interval));
    snr(3,i) = mean(abs(filtered_signal(signal_interval)))/mean(abs(filtered_signal(noise_interval)));
    if snr(1,i) <= 5 && snr(1,i) >= 1.4*snr(2,i) && snr(3,i) > 2.2
        clean_signal{1,end+1} = [files(i).folder '\' files(i).name];
    elseif snr(1,i) > 5 && snr(1,i) >= 1.2*snr(2,i) && snr(3,i) >2.2
        clean_signal{1,end+1} = [files(i).folder '\' files(i).name];
    elseif snr(1,i) > 15
        clean_signal{1,end+1} = [files(i).folder '\' files(i).name];
    end
end

%% Pick P arrivals for clean signals
for ii = 1:length(clean_signal)
    load(clean_signal{1,ii})
    filtered_signal = firfilter2023_RT(data.data,data.dt,1,'low');
    filtered_signal = firfilter2023_RT(filtered_signal,data.dt,.05,'high');
    signal_interval = (8/data.dt):(13/data.dt);
    noise_interval = (2/data.dt):(8/data.dt);
    snr3 = mean(abs(filtered_signal(signal_interval)))/mean(abs(filtered_signal(noise_interval)))
    PPick = PPicker_RT_2023Z(0.5,filtered_signal(signal_interval),'k');
    PPick = PPick+signal_interval(1)-1;
    plot(detrend(data.data)*(1/max(abs(detrend(data.data)))),'b')
    hold on
    plot(filtered_signal,'r','LineWidth',2)
    xline(PPick,'LineWidth',3)
    question = input('Are you okay with this pick? Press enter if okay, press 1 to change')
    if isempty(question)
        hold off
    else
        [newPPick,y] = ginput(1)
        xline(newPPick, 'r', 'LineWidth', 3);
        hold off
        PPick = newPPick
    end
    clean_signal{2,ii} = PPick;
    data.P_Arrival = PPick; 
    save([clean_signal{1,ii}]);
end