%We can take this short bursts at a time. As long as I don't save the clean
%signal back to the directory, we're chilling

%Lets first create an array of PPicks, then later add them to the
%file array.
p_arrival = clean_picks_above80km.p_arrival;
path = clean_picks_above80km.path;
for i = 17846:18845
    i;
    load(path{i}, 'data')
    k = data.data*(1/max(data.data));
    k = k-mean(k);
    cut1 = 8/data.dt;
    cut2 = 13/data.dt;
    plot(k)
    axis([cut1 cut2 -1.1 1.1])
    title(i)
    hold on
    filtered_signal = firfilter2023_RT(data.data,data.dt,2,'low');
    filtered_signal = firfilter2023_RT(filtered_signal,data.dt,.05,'high');
    plot(filtered_signal, 'r')
    PP = PPicker2(data.data,data.dt);
    try
        xline(PP,'k')
    catch
        disp('No PP value');

    end
    %also Plot the filtered data
    [pick,~] = ginput(1);
    if pick < cut1
        p_arrival(i) = PP;
    elseif pick > cut2
        p_arrival(i) = 0; %0 means it's not pickable.
    else
        p_arrival(i) = pick;
    end
    clf
end

index = find(p_arrival==0);
p_arrival(index) = [];
path(index) = [];
clean_picks_above80km.path = path;
clean_picks_above80km.p_arrival = p_arrival;
save('clean_picks_above80km9_4.mat','clean_picks_above80km')