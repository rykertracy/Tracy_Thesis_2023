function [y]=firfilter2023_RT(data, dt, fc, type)
data = detrend(data);
[data]=taper_2015(data,dt,2,2);
% Define filter specifications
% fc = Cutoff frequency (Hz)
filter_length = 51; % Filter length (odd number)

% Compute filter coefficients
b = fir1(filter_length-1, fc*2*dt, type);

y = filter(b, 1, data);
y = flipud(filter(b,1,flipud(y)));

y = y*(1/max(y)); %Normalize the data
end