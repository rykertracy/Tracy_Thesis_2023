function y=butterfilter(data, dt, fc, t)
data = detrend(data)
[data]=taper_2015(data,dt,t,t);

% Define filter specifications
% fc=Cutoff frequency (Hz)
order = 4; % Filter order

% Compute filter coefficients
[b, a] = butter(order, fc*2*dt, 'low');

y = filter(b, a, data);
end