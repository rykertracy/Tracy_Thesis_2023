function  FindHZ
%This code was used to determine which '*HZ' we want to use.
fs=dir('*HZ_*')
for n=1:length(fs)
    a=find(fs(n).name=='_',1,'first');
    xhz(n)=fs(n).name(a-3);
end
unique(xhz)
end