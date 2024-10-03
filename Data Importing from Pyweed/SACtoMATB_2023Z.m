function SACtoMATB_2023Z(outP)
%This program will convert your SAC files into single component mat files and
%this is our hope that it works with pyweed data

A=dir('*HZ*.SAC'); %searches current directory for all BHZ files, sets up structural array which holds info for each file
O=length(A); %number of BHZ files found
mkdir MAT
for i=1:O; %runs loop for all BHZ files
    fnameZ.name=A(i).name; %searches the structural array 'A' for the file name
    SACTOOLB_2023Z(fnameZ,'./'); %runs SACTOOLB_2023Z program
    i=i;
end

