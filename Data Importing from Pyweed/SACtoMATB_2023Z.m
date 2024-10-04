function SACtoMATB_2023Z(~)
%This program will convert your SAC files into single component mat files and
%this is our hope that it works with pyweed data

A=dir('*HZ*.SAC'); %searches current directory for all BHZ files, sets up structural array which holds info for each file
O=length(A); %number of BHZ files found
mkdir MAT
for i=1:O %runs loop for all BHZ files
    fnameZ.name=A(i).name; %searches the structural array 'A' for the file name
    SACTOOLB_2023Z(fnameZ,'./'); %runs SACTOOLB_2023Z program
    i=i;
end

function SACTOOLB_2023Z(filename,~)
%% Introduction Comments
%{
Program altered for Ryker Tracy thesis 2023 "Imaging the MOHO beneath New
Mexico, Texas, and Oklahoma"

This code is used for reading SAC files produced by IRIS' PyWEED. This
particular program will NOT prep the data for receiver functions. That
portion was removed to better code efficiency, but may be present in
previous versions. The "Z" in the funciton name indicates that this code is
to be used only for vertical component seismograms (in the Tracy thesis is
it used to assist in picking P arrivals for Pn tomography).

PyWEED organizes the files by network and station, whereas the previously
used JWEED (disconinued) organized the files by events. It is no longer
necessary to create directories of data because all of PyWEED's data from
any given range of parameters is stored into a single folder of the user's
choice. If the user's specific needs require the creation of a directory,
please see the comment section of previous versions of SACTOOLB.


Also ommitted from previous versions is the TTU format flag. See previous
versions if the format flag is necessary.

IMPORTANT ALTERATIONS from previous versions:
 1) PyWEED produces different SAC file format than JWEED did. The fields
 produced later in this code are different.

 2) The time of the data is stored in days since Jan 1, 0000. However, the
 header information is Julian day time. 

 3) PyWEED stores the initial time of recording in the title of the file.
 The last time of recording is also stored there.

 4) Previous versions of this program call 'READSAC'. That does not appear
 to function well with PyWEED. A function called 'rdsac' was downloaded
 from MathWorks and stores a structured array of time, data, and header
 information.
%}

%% Initialization
[a,b]=size(filename);
jout=0;
%streamtime = 0;

for i=1:a;
    for j=1:b;
        %determine string file name
        name=[filename(i,j).name];
        fname=[name];
  
%% Read SAC
        sacdata=rdsac(name); %This uses a rdsac (Read Sac) feature downloaded from Mathworks. It was not built by Dr. Gurrola or a student of Dr. Gurrola.

       %filetime=jconvnamehg_2023Z(name); %This is not necessary with PyWEED. The event time used to be stored in the title of JWEED files, and this program is altered from JWEED use.
        %In case I need to use it again, filetime was listed between nztime and streamtime in the title

%% Extract first sample date from file title
        T=extractAfter(name,'HZ_'); 
        T= extractBefore(T,['_' T(1:4)]);
        year=str2num(T(1:4));
        hour=str2num(T(12:13));
        min=str2num(T(15:16));
        sec=str2num(T(18:end));
        T2=extractAfter(sacdata.HEADER.KZDATE,'(');
        T2=extractBefore(T2,')');
        JDAY=str2num(T2);
        if sacdata.HEADER.NZHOUR==0 && hour==24
            EVJDAY=JDAY+1;
        else
            EVJDAY=JDAY;
        end

%% Construct MAT file with header information and data        
        begintime=[year,EVJDAY,hour,min,sec];
        streamtime=(sacdata.t(end,:)*24*60*60)-(sacdata.t(1,:)*24*60*60);
        h2=sacdata.HEADER;
        if isfield(h2,'EVLA')==0
            sacdata.HEADER.EVLA = -1234
            sacdata.HEADER.EVLO = -1234
            sacdata.HEADER.EVDP = -1234
            sacdata.HEADER.MAG = -1234
        end
        
        [arc_distance,azimuth] = distance(sacdata.HEADER.EVLA,sacdata.HEADER.EVLO,sacdata.HEADER.STLA,sacdata.HEADER.STLO,[6378.137,0.081819]);
        [arc_distance,backazimuth] = distance(sacdata.HEADER.STLA,sacdata.HEADER.STLO,sacdata.HEADER.EVLA,sacdata.HEADER.EVLO,[6378.137,0.081819]);
        gcarc = arc_distance*180/pi/6371;
        fields={'name','import','evtime','evdate_string','begintime','streamtime',...
            'station_lat','station_lon','station_elevation','borehole_depth','compass_az',...
            'compass_inc','event_lat','event_lon','event_depth','magnitude',...
            'azimuth','backazimuth','arc_distance','gcarc','SACT0',...
            'dt','data'};
 
 %The following ensures that if the data is Z component, it will set compas
 %azimuth to 0.
 CMPAZt = extractBefore(fname,'_20');
 CMPAZt = CMPAZt(end);
 if CMPAZt == 'Z', sacdata.HEADER.CMPAZ = 0; end
 if CMPAZt == 'N', sacdata.HEADER.CMPAZ = 0; end  %North Component
 if CMPAZt == 'E', sacdata.HEADER.CMPAZ = 90; end %East Component 
 
        fill={fname,datestr(datenum(now)),sacdata.HEADER.NZDTTM,[sacdata.HEADER.KZDATE, sacdata.HEADER.KZTIME],begintime,streamtime,...
            sacdata.HEADER.STLA,sacdata.HEADER.STLO,sacdata.HEADER.STEL,sacdata.HEADER.STDP,sacdata.HEADER.CMPAZ,...
            sacdata.HEADER.CMPINC,sacdata.HEADER.EVLA,sacdata.HEADER.EVLO,sacdata.HEADER.EVDP,sacdata.HEADER.MAG,...
            azimuth,backazimuth,arc_distance,gcarc,10,...
            round(sacdata.HEADER.DELTA,3),sacdata.d};
        
        d=cell2struct(fill,fields,2);
        data(j)=d;     
    end
end
fnameoutb=['MAT/' fname(1:end-3) 'mat'];
save(fnameoutb,'data');