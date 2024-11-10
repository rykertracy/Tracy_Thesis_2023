
for m=1:length(usgs_mods)
    m=m
    if min(usgs_mods(m).vps)==0 || max(usgs_mods(m).tops)<25
        disp('skipping')
    else
        clear newzs newvps tvps tzs
        hgvp=usgs_mods(m).vps;
        hgtops=usgs_mods(m).tops;
        k=1;
        tvps(k)=hgvp(k);
        tzs(k)=hgtops(k);
        for n=2:length(hgvp)
            k=k+1;
            tvps(k)=hgvp(n-1);
            tzs(k)=hgtops(n)-0.1;
            k=k+1;
            tvps(k)=hgvp(n);
            tzs(k)=hgtops(n);
        end
        k=k+1;
        tvps(k)=8.5;
        tzs(k)=120;
        newzs=[0:0.1:120];
        newvps=interp1(tzs,tvps,newzs);

        plot(tvps,-tzs)
        hold on
        plot(newvps,-newzs,'r.')
        stVp_all=[stVp_all newvps'];
        pause
        % LATS=[LATS lat(m)];
        % LONS=[LONS lon(m)];
    end
end