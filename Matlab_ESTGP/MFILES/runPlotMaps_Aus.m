%% runPlotMaps creates the maps in Ashe et al., 2018

    minlat=-36;
    maxlat=-13;
    minlong=112;
    maxlong=128;
    [minlat maxlat minlong maxlong]
    lat_incr = 1;
    long_incr = 1;
    Flat=minlat:lat_incr:maxlat;
    Flong=minlong:long_incr:maxlong;
    
    %runMapHeight_SD;

    runMapLocalRemoved;

    %runMapHighstandProb;