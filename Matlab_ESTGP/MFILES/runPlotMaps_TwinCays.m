%% runPlotMaps creates the maps in Ashe et al., 2018

    minlat=16.80;
    maxlat=16.85;
    minlong=-88.15;
    maxlong=-88.10;
    [minlat maxlat minlong maxlong]
    lat_incr = .05;
    long_incr = .05;
    Flat=minlat:lat_incr:maxlat;
    Flong=minlong:long_incr:maxlong;
    
    runMapHeight_SD;

    %runMapLocalRemoved;

    %runMapHighstandProb;