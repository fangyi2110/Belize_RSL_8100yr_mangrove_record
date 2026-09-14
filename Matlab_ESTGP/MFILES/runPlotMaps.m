%% runPlotMaps creates the maps in Ashe et al., 2018

    minlat=;
    maxlat=27;
    minlong=-84;
    maxlong=-79.5;
    [minlat maxlat minlong maxlong]
    lat_incr = 0.05;
    long_incr = 0.05;
    Flat=minlat:lat_incr:maxlat;
    Flong=minlong:long_incr:maxlong;
    
    runMapHeight_SD;

    runMapRegional;
