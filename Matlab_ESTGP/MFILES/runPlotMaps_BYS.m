%% runPlotMaps creates the maps in Ashe et al., 2018

    minlat=30;
    maxlat=45;
    minlong=110;
    maxlong=130;
    [minlat maxlat minlong maxlong]
    lat_incr = .25;
    long_incr = .25;
    Flat=minlat:lat_incr:maxlat;
    Flong=minlong:long_incr:maxlong;
    
    runMapHeight_SD;

    %runMapLocalRemoved;

    %runMapHighstandProb;