%% runPlotMaps creates the maps in Ashe et al., 2018

    minlat=15;
    maxlat=22;
    minlong=-90;
    maxlong=-85;
    [minlat maxlat minlong maxlong]
    lat_incr = .1;
    long_incr = .1;
    Flat=minlat:lat_incr:maxlat;
    Flong=minlong:long_incr:maxlong;

    % Create noise masks for different model variants (Full, Global, Regional Local)
    noiseMasks = ones(4, length(thetTGG{jj}));
    noiseMasks(1,[9]) = 0;     % Full model without white noise
    noiseMasks(2,[3 6 9]) = 0;   % Global only
    noiseMasks(3,[1 6 9]) = 0;     % Regional non-linear only
    noiseMasks(4,[1 3 9]) = 0;   % Local only
   
    clear wdataset;
    wdataset=datasets{1};

    runMapHeight_SD;

    %runMapLocalRemoved;

    %runMapHighstandProb;