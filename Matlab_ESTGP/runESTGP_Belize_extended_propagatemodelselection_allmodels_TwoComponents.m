%% Empirical Spatio-Temporal Gaussian Process analysis for Tan et al. 2025 (submitted to Comms Earth and Env)
% Uses Ashe et al. (2019) model.

% change the '~/' to the directory where you have downloaded or cloned the main files.
% Comment out section 1 to skip model regression if regression had
% previously been done.

clear all
clc

% set up file paths 
cd ('/Users/fangyi/Desktop/Belize_RSL_8100yr_mangrove_record/Matlab_ESTGP');
addpath('MFILES');
addpath('IFILES');
IFILES=[pwd '/IFILES'];

% create directory to store outputs
date_field='Belize';
label='_EST_GP_extended_modelensemble_modeofallmodels_TwoComponents';
% label='_EST_GP_publishedonly';
WORKDIR=[date_field label];
if ~exist(WORKDIR,'dir')
    mkdir(WORKDIR);
end
cd(WORKDIR);

% import proxy data 
idHolo = 3e4;
datPX=importdata(fullfile(IFILES,'Belize_6Aug25_extended.csv')); % Input age in yrs BP
% This file segregates sites by tectonic region.
% We incorporate pre-Holocene RSL data from Barbados to constrain early
% Holocene RSL and prevent the model from going to zero at the beginning of
% the record (zero mean function). However, we only use Belize data to
% train the model. Barbados data are just used to extend the prediction.

% defines the period of time of the training dataset
HolID = find(datPX.data(:,11) ~= 999); % find ID for Belize (excluding Barbados)
oldest=-max(datPX.data(HolID,6)+datPX.data(HolID,7)); % define training time period based on Belize only.

%%%%%%
% prepare data
%%%%%%

prep_data_ST;

%%%%%%
% define covariance function
%%%%%%

DefCovST_Belize_TwoComponents;

%%%%%%
%% 1. Optimise covariance and do GP regression (condition prior GP with optimized hyperparameters on the data) and plot RSL and rates
%%%%%%

n_models = size(modelspec,2); % specify number of model spec variants here
thetTGG = cell(1, n_models);
logp = zeros(1, n_models);
trainsubset = cell(1, n_models);
opt_time = zeros(1, n_models);

% --- Prepare cleaned site names ---
Loc = cell(1, length(datasets{1}.sitenames));
for ii = 1:length(datasets{1}.sitenames)
    Loc{ii} = datasets{1}.sitenames{ii};
    % Remove spaces
    sub = strfind(Loc{ii}, ' ');
    Loc{ii} = Loc{ii}(setdiff(1:length(Loc{ii}), sub));
    % Remove slashes
    sub = strfind(Loc{ii}, '/');
    Loc{ii} = Loc{ii}(setdiff(1:length(Loc{ii}), sub));
end

% --- Prepare oldest and youngest age arrays ---
oldest = zeros(1,length(Loc));
youngest = zeros(1,length(Loc));
for ii = 1:length(Loc)
    sub = find(datasets{1}.datid == datasets{1}.siteid(ii));
    if length(sub) > 0
        oldest(ii) = floor(min(datasets{1}.time1(sub)));
        youngest(ii) = ceil(max(datasets{1}.time2(sub)));
    else
        youngest(ii) = 2000;
        oldest(ii) = 2010;
    end
end

% --- Create testsitedef structure ---
clear testsitedef;
testsitedef.sites = [];
testsitedef.names = {};
testsitedef.names2 = {};
testsitedef.firstage = [];

for ii = 1:length(Loc)
    si = find(datasets{1}.datid == datasets{1}.siteid(ii));
    site_lat = datasets{1}.sitecoords(ii,1);
    site_long = datasets{1}.sitecoords(ii,2);
    si = si(1); % take first index
    
    testsitedef.sites(end+1,:) = [datasets{1}.datid(si), site_lat, site_long];
    testsitedef.names2 = {testsitedef.names2{:}, datasets{1}.sitenames{ii}};
    testsitedef.names = {testsitedef.names{:}, Loc{ii}};
    testsitedef.firstage = [testsitedef.firstage min(oldest(ii), -17000)];

end

% specify prediction time frame to regress data (predict to time span
% including Barbados record)
maxage=1950-max(datPX.data(:,6)+datPX.data(:,7));
minage=1950-min(datPX.data(:,6)+datPX.data(:,7));
testt=[maxage:100:minage]; % in years CE

% For each model spec, optimise hyperparameters and run regression

for m = 1:length(modelspec)
    fprintf('Optimizing and regressing with model %d...\n', m);
    rng(20); % set seed

    % Optimize for each model
    [thetPX, trainsubset{m}, logp(m)] = OptimizeHoloceneCovariance( ...
        datasets{1}, modelspec(m), modelspec(m).thet0, [], [], 0, 0, 0); % SPECIFY TRAINING DATASET TIME LIMIT HERE (FOURTH VARIABLE)
    thetTGG{m} = thetPX(1:end-1);

    % Run regression with optimized parameters
    [f2s{m}, sd2s{m}, V2s{m}, testlocs, logp_out(m,:)] = regress_data_ST_Belize_TwoComponents( ...
        m, datasets, testsitedef, modelspec, thetTGG, trainsubset{m}, ...
        testt, refyear, []);
end

histogram(logp,binwidth=0.5); % based on all models, choose "optimal" models as mode

% Visualise logp values
% assign most stable models as optimal models, and the highest logp of
% those as the chosen model
optimalmodelIDs=find(logp >= 3921.5 & logp <= 3922); % histogram mode based on all models
allmodelindex = 1:length(logp);
nonoptimalmodelIDs=allmodelindex(~ismember(allmodelindex,optimalmodelIDs));
chosenmodelID = find(logp== max(logp(optimalmodelIDs)));
figure;
% yyaxis left
scatter(allmodelindex(nonoptimalmodelIDs), logp(nonoptimalmodelIDs), 'x','DisplayName',sprintf('Non-optimal model'));
hold on
scatter(allmodelindex(optimalmodelIDs), logp(optimalmodelIDs), 'o','DisplayName',sprintf('Optimal model'));
hold on
scatter(allmodelindex(chosenmodelID), logp(chosenmodelID),'filled','LineWidth',2,'DisplayName',sprintf('Chosen model'));
ylabel('logp');
% yyaxis right
% scatter(allmodelindex(nonoptimalmodelIDs), abs(logp_out(nonoptimalmodelIDs,1)), 'x','DisplayName',sprintf('Non-optimal model'));
% hold on
% scatter(allmodelindex(optimalmodelIDs), abs(logp_out(optimalmodelIDs,1)), 'o','DisplayName',sprintf('Optimal model'));
% hold on
% scatter(allmodelindex(chosenmodelID), abs(logp_out(chosenmodelID,1)),'filled','LineWidth',2,'DisplayName',sprintf('Chosen model'));
% ylabel('abs ( logp\_out )');
legend('Location','northoutside','Orientation','horizontal');
fig = gcf;
set(fig, 'PaperUnits', 'inches');
set(fig, 'PaperSize', [11 6]);           % width x height
print('model_logps','-dpdf','-bestfit')

% Export optimised hyperparameters from model variants
thetMatrix = cell2mat(thetTGG');  % Transpose to get models as rows
thetMatrix = [(1:size(thetMatrix,1))' thetMatrix]; % add model ID
thetMatrix = [thetMatrix logp' logp_out(:,1)]; % add logp
param_names = {'modelID','global_amp', 'global_tscale','local_amp', 'local_tscale', 'local_sscale', 'noise', 'logp', 'logp_out'};
thetTable = array2table(thetMatrix, 'VariableNames', param_names);
plotstatus = repmat("non-optimal model",size(thetMatrix,1),1);
% specify which models are optimal
    for n=1:length(optimalmodelIDs)
        optimalindex = find(allmodelindex==optimalmodelIDs(n));
        plotstatus(optimalindex)="optimal model";
    end
% specify chosen model 
plotstatus(find(allmodelindex==chosenmodelID))="chosen model";
% add plot status
thetTable.plotstatus = plotstatus;
cd ('/Users/fangyi/Desktop/Belize_RSL_8100yr_mangrove_record/Matlab_ESTGP/Belize_EST_GP_extended_modelensemble_modeofallmodels_TwoComponents/');
writetable(thetTable, 'optimised_hyperparameters.csv');

% Query the bounds for the chosen model
chosenmodel_upr = modelspec(chosenmodelID).ub;
chosenmodel_lwr = modelspec(chosenmodelID).lb;
sprintf('%s %d %d %d %d %.2f %.1f \n %s %d %d %d %d %.2f %d',"Lower bound", chosenmodel_lwr,"Upper bound:", chosenmodel_upr)

%% 2. Export RSL results 

RSLMatrix = [];
pos_RSL_samples = [];
for c=1:3
    for m=1:length(modelspec)

        clear RSL RSLsd yearBP siteID;

        % Below, {m}{c}, 
            % {m} is model number
            % {1} extracts the full posterior output without white noise;
            % {2} extracts global component only;
            % {3} extracts local component only
        RSL = f2s{m}{c};
        RSLsd = sd2s{m}{c};
        yearBP = 1950-testlocs.X(:,3);
        siteID = (testlocs.reg-30000)/1000;

         % add model ranking
        if strcmp(thetTable.plotstatus(m), 'non-optimal model')
            model_ranking = 2;
        elseif strcmp(thetTable.plotstatus(m), 'optimal model')
            model_ranking = 1;
        elseif strcmp(thetTable.plotstatus(m), 'rejected model')
            model_ranking = -1;
        else
            model_ranking = 0;
        end

        if c==2 % for global RSL

            % filter to one (first) site only, as global component is common to all sites
            idx = find(testlocs.reg==testlocs.sites(1,1));
            RSL = RSL(idx);
            RSLsd = RSLsd(idx);
            yearBP = yearBP(idx);
            siteID = repmat(0,length(RSL),1); % set unique identifier for global rates 

            % sample and append global RSL 
            clear samples
            nsamps = 1000; % number of posterior samples to generate
            rng(20); % set seed
            RSLcov = V2s{m}{2}(idx,idx);
            samples = mvnrnd(RSL, RSLcov, nsamps)'; % sample multivariate normal distribution from mean and covariance function
            samples=samples(:); % reshape samples to have one column
            pos_RSL_samples = [pos_RSL_samples; repmat(m,length(samples),1) repmat(yearBP,nsamps,1) repmat(model_ranking,length(samples),1) samples];
        
        else 

        % Exclude barbados from output csv
        idx = find(siteID~=999); 
        RSL = RSL(idx);
        RSLsd = RSLsd(idx);
        yearBP = yearBP(idx);
        siteID=siteID(idx);

        end

        % append RSL 
        RSLMatrix = [RSLMatrix; repmat(c,length(RSL),1) repmat(m,length(RSL),1) RSL RSLsd yearBP siteID repmat(model_ranking,size(RSL))];
    
    end
end

% Format RSL as table and export as .csv
RSL_headers = {'component', 'modelID', 'RSL','RSLsd','yearBP','siteID','model_ranking'};
RSLTable = array2table(RSLMatrix,'VariableNames',RSL_headers);
cd ('/Users/fangyi/Desktop/Belize_RSL_8100yr_mangrove_record/Matlab_ESTGP/Belize_EST_GP_extended_modelensemble_modeofallmodels_TwoComponents');
writetable(RSLTable, 'RSL.csv'); 

% Format posterior RSL samples as table and export as .csv
RSL_headers = {'modelID','yearBP','model_ranking','RSL'};
RSLTable = array2table(pos_RSL_samples,'VariableNames',RSL_headers);
writetable(RSLTable, '/Users/fangyi/Desktop/Belize_RSL_8100yr_mangrove_record/Matlab_ESTGP/Belize_EST_GP_extended_modelensemble_modeofallmodels_TwoComponents/GlobalRSL_samples.csv');

%% 3. Calculate RSL rates 

difftimestep=100; % Temporal spacing (in years) between time steps for computing and plotting rate estimates.

% Find time-paired rows separated by exactly 'difftimestep' units
% For all sites, and all times, taken from <testlocs.X(:,3)>
% : represents all sites (lat, lon) 
% 3 is all subtimes, where subtimes=find(testt>=firstage(i));
% Mdiff returns 1 if time(i) == time(j); -1 if time(i) == time(j) +
% difftimestep; and 0 otherwise
Mdiff = bsxfun(@eq,testlocs.X(:,3),testlocs.X(:,3)')-bsxfun(@eq,testlocs.X(:,3),testlocs.X(:,3)'+ difftimestep);

% Make sure that Mdiff only considers difftimestep pairs from within the
% same region (not across regions)
% testreg lists unique site IDs (from testlocs.sites)
% 31000 to 38000 are our data sites (defined regions)
% 1029000 is Barbados (pre-Holocene data)
Mdiff = Mdiff .* bsxfun(@eq,testlocs.reg,testlocs.reg');

% Remove rows without any matching time-difference pairs (i.e., only want rows with difftimestep)
sub=find(sum(Mdiff,2)==0);
Mdiff=Mdiff(sub,:);

% Compute mean age of time difference pair 
% abs(Mdiff) gives binary whether or not the data is a valid time pair
% (turns -1 and 1 into just 1), versus 0 for not a timewise pair.
% multiplying abs(Mdiff) sums the two times involved in the pair (matrix multiplication)
% sum(abs(Mdiff),2) finds the number of times, in this case 2 for start and end time
% bsxfun(@rdivide ...) then divides sum of times by number of times to get mean age of time-difference pair
difftimes=bsxfun(@rdivide,abs(Mdiff)*testlocs.X(:,3),sum(abs(Mdiff),2)); 
diffreg=bsxfun(@rdivide,abs(Mdiff)*testlocs.reg,sum(abs(Mdiff),2)); % similarly this gives the region's index (average is the same as this is all from same region)

% Calculate finite-difference operator to convert RSL to rates using time relationship between adjacent pairs
Mdiff=bsxfun(@rdivide,Mdiff,Mdiff*testlocs.X(:,3)); % testlocs.X(:,3) is test times for all sites

RateMatrix=[]; % To store rates (df2s, dsd2s)
pos_rate_samples = []; % To store posterior global rate samples
% Calculate rates
for c=1:4
    for m=1:length(modelspec)
    
        clear df2s dV2s dsd2s;
            % Below, {m}{c}, 
            % {1} extracts the full posterior output without white noise;
            % {2} extracts global component only;
            % {3} extracts regional non-linear component only; 
            % {4} extracts local component only
            for n=1:size(f2s{m}{c},2)
                df2s(:,n)=Mdiff*f2s{m}{c}(:,n); % mean function for rates
                dV2s(:,:,n)=Mdiff*V2s{m}{c}(:,:,n)*Mdiff'; % covariance matrix for rates
                dsd2s(:,n)=sqrt(diag(dV2s(:,:,n))); % sd of rates
            end
        
        % convert difftimes to BP instead of CE
        difftimesBP = 1950 - difftimes;
    
        % add model ranking
        if strcmp(thetTable.plotstatus(m), 'non-optimal model')
            model_ranking = 2;
        elseif strcmp(thetTable.plotstatus(m), 'optimal model')
            model_ranking = 1;
        elseif strcmp(thetTable.plotstatus(m), 'rejected model')
            model_ranking = -1;
        else
            model_ranking = 0;
        end
        
        if c==2 % for global rates

            % filter to one (first) site only, as global component is common to all sites
            idx = find(diffreg==testlocs.sites(1,1));
            df2s = df2s(idx);
            dsd2s = dsd2s(idx);
            difftimesBP = difftimesBP(idx);
            dV2s = dV2s(idx,idx);
            siteID = repmat(0,length(df2s),1); % set unique identifier for global rates 

            % sample and append posterior global rates 
            clear samples
            nsamps = 1000; % number of posterior samples to generate
            rng(20); % set seed
            samples = mvnrnd(df2s, dV2s, nsamps)'; % sample multivariate normal distribution from mean and covariance function
            samples=samples(:); % reshape samples to have one column
            pos_rate_samples = [pos_rate_samples; repmat(m,length(samples),1) repmat(difftimesBP,nsamps,1) repmat(model_ranking,length(samples),1) samples];
        
        else

        % get site ID
        siteID=(diffreg-30000)/1000;
        idx = find(siteID~=999); % exclude barbados from output csv
        df2s = df2s(idx);
        dsd2s = dsd2s(idx);
        difftimesBP = difftimesBP(idx);
        dV2s = dV2s(idx,idx);
        siteID=siteID(idx);

        end

        % append rates 
        RateMatrix = [RateMatrix; repmat(c,length(df2s),1) repmat(m,length(df2s),1) df2s dsd2s difftimesBP siteID repmat(model_ranking,size(df2s))];
    end
end

% Format rates as table and export as .csv
rate_headers = {'component', 'modelID', 'df2s','dsd2s','yearBP','siteID','model_ranking'};
RateTable = array2table(RateMatrix,'VariableNames',rate_headers);
cd ('/Users/fangyi/Desktop/Belize_RSL_8100yr_mangrove_record/Matlab_ESTGP/Belize_EST_GP_extended_modelensemble_modeofallmodels_TwoComponents');
writetable(RateTable, 'RSLRates.csv');

% Format posterior rate samples as table and export as .csv
rate_headers = {'modelID','yearBP','model_ranking','rate'};
RateTable = array2table(pos_rate_samples,'VariableNames',rate_headers);
writetable(RateTable, '/Users/fangyi/Desktop/Belize_RSL_8100yr_mangrove_record/Matlab_ESTGP/Belize_EST_GP_extended_modelensemble_modeofallmodels_TwoComponents/GlobalRates_samples.csv');

%% Plot RSL and rate maps for chosen model %%

% set limits for RSL plots
ylims=[ -80e3 20e3;   % Full
        -80e3 20e2;   % Global
        -30e3 10e3];    % Local 
lab=[{'Full'} {'Global'} {'Local'}];
xlim = [-6550 -69]; % x axis limits in years CE

% specify chosen model
% chosenmodelID = 37;

% Define output folder for RSL plots and results
outFolder = '/Users/fangyi/Desktop/Belize_RSL_8100yr_mangrove_record/Matlab_ESTGP/Belize_EST_GP_extended_modelensemble_modeofallmodels_TwoComponents/ChosenModelResults';

% Create folder if it does not exist
if ~exist(outFolder, 'dir')
    mkdir(outFolder);
end

% Change to output folder
cd(outFolder);

for iii=1:3
    if iii==1
        plot_dat=1;
    else
        plot_dat=0;
    end
    makeplot_slrate(datasets{1},f2s{chosenmodelID}{iii},sd2s{chosenmodelID}{iii},V2s{chosenmodelID}{iii},testlocs,lab{iii},2,100,22,xlim,meanSL,plot_dat,[],[],[],[],ylims(iii,:));
end

testX=[];
    testreg=testlocs.reg;
    testsites=testlocs.sites;
    testX(:,1:2)=testlocs.X(:,1:2);
    testX(:,3)=1950-testlocs.X(:,3);
    firstyr = -3000:500:1000;
    lastyr = -2000:500:2000;

jj = chosenmodelID; % Use the passed model index to set jj
runPlotMaps_Belize_TwoComponents;

%% Generate predictions at data time points to calculate residuals for chosen model

% specify prediction time frame to regress data (predict to time span
% including Barbados record)
testt_check = 1950-datPX.data(find(datPX.data(:,9)==0),6); % use SLIPs only
% ages are in years CE

% Run regression with optimized parameters
rng(20); % set seed
[f2s_check, sd2s_check, V2s_check, testlocs_check, logp_out_check(:)] = regress_data_ST_Belize_TwoComponents( ...
chosenmodelID, datasets, testsitedef, modelspec, thetTGG, trainsubset{chosenmodelID}, ...
testt_check, refyear, []);

% Export RSL
RSLMatrix_check = [];
clear RSL_check RSLsd_check yearBP_check siteID_check;

% {1} extracts the full posterior output without white noise;
RSL_check = f2s_check{4};
RSLsd_check = sd2s_check{4};
yearBP_check = 1950-testlocs_check.X(:,3);
siteID_check = (testlocs_check.reg-30000)/1000;

% add model ranking
if strcmp(thetTable.plotstatus(chosenmodelID), 'non-optimal model')
model_ranking = 2;
elseif strcmp(thetTable.plotstatus(chosenmodelID), 'optimal model')
model_ranking = 1;
elseif strcmp(thetTable.plotstatus(chosenmodelID), 'rejected model')
model_ranking = -1;
else
model_ranking = 0;
end

% Exclude barbados from output csv
idx = find(siteID_check~=999); 
RSL_check = RSL_check(idx);
RSLsd_check = RSLsd_check(idx);
yearBP_check = yearBP_check(idx);
siteID_check=siteID_check(idx);

% append RSL 
RSLMatrix_check = [RSLMatrix_check; repmat(1,length(RSL_check),1) RSL_check RSLsd_check yearBP_check siteID_check repmat(model_ranking,size(RSL_check))];

% Format RSL as table and export as .csv
RSL_headers = {'component', 'RSL','RSLsd','yearBP','siteID','model_ranking'};
RSLTable_check = array2table(RSLMatrix_check,'VariableNames',RSL_headers);
cd ('/Users/fangyi/Desktop/Belize_RSL_8100yr_mangrove_record/Matlab_ESTGP/Belize_EST_GP_extended_modelensemble_modeofallmodels_TwoComponents');
writetable(RSLTable_check, 'RSL_check.csv'); 