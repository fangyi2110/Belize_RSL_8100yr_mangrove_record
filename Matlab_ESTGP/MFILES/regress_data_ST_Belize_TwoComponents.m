function [f2s, sd2s, V2s, testlocs, logp_out] = regress_data_ST_TwinCays( ...
    model_index, datasets, testsitedef, modelspec, thetTGG, trainsub, ...
    testt, refyear, collinear)

% Regresses Holocene RSL data using the selected GP model.
% Inputs:
%   model_index: Index into modelspec and thetTGG for parameter choice
%   datasets: Data structure array with RSL observations
%   testsitedef: Prediction site definitions
%   modelspec: Struct array with kernel/covariance model specifications
%   thetTGG: Cell array of optimized hyperparameters
%   trainsub: Indices of data points used for training
%   testt: Time vector for prediction (e.g., -8000:100:2000)
%   refyear: Reference year (e.g., 2010)
%   collinear: Optional collinearity input (can be empty)

disp(['Running regression with model index: ' num2str(model_index)]);

% Use the passed model index to set jj
jj = model_index;

% Create noise masks for different model variants (Full, Global, Regional Local)
% noiseMasks = ones(3, length(thetTGG{jj}));
% noiseMasks(1,[6]) = 0;     % Full model without white noise
% noiseMasks(2,[3 6]) = 0;   % Global only
% noiseMasks(3,[1 6]) = 0;   % Local only

noiseMasks = ones(4, length(thetTGG{jj}));
noiseMasks(1,[6]) = 0;     % Full model without white noise
noiseMasks(2,[3 6]) = 0;   % Global only
noiseMasks(3,[1 6]) = 0;   % Local only
noiseMasks(4,[]) = 0;   % Full model with white noise included

% Prepare output containers
f2s = cell(1, size(noiseMasks, 1));
sd2s = cell(1, size(noiseMasks, 1));
V2s = cell(1, size(noiseMasks, 1));
logp_out = zeros(1, size(noiseMasks, 1));

% Run regression for each mask variant
for iii = 1:size(noiseMasks, 1)
    [f2s{iii}, sd2s{iii}, V2s{iii}, testlocs, logp_out(iii), ...
     passderivs, invcv] = RegressHoloceneDataSets( ...
        datasets{1}, ...
        testsitedef, ...
        modelspec(jj), ...
        thetTGG{jj}, ...
        trainsub, ...
        noiseMasks(iii,:), ...
        testt, ...
        refyear, ...
        collinear);
end

end
