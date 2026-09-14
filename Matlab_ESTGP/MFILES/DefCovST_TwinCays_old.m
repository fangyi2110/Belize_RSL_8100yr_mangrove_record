% DefCovST defines the covariance functions and bounds on the
% hyperparameters for the analysis of the spatio-temporal dataset in Ashe et al., 2018

refyear=1950; % reference year for GIA calculations (see RegressHoloceneDataSets.m)

%% define covariance functions 
kMat1 = @(dx,thetas) thetas(1).^2 .* (1).*exp(-dx/thetas(2)); % Squared Exponential (Ornstein–Uhlenbeck) kernel
kMat3 = @(dx,thetas) thetas(1).^2 .* (1 + sqrt(3)*dx/thetas(2)).*exp(-sqrt(3)*dx/thetas(2)); % Matern 3/2 kernel
kMat5 = @(dx,thetas) thetas(1).^2 .* (1 + (sqrt(5)*dx/thetas(2)).*(1 + sqrt(5)*dx/thetas(2)/3)).*exp(-sqrt(5)*dx/thetas(2)); % Matern 5/2 kernel
kDELTA = @(dx,thetas) thetas(1).^2 .* (dx==0); % Kronecker Delta kernel (for noise terms)
kDP = @(years1,years2,thetas) thetas(1).^2 * bsxfun(@times,(years1-refyear)',(years2-refyear)); % Deterministic Polynomial kernel (Linear time covariance)
kDELTAG = @(ad,thetas)thetas(1).^2.*(abs(ad)<1e-4).*(ad<360); % Geographic delta kernel (near zero angular distances; to detect localised signals)

kMat3d = @(years1,years2,dx,thetas)thetas(1).^2.*(-3/(thetas(2).^2)).*dx.*exp(-sqrt(3)*dx/thetas(2)).*(-1+2*bsxfun(@ge,years1',years2)); % first derivative of Matern 3/2
kMat5d = @(years1,years2,dx,thetas)thetas(1).^2.*(-5*dx.*exp(-sqrt(5)*dx/thetas(2))).*(thetas(2)+sqrt(5)*dx)/(3*thetas(2).^3); % first derivative of Matern 5/2
kDPd = @(years1,years2,thetas) thetas(1).^2 * repmat((years1-refyear)',length(years2),1); % first derivative of Deterministic Polynomial
kMat3dd = @(dx,thetas)thetas(1).^2.*(3/(thetas(2).^2)).*(1-sqrt(3)*dx/thetas(2)).*exp(-sqrt(3)*dx/thetas(2)); % second derivative of Matern 3/2
kMat5dd = @(dx,thetas)thetas(1).^2.*-(5*exp(-(sqrt(5)*dx)/thetas(2)).*(thetas(2).^2+sqrt(5)*thetas(2).*dx-5*dx.^2))/(3*thetas(2).^4); % second derivative of Matern 5/2
kDPdd = @(years1,years2,thetas) thetas(1).^2 * ones(length(years2),length(years1)); % second derivative of Deterministic Polynomial

clear modelspec;

%% defines higher level global, regional-linear, regional non-linear, local, and white noise covariance functions
 cvfunc.G = @(dt1t2,thetas) kMat3(dt1t2,thetas(1:2)); % global kernel
 cvfunc.L = @(dt1t2,ad,thetas) (kMat3(dt1t2,thetas(1:2))).*(kMat3(ad,[1 thetas(3)])).*(ad<360); % local non-linear
 cvfunc.W = @(dt1t2,ad,thetas) kDELTAG(ad,1).*kDELTA(dt1t2,thetas(1)); % white noise

%% defines first derivatives
 dcvfunc.G = @(t1,t2,dt1t2,thetas) kMat3d(t1,t2,dt1t2,thetas(1:2)); % global
 dcvfunc.L = @(t1,t2,dt1t2,ad,thetas) (kMat3d(t1,t2,dt1t2,thetas(1:2))).*(kMat3(ad,[1 thetas(3)])).*(ad<360); % local non-linear
 dcvfunc.W = 0; % white noise

%% defines second derivatives
 ddcvfunc.G = @(dt1t2,thetas) kMat3dd(dt1t2,thetas(1:2)); % global
 ddcvfunc.L = @(dt1t2,ad,thetas) (kMat3dd(dt1t2,thetas(1:2))).*(kMat3(ad,[1 thetas(3)])).*(ad<360); % local non-linear
 ddcvfunc.W = 0; % white noise
 
%% define parameter ranges for kernels (starting point, lower bound of search, upper bound of search)
% We iterate in a sysetmatic manner to find the best-fitting parameters

% Base theta configuration:
% [global amplitude, global t lengthscale, local amplitude, local t lengthscale, local s lengthscale, whitenoise]
base_thet0 = [500 1 100 100 0.01 100]; % starting point
base_lb    = [1 1 0.01 1 0.001 0.1]; % lower bound
base_ub    = [120000 80000 40000 40000 2 1e6]; % upper bound

% Define parameter indices and scaling parameters to loop
param_index = 1:6;
param_names = {'global_amp', 'global_tscale', 'local_amp', 'local_tscale', 'local_sscale', 'noise'};
scales = [0.01, 0.1, 1, 10];

% set up array to store theta_sets
theta_sets = {};
counter = 1;

for s = scales
    for p = param_index

        % scale thet0 only
        thet0 = base_thet0 * s;
        theta_sets{counter,1} = [thet0; base_lb; base_ub];
        theta_sets{counter,2} = sprintf('vary-thet0-%s-x%.2f', param_names{p}, s);
        counter = counter+1;

        % scale lb only
        lb = base_lb * s;
        theta_sets{counter,1} = [base_thet0; lb; base_ub];
        theta_sets{counter,2} = sprintf('vary-lb-%s-x%.2f', param_names{p}, s);
        counter = counter+1;

        % scale lb only
        ub = base_ub * s;
        theta_sets{counter,1} = [base_thet0; base_lb; ub];
        theta_sets{counter,2} = sprintf('vary-ub-%s-x%.2f', param_names{p}, s);
        counter = counter+1;

    end
end

%% build model specs based on the above parameter sets

clear modelspec

for i = 1:size(theta_sets,1) % loop through all parameter sets (theta_sets)
    theta_values = theta_sets{i,1}; % extract the paramater values
    label = theta_sets{i,2}; % extract the parameter scaling descriptions

    modelspec(i).thet0 = theta_values(1,:); % starting point
    modelspec(i).lb = theta_values(2,:); % lower bound
    modelspec(i).ub = theta_values(3,:); % upper bound

    % covariance functions to use in training model
    modelspec(i).cvfunc = @(t1,t2,dt1t2,thetas,ad,fp1fp2)  cvfunc.G(dt1t2,thetas(1:2)) + cvfunc.L(dt1t2,ad,thetas(3:5)) + cvfunc.W(dt1t2,ad,thetas(6));
    % specify training model based on above covariance functions
    modelspec(i).traincv = @(t1,t2,dt1t2,thetas,errcv,ad,fp1fp2) modelspec(1).cvfunc(t1,t2,dt1t2,thetas,ad,fp1fp2) + errcv;
    % first derivatives: for model prediction from trained model  
    modelspec(i).dcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) dcvfunc.G(t1,t2,dt1t2,thetas(1:2)) + dcvfunc.L(t1,t2,dt1t2,ad,thetas(3:5));    
    % second derivatives: for model prediction from trained model
    modelspec(i).ddcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) ddcvfunc.G(dt1t2,thetas(1:2)) + ddcvfunc.L(dt1t2,ad,thetas(3:5)); 

    modelspec(i).subfixed=[];
    modelspec(i).sublength=[5]; % spatial lengthscale
    modelspec(i).subamp = [1 3 6]; % global, local and white noise amplitudes
    modelspec(i).subamplinear = [];
    modelspec(i).subampnonlinear = [];
    modelspec(i).subampoffset = [];
    modelspec(i).subampnoise = [6]; % white noise
    modelspec(i).subampHF = [];
    modelspec(i).label=label;

end

save('modelspecs.mat', 'modelspec');


