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
kLinPoly = @(years1, years2, thetas) ... % Second order Polynomial (Quadratic) + linear
    thetas(1).^2 * ...
    bsxfun(@times, ...
        ((years1 - refyear) .* ((years1 - refyear) ./ ((years1 - refyear) + thetas(2))).^2)', ...
        ((years2 - refyear) .* ((years2 - refyear) ./ ((years2 - refyear) + thetas(2))).^2));
kLinPolyDecay = @(years1, years2, thetas) ... % Second order Polynomial (Quadratic) + with decay
    thetas(1).^2 .* ...
    bsxfun(@times, ...
        ((years1 - refyear) ./ (years1 - refyear + thetas(2)).^2)', ...
        ((years2 - refyear) ./ (years2 - refyear + thetas(2)).^2));


% first derivatives
kMat1d = @(dx,thetas) - (thetas(1).^2 ./ thetas(2)) .* exp(-dx ./ thetas(2));
kMat3d = @(years1,years2,dx,thetas)thetas(1).^2.*(-3/(thetas(2).^2)).*dx.*exp(-sqrt(3)*dx/thetas(2)).*(-1+2*bsxfun(@ge,years1',years2)); % first derivative of Matern 3/2
kMat5d = @(years1,years2,dx,thetas)thetas(1).^2.*(-5*dx.*exp(-sqrt(5)*dx/thetas(2))).*(thetas(2)+sqrt(5)*dx)/(3*thetas(2).^3); % first derivative of Matern 5/2
kDPd = @(years1,years2,thetas) thetas(1).^2 * repmat((years1-refyear)',length(years2),1); % first derivative of Deterministic Polynomial
kLinPolyd = @(years1, years2, thetas) ... % First derivative of Second order Polynomial (Quadratic) + linear
    thetas(1).^2 * ...
    bsxfun(@times, ...
        (((years1 - refyear).^3 + 3 * thetas(2) * (years1 - refyear).^2) ...
        ./ ((years1 - refyear) + thetas(2)).^3)', ...
        ((years2 - refyear).^3 ./ ((years2 - refyear) + thetas(2)).^2));
kLinPolyDecayd = @(years1, years2, thetas) ... % First derivative of linear-decaying kernel
    thetas(1).^2 * ...
    bsxfun(@times, ...
        (((1 - (2 * (years1 - refyear) ./ (years1 - refyear + thetas(2)))) ./ (years1 - refyear + thetas(2))).^2)', ...
        ((years2 - refyear) ./ ((years2 - refyear + thetas(2)).^2)) ...
    );

% second derivatives
kMat1dd = @(dx,thetas) (thetas(1).^2 ./ thetas(2).^2) .* exp(-dx ./ thetas(2));
kMat3dd = @(dx,thetas)thetas(1).^2.*(3/(thetas(2).^2)).*(1-sqrt(3)*dx/thetas(2)).*exp(-sqrt(3)*dx/thetas(2)); % second derivative of Matern 3/2
kMat5dd = @(dx,thetas)thetas(1).^2.*-(5*exp(-(sqrt(5)*dx)/thetas(2)).*(thetas(2).^2+sqrt(5)*thetas(2).*dx-5*dx.^2))/(3*thetas(2).^4); % second derivative of Matern 5/2
kDPdd = @(years1,years2,thetas) thetas(1).^2 * ones(length(years2),length(years1)); % second derivative of Deterministic Polynomial
kLinPolydd = @(years1, years2, thetas) ... % second derivative of Second order Polynomial (Quadratic) + linear
    thetas(1).^2 * ...
    bsxfun(@times, ...
        (((years1 - refyear).^3 + 3 * thetas(2) * (years1 - refyear).^2) ...
        ./ ((years1 - refyear) + thetas(2)).^3)', ...
        (((years2 - refyear).^3 + 3 * thetas(2) * (years2 - refyear).^2) ...
        ./ ((years2 - refyear) + thetas(2)).^3) ...
    );
kLinPolyDecaydd = @(years1, years2, thetas) ... % second derivative of Second order Polynomial (Quadratic) + linear + decay
    thetas(1).^2 * ...
    bsxfun(@times, ...
        ((2 * (years1 - refyear) - 4 * thetas(2)) ./ ((years1 - refyear + thetas(2)).^4))', ...
        ((2 * (years2 - refyear) - 4 * thetas(2)) ./ ((years2 - refyear + thetas(2)).^4)) ...
    );

clear modelspec;

%% defines higher level global, regional-linear, regional non-linear, local, and white noise covariance functions
 cvfunc.G = @(dt1t2,thetas) kMat5(dt1t2,thetas(1:2)); % global kernel
 % cvfunc.RN = @(t1,t2,ad,thetas) kLinPolyDecay(t1, t2, thetas(1:2)).* kMat3(ad, [1 thetas(3)]) .* (ad < 360); % regional non-linear
 cvfunc.RN = @(t1,t2,dt1t2,ad,thetas) (kMat5(dt1t2,thetas(1:2))).*(kMat3(ad,[1 thetas(3)])).*(ad<360); % regional non-linear
 cvfunc.L = @(dt1t2,ad,thetas) (kMat3(dt1t2,thetas(1:2))).*(kMat3(ad,[1 thetas(3)])).*(ad<360); % local non-linear
 cvfunc.W = @(dt1t2,ad,thetas) kDELTAG(ad,1).*kDELTA(dt1t2,thetas(1)); % white noise

%% defines first derivatives
 dcvfunc.G = @(t1,t2,dt1t2,thetas) kMat3d(t1,t2,dt1t2,thetas(1:2)); % global
 % dcvfunc.RN = @(t1,t2,ad,thetas)(kLinPolyDecayd(t1, t2,thetas(1:2))).*(kMat3(ad,[1 thetas(3)])).*(ad<360); % kLinPoly RN
 dcvfunc.RN = @(t1,t2,dt1t2,ad,thetas)(kMat5d(t1,t2,dt1t2,thetas(1:2))).*(kMat3(ad,[1 thetas(3)])).*(ad<360); % regional non-linear
 dcvfunc.L = @(t1,t2,dt1t2,ad,thetas) (kMat3d(t1,t2,dt1t2,thetas(1:2))).*(kMat3(ad,[1 thetas(3)])).*(ad<360); % local non-linear
 dcvfunc.W = 0; % white noise

%% defines second derivatives
 ddcvfunc.G = @(dt1t2,thetas) kMat3dd(dt1t2,thetas(1:2)); % global
 % ddcvfunc.RN = @(t1,t2,ad,thetas) (kLinPolyDecaydd(t1, t2, thetas(1:2))).*(kMat3(ad,[1 thetas(3)])).*(ad<360); % kLinPoly RN
 ddcvfunc.RN = @(dt1t2,ad,thetas) (kMat5dd(dt1t2,thetas(1:2))).*(kMat3(ad,[1 thetas(3)])).*(ad<360); % regional non-linear
 ddcvfunc.L = @(dt1t2,ad,thetas) (kMat3dd(dt1t2,thetas(1:2))).*(kMat3(ad,[1 thetas(3)])).*(ad<360); % local non-linear
 ddcvfunc.W = 0; % white noise
 
%% define parameter ranges for kernels (starting point, lower bound of search, upper bound of search)
% We iterate in a sysetmatic manner to find the best-fitting parameters

% Base theta configuration:
base_thet0 = [3000 3000 500 1000 0.25 500 500 0.01 100]; % starting point
base_lb    = [ 1000  1000   10  500 0.05   10  300 0.005  0.1]; % lower bound
base_ub    = [100000 10000  8000 2000 0.60  5000 1000 0.1 10000]; % upper bound

% Define parameter indices and scaling parameters to loop
param_index = 1:9;
param_names = {'global_amp', 'global_tscale', 'regional_nl_amp', 'regional_nl_tscale','regiona_nl_sscale','local_amp', 'local_tscale', 'local_sscale', 'noise'};
scales = [0.1, 1, 10]; % 0.01 0.1 1 10

% set up array to store theta_sets
theta_sets = {};
counter = 1;

% combinations to scale bounds for: binary classifications
% [thet0 lb ub]
% 1 for scale, 0 for don't scale
modes = {
    % [1 0 0], 'thet0';
    % [0 1 0], 'lb';
    % [0 0 1], 'ub';
    [1 1 0], 'thet0+lb';
    [1 0 1], 'thet0+ub';
    [0 1 1], 'lb+ub';
    % [1 1 1], 'thet0+lb+ub';
};

for s = scales %3
    for p = param_index %9
        for m = 1:length(modes) %3-4

        mode = modes{m,1}; % extract the mode to scale
        tag = modes{m,2}; % label for warning message if any

        % always start with the base
        thet0 = base_thet0;
        lb = base_lb;
        ub = base_ub;

        if mode(1) % if mode(1) == 1 (TRUE), scale thet0
            thet0(p) = base_thet0(p) * s;
        end

        if mode(2) % if mode(2) == 1 (TRUE), scale lb
            lb(p) = base_lb(p) * s;
        end 

        if mode(3) % if mode(2) == 1 (TRUE), scale ub
            ub(p) = base_ub(p) * s;
        end

        % ignore run invalid run (i.e., if lb > ub)
        if any(lb > ub)
            warning('Skipping invalid run: scaling %s %s by %0.2f resulted in lb>ub',param_names{p},tag,s);
            continue; % skip this iteration if lower bound exceeds upper bound
        end

        % store current run
        theta_sets{counter,1} = [thet0; lb; ub];
        theta_sets{counter,2} = sprintf('vary-%s-%s-x%.2f', tag,param_names{p}, s);
        counter = counter+1;
        end
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
    % modelspec(i).cvfunc = @(t1,t2,dt1t2,thetas,ad,fp1fp2)  cvfunc.G(dt1t2,thetas(1:2)) + cvfunc.RN(t1,t2,ad,thetas(3:5)) + cvfunc.L(dt1t2,ad,thetas(6:8)) + cvfunc.W(dt1t2,ad,thetas(6));
    modelspec(i).cvfunc = @(t1,t2,dt1t2,thetas,ad,fp1fp2)  cvfunc.G(dt1t2,thetas(1:2)) + cvfunc.RN(t1,t2,dt1t2,ad,thetas(3:5)) + cvfunc.L(dt1t2,ad,thetas(6:8)) + cvfunc.W(dt1t2,ad,thetas(6));
    % specify training model based on above covariance functions
    modelspec(i).traincv = @(t1,t2,dt1t2,thetas,errcv,ad,fp1fp2) modelspec(1).cvfunc(t1,t2,dt1t2,thetas,ad,fp1fp2) + errcv;
    % first derivatives: for model prediction from trained model  
    % modelspec(i).dcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) dcvfunc.G(t1,t2,dt1t2,thetas(1:2)) + dcvfunc.RN(t1,t2,ad,thetas(3:5)) + dcvfunc.L(t1,t2,dt1t2,ad,thetas(6:8));    
    modelspec(i).dcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) dcvfunc.G(t1,t2,dt1t2,thetas(1:2)) + dcvfunc.RN(t1,t2,dt1t2,ad,thetas(3:5)) + dcvfunc.L(t1,t2,dt1t2,ad,thetas(6:8));    
    % second derivatives: for model prediction from trained model
    % modelspec(i).ddcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) ddcvfunc.G(dt1t2,thetas(1:2)) + ddcvfunc.RN(t1,t2,ad,thetas(3:5)) + ddcvfunc.L(dt1t2,ad,thetas(6:8)); 
    modelspec(i).ddcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) ddcvfunc.G(dt1t2,thetas(1:2)) + ddcvfunc.RN(dt1t2,ad,thetas(3:5)) + ddcvfunc.L(dt1t2,ad,thetas(6:8)); 
    modelspec(i).subfixed=[];
    modelspec(i).sublength=[5 8]; % spatial lengthscale
    modelspec(i).subamp = [1 3 6 9]; % global, regional, local and white noise amplitudes
    modelspec(i).subamplinear = [];
    modelspec(i).subampnonlinear = [];
    modelspec(i).subampoffset = [];
    modelspec(i).subampnoise = [9]; % white noise
    modelspec(i).subampHF = [];
    modelspec(i).label=label;

end

save('modelspecs.mat', 'modelspec');


