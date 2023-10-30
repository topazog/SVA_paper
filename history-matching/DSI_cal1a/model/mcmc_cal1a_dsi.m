%% obtain the posterior of standard deviates from DSI

%first we need to defined the number of standard deviates depending on the energy level from SVD

%load data realizations
d_full = readmatrix("cal1a_ies.0.obs.csv");
Ne = size(d_full,2);
Nd = size(d_full,1);

%calculate the mean of d_full for each data point
d_prior = repmat(mean(d_full,2),1,Ne);

%calculate the centered matrix X
X = (d_full-d_prior)/((Ne-1)^0.5);

%perform svd to obtain the number of singular values and phi matrix given
%certain level of energy
le = 0.995;
[U,S,V]=svd(X);
sv_squared = diag(S).^2;
total_energy = sum(sv_squared);
cum_energy = cumsum(sv_squared);
max_sv = find(cum_energy >= le * total_energy, 1);

%% Problem settings defined by user
DREAMPar.d = max_sv;                         % Dimension of the problem
DREAMPar.lik = 2;                      % log likelihood
DREAMPar.T = 20000;
DREAMPar.mt = 5;
%DREAMPar.nCR = max(floor(DREAMPar.d/5),1);
%DREAMPar.nCR = 10;
%DREAMPar.thinning = 5;                  % Only store each 5th sample
DREAMPar.N = 5;
%DREAMPar.psnooker = 0.1;

%% Provide information parameter space and initial sampling
Par_info.initial = 'prior';             % Initial sample from prior distribution
priors = cell(1, max_sv);

for i=1:max_sv
    priors{i}='normpdf(0.0,1.0)'; % Marginal prior of z
end

Par_info.prior = priors;    
Par_info.boundhandling = 'reflect';     % Enforce boundaries of parameter space

%% Define feasible parameter space (minimum and maximum values)
% parname		[ log10(range_m) af_m	aa_m]
%Par_info.zmin = -3 *ones(1,900);
%Par_info.zmax = 3 *ones(1,900);

Par_info.min =	-5*ones(1,max_sv);  % For boundary handling
Par_info.max =	5*ones(1,max_sv);  % For boundary handling

%% Define name of function (.m file) for posterior exploration
Func_name = 'MODEL2';

%% Optional settings
options.parallel = 'yes';              % Run chains in parallel
options.IO = 'yes';                    % Input-output writing of model files (for parallel setting only!)
options.save = 'yes';                  % Save memory of DREAM during trial
options.print = 'yes';

%% Run DREAM package
tic
[chain,output,FX,Z] = MTDREAM_ZS(Func_name,DREAMPar,Par_info,[],options);
elapsed_time_min = toc/60;
disp(elapsed_time_min);