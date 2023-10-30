%load data realizations
d_full = readmatrix("cal1a_ies.0.obs.csv");
Ne = size(d_full,2);
Nd = size(d_full,1);

%load obs data
data_model = readmatrix("obs_data.csv");
nobs = size(data_model,1);

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
phi = U*S;
tphi = phi(1:end,1:max_sv);
d_PCA_cov = tphi*tphi';
d_PCA_sigma = diag(d_PCA_cov);

%load MCMC results and get 50% of the chains
load('MT-DREAM_ZS.mat')
nchains = 5;
len_chain = 20000;
ParSet = GenParSet(chain)';
ParSubset = zeros(max_sv,nchains*len_chain/2);
chain_size = size(ParSubset,2);
for i=1:nchains
   psss = (i-1)*len_chain/2+1;
   psse = i*len_chain/2;
   pss = (i-1)*len_chain+len_chain/2+1;
   pse = i*len_chain;
   ParSubset(:,psss:psse)=ParSet(1:max_sv,pss:pse);
end

% Compute the empirical CDF for each observation
[f_data,x_data] = arrayfun(@(row) ecdf(d_full(row,:)), 1:Nd, 'UniformOutput', false);

%define the mean vector with size equals to chain size
d_prior = repmat(mean(d_full,2),1,chain_size);

%calculate d_fPCA
d_PCA = d_prior+tphi*ParSubset(1:max_sv,1:chain_size);

%calculate results (d_f_m)
d_f_m = zeros(Nd,chain_size);
d_prior=d_prior(:,1);

f_I_d_PCA = normcdf(d_PCA(:,1),d_prior,d_PCA_sigma.^0.5);
f_I_d_PCA_c = num2cell(f_I_d_PCA');
d_f = cellfun(@(a,b,c) interp1qr(a(1:end),b(1:end),c),f_data,x_data,f_I_d_PCA_c,'UniformOutput', false);
d_f = reshape(cell2mat(d_f),[Nd 1]);

res = d_f(1:nobs,1)-data_model;
C  = (0.01)^(2)*eye(nobs,nobs);
invC = inv(C);
log_F = log ( ( ( 2 * pi )^( - nobs / 2 ) ) * det(C)^( - 1 / 2 ) );
if ( nobs > 150 )
    log_F = 0;
end
log_L = log_F - 1/2 * sum((res'/C)*res,1);

sim = zeros(300,1);
sim(1:100,:)=d_f(301:400);
sim(101:300)=d_f(1601:1800);
hold on
plot(data_model(301:400));
plot(sim);
hold off
parfor i=1:chain_size
    %calculate the initial cdf for d_PCA 
    f_I_d_PCA = normcdf(d_PCA(:,i),d_prior,d_PCA_sigma.^0.5);
    f_I_d_PCA_c = num2cell(f_I_d_PCA');
    %calculate the transformed value for the given initial cdf of d_PCA
    d_f = cellfun(@(a,b,c) interp1qr(a(2:end),b(2:end),c),f_data,x_data,f_I_d_PCA_c,'UniformOutput', false);
    d_f = reshape(cell2mat(d_f),[Nd 1]);
    %assign resultant vector to matrix
    d_f_m(:,i)=d_f;
end

writematrix(d_f_m,'dsi_results.csv');