function [log_L] = MODEL2(y)

persistent data_model nobs d_prior Nd max_sv tphi d_PCA_sigma f_data x_data invC log_F

if isempty(data_model)
    data_model = readmatrix("obs_data.csv");
    nobs = size(data_model,1);
end

if isempty(max_sv)
    %define the level of energy for TSVD
    le = 0.995;
    %load realizations
    d_full = readmatrix("cal1a_ies.0.obs.csv");
    %identify ensemble size (Ne) and number of data points(Nd)
    Ne = size(d_full,2);
    Nd = size(d_full,1);
    %calculate the mean for each data point
    d_prior = repmat(mean(d_full,2),1,Ne);
    %perform SVD
    X = (d_full-d_prior)/((Ne-1)^0.5);
    [U,S,~]=svd(X);
    sv_squared = diag(S).^2;
    total_energy = sum(sv_squared);
    cum_energy = cumsum(sv_squared);

    %calculate the number of svs
    max_sv = find(cum_energy >= le * total_energy, 1);

    %calculate truncated phi (tphi) and d_PCA standard deviations
    %(d_PCA_sigma)
    phi = U*S;
    tphi = phi(1:end,1:max_sv);
    d_PCA_cov = tphi*tphi';
    d_PCA_sigma = diag(d_PCA_cov).^0.5;

    % Compute the empirical CDF for each observation
    [f_data,x_data] = arrayfun(@(row) ecdf(d_full(row,:)), 1:Nd, 'UniformOutput', false);
end

% Specify the covariance matrix --> do only once
if isempty(invC)
    
    % How many dimensions?
    %d = size(x,2);
    
    % Target covariance
    C  = (0.01)^(2)*eye(nobs,nobs);
    invC = inv(C);
  
    % Calculate integration constant
    log_F = log ( ( ( 2 * pi )^( - nobs / 2 ) ) * det(C)^( - 1 / 2 ) );
    
    % log_F becomes -Inf for large d --> hence we set log_F to zero
    % Also need to resolve this in importance distribution as well!!
    if ( nobs > 150 )
        log_F = 0;
    end

end

%calculate d_PCA with standard variates y
d_PCA = d_prior(:,1)+tphi*y';
%calculate the cum prob
f_I_d_PCA = normcdf(d_PCA,d_prior(:,1),d_PCA_sigma);
f_I_d_PCA_c = num2cell(f_I_d_PCA');

d_f = cellfun(@(a,b,c) interp1qr(a(1:end),b(1:end),c),f_data,x_data,f_I_d_PCA_c,'UniformOutput', false);
d_f_m = reshape(cell2mat(d_f),[Nd 1]);

res = d_f_m(1:nobs,1)-data_model;

log_L = log_F - 1/2 * sum(res'*invC*res,1);

%log_L = sum(log(normpdf(results_model,data_model,0.01*ones(nobs))));