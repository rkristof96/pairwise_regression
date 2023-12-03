%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code serves to calculate the appropriate correlation values to
% generate correlated x and u values (minimum correlation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% you can play with the parameters one-by-one
% T = 50 is enough to determine the correlation values
T = 50; % number of observations
reps = 1000; % number of Monte Carlo repetitions

 sigma = 0.8;
 % sigma = 0.1
 %cov = 0
 % this values are appropritae for the normal case, but not for the uniform
 %[0,0.2,0.5,0.8];
 %[0,0.0071,0.0491,0.1272]; 
 %sigma
 %0,0.2,0.5,0.8

 dist ='uniform'; % or you can write "norm"
 %uniform distribution is the bottleneck
 %corr(x,u) - correlation for data generation
 %0.8 - 0.6397
 %0.5 - 0.2493
 %0.2 - 0.0391
 %0.1 - 0.0091
 cov_list = 0.8;

%rand('seed',20230318);
% generate x: (Tx1) u: (Txreps) vector of bi-variate normal distributed
% random variables
% avarages
mu = zeros(1,reps+1);
cov_matrix = zeros(reps+1);
% reps+1 is needed, first data column includes x values

% iterate over the covariance list
while cov_list > 0
   % cov_matrix's first row + column describes the relationship between x
   % and u
cov_matrix(:,:) = cov_list;
cov_matrix(1,:) = sigma; %correlation between x and u is x
cov_matrix(:,1) = sigma; %correlation between x and u is x

for diag_i=(1:1:reps+1)
    cov_matrix(diag_i,diag_i) = 1; % diagonal correlations are 1
end

% case selection for different DGPs
if strcmp(dist,'uniform') % uniform DGP
    % generate correlated values
    R = mvnrnd(mu,cov_matrix,T);
    x = R(:,1);
    x = normcdf(x);
    x = unifinv(x,-5,5); % U(-5,5)
else
    % normal DGP
    cov_matrix(1,1) = 5; %variance is 5
    
    R = mvnrnd(mu,cov_matrix,T);
    x = R(:,1);
end

% correlation values
cov_list = cov_list - 0.0001;

end


% you have to print out it manually
% code runs until it reaches an error (non-positive semi-definite
% covariance matrix)
% last value of covariance
cov_list
