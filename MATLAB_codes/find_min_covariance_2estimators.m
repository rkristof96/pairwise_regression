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

 dist ='norm'; % or you can write "student-t"
 %uniform distribution is the bottleneck
 %corr(x,u) - correlation for data generation
 %0.8 - 0.6397 - 0.9594
 %0.5 - 0.2493 - 0.3738
 %0.2 - 0.0391 - 0.0585
 %0.1 - 0.0091 - 0.0135
 cov_list = 1;
 
 % parameters:
 mu_x = 1;
 sigma_x = sqrt(2);
 sigma_u = sqrt(1.5);

%rand('seed',20230318);
% generate x: (Tx1) u: (Txreps) vector of bi-variate normal distributed
% random variables
% avarages
mu = zeros(1,reps+1);
mu(1) =mu_x;
cov_matrix = zeros(reps+1);
% reps+1 is needed, first data column includes x values

% iterate over the covariance list
while cov_list > 0
   % cov_matrix's first row + column describes the relationship between x
   % and u
cov_matrix(:,:) = cov_list;
cov_matrix(1,:) =  sigma*sigma_x*sigma_u; %correlation between x and u is x
cov_matrix(:,1) =  sigma*sigma_x*sigma_u; %correlation between x and u is x

% diagonal of the covariance matrix includes ones
cov_matrix(1,1) = sigma_x^2;
for diag_i=(2:1:reps+1)
    cov_matrix(diag_i,diag_i) = sigma_u^2;
end

% case selection for different DGPs
if strcmp(dist,'student-t') % student-t DGP
    %student-t case
    rng('default')  % For reproducibility
    R = mvnrnd(mu,cov_matrix,T); % generate data
    x = R(:,1);
    x = normcdf(x);
    x = tinv(x,5); %transform data to t-distribution data
end
if strcmp(dist,'norm') % normal DGP
    % normal DGP
    %cov_matrix(1,1) = sigma_x^2; %variance is 5
    
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
