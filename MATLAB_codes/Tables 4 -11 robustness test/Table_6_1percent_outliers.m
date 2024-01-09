close all;
clear all;
clc;

% true parameters


alpha = 1;
beta  = 0.5;
sigma = 1;
xi = -sqrt(2/pi);

b_true = [alpha;beta;sigma];

T = 5000; % number of observations
reps = 1000; % number of Monte Carlo repetitions

% explanatory variable
rand('seed',202101);
% generate x: (Tx1) vector of uniformly distributed random
%    variables on the interval (-1;+1) 
%x = rand(T,1)*2-1;
x = rand(T,1)*20-10;
%x = normrnd(0,5, [T,1]);

%Z = normrnd(0,1, [T,1]);
%tau = abs(Z);
%rand('seed',202020);
%U = normrnd(0,1, [T,1]);
%x = xi + tau + U;

% error terms

randn('seed',202101);
%eps = normrnd(0,sigma, [T,reps]);  %generate (T x reps) matrix of normally distributed i.i.d. errors,
    %with mean 0 and variance sigma^2

Z_eps = normrnd(0,1, [T,reps]);
tau_eps = abs(Z_eps);
rand('seed',222022);
U_eps = normrnd(0,1, [T,reps]);
% generate outliers
% 1% of the error terms are multiplied by 10
U_eps_outliers = U_eps;
num_outliers = ceil(T/100);

for i = 1:reps
   [Bsort Bidx] = getNElements(U_eps(:,i), num_outliers);
   %[Bsort Bidx] = getNElements_abs(U_eps(:,i), num_outliers);
   U_eps_outliers(Bidx,i) = U_eps_outliers(Bidx,i)*10;
end

eps = xi*0 + tau_eps*0 + U_eps_outliers;
    
% dependent variables, in each of the repetitions

y = alpha+beta*x+eps;  % (T x reps) matrix of dependent variables

% sort
xy = [x y];

xy = sortrows(xy,1);

x = xy(:,1);
y = xy(:,2:reps+1);


%%%%%%%%%%%%%%
% OLS ESTIMATION %
%%%%%%%%%%%%%%

x_matr = [ones(T,1) x];  % this is matrix X in betahat = (X'X)^(-1)*(X'y)
var_true = (sigma^2)*inv(x_matr'*x_matr);  % true variance-covariance matrix

b_hat_all = zeros(2,reps);  % store estimated betahats, r-th repetition in r-th column

r = 1;
while r < reps+0.5
    x_avg     = mean(x);
    y_avg     = mean(y);
    y_avg_r   = y_avg(r);
    numerator = 0;
    denominator = 0;
    for i=(1:1:T)
        x_dev = x(i,1)-x_avg;
        y_dev = y(i,r)-y_avg_r;
        numerator = numerator + x_dev*y_dev;
        denominator = denominator + x_dev*x_dev;
    end;
    b_hat              = numerator/denominator;
    b_hat_all(2,r)     = b_hat;
    b_hat_all(1,r)     = y_avg_r - b_hat*x_avg;

    r = r + 1;
end

standard_dev1=std(b_hat_all(1,:));

standard_dev2=std(b_hat_all(2,:));

%%%%%%%%%%%%
% PRINTING %
%%%%%%%%%%%%

fprintf('\nTrue parameters\n');
fprintf('Alpha:%8.4f',b_true(1));
fprintf('  Beta:%8.4f',b_true(2));
fprintf('  Sigma:%8.4f\n',b_true(3));
fprintf('Se(a):%8.4f',var_true(1,1)^0.5);
fprintf('  Se(b):%7.4f\n',var_true(2,2)^0.5);
% print your results: means across Monte-Carlo repetitions
fprintf('\n');
fprintf('OLS Estimation\n');
fprintf('Estimated parameters (mean of Monte Carlo repetitions)\n');
fprintf('Alpha:%8.4f',mean(b_hat_all(1,:),2));
fprintf('  Beta:%8.4f\n',mean(b_hat_all(2,:),2));
fprintf('Standard errors (standard deviation of estimates at Monte Carlo repetitions)\n');
fprintf('Alpha:%8.4f',standard_dev1);
fprintf('  Beta:%8.4f',standard_dev2);


%%%%%%%%%%%%%%
% ADJACENT PAIRWISE ESTIMATION (WITHOUT CONNECTING FIRST AND LAST) %
%%%%%%%%%%%%%%

b_hat_all = zeros(2,reps);  % store estimated betahats, r-th repetition in r-th column

r = 1;
while r < reps+0.5
    pairwise_betas=zeros(2,T-1);
    
    for i=(1:1:T-1)
        % calculate betahat
        x_avg     = mean(x(i:i+1));
        y_avg     = mean(y(i:i+1, r));
        numerator = y(i+1,r) - y(i,r);
        denominator = x(i+1,1) - x(i,1);
        b_hat_i     = numerator/denominator;
        alpha_hat_i = y_avg - b_hat_i*x_avg;
        pairwise_betas(1,i)=alpha_hat_i;
        pairwise_betas(2,i)=b_hat_i;
    end
    % Obtain the delta-x weighted average of pairwise betas
    
    delta_x = diff(x);
    delta_x = abs(delta_x);
    sum_delta_x = sum(delta_x);
    weighted_parwise_betas = pairwise_betas*delta_x;
    weighted_average_parwise_betas = weighted_parwise_betas./sum_delta_x;
    %weighted_average_parwise_betas = weighted_parwise_betas./(T-1);
    
    %inv_delta_x = 1./delta_x;
    %sum_inv_delta_x = sum(inv_delta_x);
    %weighted_parwise_betas = pairwise_betas*inv_delta_x;
    %weighted_average_parwise_betas = weighted_parwise_betas./sum_inv_delta_x;
    %weighted_average_parwise_betas = weighted_parwise_betas./(T-1);
    
    % Simple average for beta_0
    %pairwise_betas = sum(pairwise_betas,2)./(T-1);
    %b_hat_all(1,r)        = pairwise_betas(1);
    
    b_hat_all(1,r)        = weighted_average_parwise_betas(1);
    b_hat_all(2,r)        = weighted_average_parwise_betas(2);

    r = r + 1;   
    
end
    
standard_dev1=std(b_hat_all(1,:));

standard_dev2=std(b_hat_all(2,:));

%%%%%%%%%%%%
% PRINTING %
%%%%%%%%%%%%

fprintf('\n');
fprintf('\n ADJACENT PAIRWISE ESTIMATION (WITHOUT CONNECTING FIRST AND LAST)\n');
fprintf('Estimated parameters (mean of Monte Carlo repetitions)\n');
fprintf('Alpha:%8.4f',mean(b_hat_all(1,:),2));
fprintf('  Beta:%8.4f\n',mean(b_hat_all(2,:),2));
fprintf('Standard errors (standard deviation of estimates at Monte Carlo repetitions)\n');
fprintf('Alpha:%8.4f',standard_dev1);
fprintf('  Beta:%8.4f',standard_dev2);
