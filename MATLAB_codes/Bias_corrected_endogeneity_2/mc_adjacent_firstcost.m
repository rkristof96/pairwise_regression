close all;
clear all;
clc;

global assigned_weight; 
global pairwise_beta1;

% true parameters

alpha = 1;
beta  = 0.5;
epsilon = 0.01;
sigma = 0.5;

b_true = [alpha;beta;sigma];

T = 50; % number of observations
reps = 1000; % number of Monte Carlo repetitions

% explanatory variable
rand('seed',202101);
% generate x: (Tx1) u: (Txreps) vector of bi-variate normal distributed
% random variables

mu = zeros(1,reps+1);
cov_matrix = zeros(reps+1);
cov_matrix(:,:) = 0.65;

cov_matrix(1,:) = sigma;
cov_matrix(:,1) = sigma;

for i=(1:1:reps+1)
    cov_matrix(i,i) = 1;
end

rng('default')  % For reproducibility

R = mvnrnd(mu,cov_matrix,T);

x = R(:,1);
eps = R(:,2:end);
    
% dependent variables, in each of the repetitions

y = alpha+beta*x+eps;  % (T x reps) matrix of dependent variables

% sort
xy = [x y];

%xy = sortrows(xy,1);

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

b_hat_all = zeros(1,reps);  % store estimated betahats, r-th repetition in r-th column

% Calculate d

x_differences = diff(x);
d = median(x_differences);
number_to_keep = 0;
total_deviation = 0;

for i=(1:1:T-1)
    absolute_deviation = abs(x_differences(i)-d);
    if absolute_deviation<epsilon
        number_to_keep = number_to_keep + 1;
        total_deviation = total_deviation + x_differences(i)-d;
    end
end

r = 1;
while r < reps+0.5

    pairwise_betas=zeros(1,number_to_keep);
    delta_x = zeros(1,number_to_keep);
    inverse_delta_y = zeros(1,number_to_keep);
    length = zeros(1,number_to_keep);
    inverse_length = zeros(1,number_to_keep);
    counter = 1;
    
    for i=(1:1:T-1)
        absolute_deviation = abs(x_differences(i)-d);
        if absolute_deviation<epsilon
           pairwise_betas(1,counter) = (y(i+1,r)-y(i,r))/x_differences(i);
           delta_x(1,counter) = x_differences(i);
           inverse_delta_y(1,counter) = 1/(y(i+1,r)-y(i,r));
           length(1,counter) = sqrt(x_differences(i)^2 + (y(i+1,r)-y(i,r))^2);
           inverse_length(1,counter) = 1/(sqrt(x_differences(i)^2 + (y(i+1,r)-y(i,r))^2));
           counter = counter + 1;
        end
    end
    
    assigned_weight = delta_x;
    %assigned_weight = inverse_delta_y;
    %assigned_weight = length;
    %assigned_weight = inverse_length;
        
    pairwise_beta1 = pairwise_betas(1, :);
    
    % Optimization part

    x0             = 5;
    
    [beta1] = fminunc(@lossfunction1_beta1,x0, optimoptions('fminunc','Display','none'));
    
    
    %estimate gamma
    selected_delta_y = zeros(number_to_keep,1);
    selected_delta_x = zeros(number_to_keep,1);
    M = 1;
    
    for i=(1:1:T-1)
        absolute_deviation = abs(x_differences(i)-d);
        if absolute_deviation<epsilon
            selected_delta_y(M,1) = y(i+1,r)-y(i,r);
            selected_delta_x(M,1) = x_differences(i);
            M = M+1;
        end
    end
    
    delta_x_avg     = mean(selected_delta_x);
    delta_y_avg     = mean(selected_delta_y);
    numerator = 0;
    denominator = 0;
    for i=(1:1:number_to_keep)
        x_dev = selected_delta_x(i,1)-delta_x_avg;
        y_dev = selected_delta_y(i,1)-delta_y_avg;
        numerator = numerator + x_dev*y_dev;
        denominator = denominator + x_dev*x_dev;
    end
    b_hat              = numerator/denominator;
    gamma_hat     = delta_y_avg - b_hat*delta_x_avg;
    
    %estimate mu_sx
    
    mu_hat = d + total_deviation/number_to_keep;
    
    %estimate beta_tilde
    
    beta_tilde = (beta1-gamma_hat)*d/mu_hat;    
    
    b_hat_all(1,r)        = beta_tilde;

    r = r + 1;   
    
end
    
standard_dev2=std(b_hat_all(1,:));

%%%%%%%%%%%%
% PRINTING %
%%%%%%%%%%%%

fprintf('\n');
fprintf('\n ADJACENT PAIRWISE ESTIMATION (WITHOUT CONNECTING FIRST AND LAST)\n');
fprintf('Estimated parameters (mean of Monte Carlo repetitions)\n');
fprintf('  Beta:%8.4f\n',mean(b_hat_all(1,:),2));
fprintf('Standard errors (standard deviation of estimates at Monte Carlo repetitions)\n');
fprintf('  Beta:%8.4f',standard_dev2);
fprintf('\n  Number of observations we keep:%8.4f',number_to_keep);
