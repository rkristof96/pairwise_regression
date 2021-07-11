close all;
clear all;
clc;

% true parameters

alpha = 1;
beta  = 0.5;
epsilon = 0.01;
sigma = 0.5;

b_true = [alpha;beta;sigma];

T = 5000; % number of observations
reps = 1000; % number of Monte Carlo repetitions

%%%%%%%%%%%%%%%%%%%
% DATA GENERATION %
%%%%%%%%%%%%%%%%%%%

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

% generate the dependent variable
y = alpha+beta*x+eps;

% sort
xy = [x y];

%xy = sortrows(xy,1);

x = xy(:,1);
y = xy(:,2:reps+1);

%%%%%%%%%%%%%%
% OLS ESTIMATION %
%%%%%%%%%%%%%%

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
% ADJACENT ALMOST EQUI PAIRWISE ESTIMATION (WITHOUT CONNECTING FIRST AND LAST)  %
%%%%%%%%%%%%%%
b_hat_all = zeros(1,reps);  % store estimated betahats, r-th repetition in r-th column

% Calculate d_1

x_differences = diff(x);
abs_x_differences = abs(x_differences);
d = median(abs_x_differences);

r = 1;
while r < reps+0.5 
    total_deviation = 0;
    
    sum_betas = 0;
    delta_x_sum_betas = 0;
    inverse_delta_y_sum_betas = 0;
    length_sum_betas = 0;
    inverse_length_sum_betas = 0;
    N = 0;

    inverse_abs_delta_y_sum_betas = 0;
    abs_delta_x_sum_betas = 0;
    
    for i=(1:1:T-1)
        absolute_deviation = abs(abs_x_differences(i)-d);
        if absolute_deviation<epsilon
            total_deviation = total_deviation + abs_x_differences(i)-d;

            sum_betas = sum_betas + (y(i+1,r)-y(i,r))/x_differences(i);
            delta_x_sum_betas = delta_x_sum_betas + (y(i+1,r)-y(i,r));
            inverse_delta_y_sum_betas = inverse_delta_y_sum_betas + 1/x_differences(i);
            length = sqrt(x_differences(i)^2 + (y(i+1,r)-y(i,r))^2);
            length_sum_betas = length_sum_betas + length * (y(i+1,r)-y(i,r))/x_differences(i);
            inverse_length_sum_betas = inverse_length_sum_betas + (1/length) * (y(i+1,r)-y(i,r))/x_differences(i);     
            N = N+1;
            
            inverse_abs_delta_y_sum_betas = inverse_abs_delta_y_sum_betas + (1/abs(y(i+1,r)-y(i,r))) * (y(i+1,r)-y(i,r))/x_differences(i);
            abs_delta_x_sum_betas = abs_delta_x_sum_betas + abs(x_differences(i))*(y(i+1,r)-y(i,r))/x_differences(i);
        end
    end
    
    %estimate beta
    %beta_hat = sum_betas/(N);
    %beta_hat = delta_x_sum_betas/N;
    %beta_hat = inverse_delta_y_sum_betas/N;
    %beta_hat = length_sum_betas/N;
    beta_hat = inverse_length_sum_betas/N;

    %beta_hat = inverse_abs_delta_y_sum_betas/N;
    %beta_hat = abs_delta_x_sum_betas/N;
   
    
    %estimate gamma
    selected_delta_y = zeros(N,1);
    selected_delta_x = zeros(N,1);
    M = 1;
    
    for i=(1:1:T-1)
        absolute_deviation = abs(abs_x_differences(i)-d);
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
    for i=(1:1:N)
        x_dev = selected_delta_x(i,1)-delta_x_avg;
        y_dev = selected_delta_y(i,1)-delta_y_avg;
        numerator = numerator + x_dev*y_dev;
        denominator = denominator + x_dev*x_dev;
    end
    b_hat              = numerator/denominator;
    gamma_hat     = delta_y_avg - b_hat*delta_x_avg;
    
    %estimate mu_sx
    
    mu_hat = d + total_deviation/N;
    
    %estimate beta_tilde
    
    beta_tilde = (beta_hat-gamma_hat)*d/mu_hat;
   
    
    b_hat_all(1,r)        = beta_tilde;

    r = r + 1;
end

standard_dev1=std(b_hat_all(1,:));

%%%%%%%%%%%%
% PRINTING %
%%%%%%%%%%%%

fprintf('\n');
fprintf('\n ADJACENT ALMOST EQUI PAIRWISE ESTIMATION (WITHOUT CONNECTING FIRST AND LAST)\n');
fprintf('Estimated parameters (mean of Monte Carlo repetitions)\n');
fprintf('  Beta:%8.4f\n',mean(b_hat_all(1,:),2));
fprintf('Standard errors (standard deviation of estimates at Monte Carlo repetitions)\n');
fprintf('  Beta:%8.4f',standard_dev1);
fprintf('\n  Number of observations we keep:%8.4f',N);
fprintf('\n  The d we set:%8.4f',d);
