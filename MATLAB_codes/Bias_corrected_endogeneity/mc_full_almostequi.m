close all;
clear all;
clc;

% true parameters

alpha = 1;
beta  = 0.5;
epsilon = 0.01;
sigma = 0.5;

b_true = [alpha;beta;sigma];

T = 500; % number of observations
reps = 1000; % number of Monte Carlo repetitions

%%%%%%%%%%%%%%%%%%%
% DATA GENERATION %
%%%%%%%%%%%%%%%%%%%

% explanatory variable
rand('seed',202101);
% generate x: (Txreps) u: (Txreps) vector of bi-variate normal distributed
% random variables

x = zeros(T,reps);
eps = zeros(T,reps);

mu = zeros(1,2);
cov_matrix = zeros(2);

cov_matrix(1,1) = 1;
cov_matrix(1,2) = sigma;
cov_matrix(2,1) = sigma;
cov_matrix(2,2) = 1;

rng('default')  % For reproducibility
for i=(1:1:reps)
    R = mvnrnd(mu,cov_matrix,T);
    x(:,i) = R(:,1);
    eps(:,i) = R(:,2);
end

% generate the dependent variable
y = alpha+beta*x+eps;

% sort
x_sorted = zeros(T,reps);
y_sorted = zeros(T,reps);

for i=(1:1:reps)
    x_to_sort = x(:,i);
    y_to_sort = y(:,i);
    xy = [x_to_sort y_to_sort];
    xy = sortrows(xy,1);
    x_sorted(:,i) = xy(:,1);
    y_sorted(:,i) = xy(:,2);
end

%x = x_sorted;
%y = y_sorted;

%%%%%%%%%%%%%%
% OLS ESTIMATION %
%%%%%%%%%%%%%%

b_hat_all = zeros(2,reps);  % store estimated betahats, r-th repetition in r-th column

r = 1;
while r < reps+0.5
    x_mc = x(:,r);
    x_avg     = mean(x_mc);
    y_avg     = mean(y);
    y_avg_r   = y_avg(r);
    numerator = 0;
    denominator = 0;
    for i=(1:1:T)
        x_dev = x_mc(i,1)-x_avg;
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
% Full PAIRWISE ESTIMATION %
%%%%%%%%%%%%%%

b_hat_all = zeros(1,reps);  % store estimated betahats, r-th repetition in r-th column

r = 1;
while r < reps+0.5
    number_of_pairs = T * (T-1) /2;
    delta_x = zeros(1,number_of_pairs);
    delta_y = zeros(1,number_of_pairs);
    counter = 1;
    
    % iterate over all pairs
    for i=(1:1:T-1)
        for j=(2:1:T)
            if i<j
                % calculate x difference
                delta_x(1, counter) = x(j,r) - x(i,r);
                delta_y(1, counter) = y(j,r) - y(i,r);
                counter = counter + 1;
            end
        end
    end
    
    % Calculate d
    
    d = median(delta_x(:,r));
    total_deviation = 0;
    
    sum_delta_y = 0;
    N = 0;
    for i=(1:1:number_of_pairs)
        absolute_deviation = abs(delta_x(1, i)-d);
        if absolute_deviation<epsilon
            total_deviation = total_deviation + delta_x(1,i)-d;
            sum_delta_y = sum_delta_y + delta_y(1, i);
            N = N+1;
        end
    end
    N
    
    %estimate beta
    beta_hat = sum_delta_y/(N*d);
    beta_hat_corrected = ((1+ total_deviation/(d*N))^(-1))*beta_hat;
    
    b_hat_all(1,r)        = beta_hat_corrected;
    
    r = r + 1;
end

standard_dev1=std(b_hat_all(1,:));

%%%%%%%%%%%%
% PRINTING %
%%%%%%%%%%%%

fprintf('\n');
fprintf('\n FULL PAIRWISE ESTIMATION\n');
fprintf('Estimated parameters (mean of Monte Carlo repetitions)\n');
fprintf('Beta:%8.4f',mean(b_hat_all(1,:),2));
fprintf('\n Standard errors (standard deviation of estimates at Monte Carlo repetitions)\n');
fprintf('  Beta:%8.4f',standard_dev1);
fprintf('\n  Number of observations we keep:%8.4f',N);