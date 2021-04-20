close all;
clear all;
clc;

global assigned_weight; 
global pairwise_beta1;

% true parameters


alpha = 1;
beta  = 0.5;
epsilon = 0.01;
sigma = 0.8;

b_true = [alpha;beta;sigma];

T = 5000; % number of observations
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
% FULL PAIRWISE PAIRWISE ESTIMATION %
%%%%%%%%%%%%%%

b_hat_all = zeros(1,reps);  % store estimated betahats, r-th repetition in r-th column

number_of_pairs = T * (T-1) /2;
delta_x = zeros(1,number_of_pairs);
counter = 1;

% iterate over all pairs
for i=(1:1:T-1)
    for j=(2:1:T)
        if i<j
            % calculate x difference
            delta_x(1, counter) = x(j,1) - x(i,1);
            counter = counter + 1;
        end
    end
end 

d_1 = median(delta_x);

number_to_keep = 0;

for i=(1:1:number_of_pairs)
    absolute_deviation = abs(delta_x(i)-d_1);
    if absolute_deviation<epsilon
        number_to_keep = number_to_keep + 1;
    end
end

r = 1;
while r < reps+0.5
    delta_y = zeros(1,number_of_pairs);
    counter = 1;

    % iterate over all pairs
    for i=(1:1:T-1)
        for j=(2:1:T)
            if i<j
                % calculate x difference
                delta_y(1, counter) = y(j,r) - y(i,r);
                counter = counter + 1;
            end
        end
    end    

    pairwise_betas=zeros(1,number_to_keep);
    delta_x_weight = zeros(1,number_to_keep);
    inverse_delta_y = zeros(1,number_to_keep);
    length = zeros(1,number_to_keep);
    inverse_length = zeros(1,number_to_keep);
    counter = 1;
    
    % iterate over all pairs
    for i=(1:1:number_of_pairs)
        absolute_deviation = abs(delta_x(i)-d_1);
        if absolute_deviation<epsilon
           pairwise_betas(1,counter) = delta_y(i)/delta_x(i);
           delta_x_weight(1,counter) = delta_x(i);
           inverse_delta_y(1,counter) = 1/delta_y(i);
           length(1,counter) = sqrt(delta_x(i)^2 + delta_y(i)^2);
           inverse_length(1,counter) = 1/(sqrt(delta_x(i)^2 + delta_y(i)^2));
           counter = counter + 1;
        end
    end
    
    %assigned_weight = delta_x_weight;
    %assigned_weight = inverse_delta_y;
    %assigned_weight = length;
    assigned_weight = inverse_length;
        
    pairwise_beta1 = pairwise_betas(1, :);
    
    % Optimization part

    x0             = 5;
    
    [beta1] = fminunc(@lossfunction1_beta1,x0, optimoptions('fminunc','Display','none'));
    
    b_hat_all(1,r)        = beta1;
    

    r = r + 1;
end

standard_dev2=std(b_hat_all(1,:));

%%%%%%%%%%%%
% PRINTING %
%%%%%%%%%%%%

fprintf('\n');
fprintf('\n FULL PAIRWISE PAIRWISE ESTIMATION\n');
fprintf('Estimated parameters (mean of Monte Carlo repetitions)\n');
fprintf('  Beta:%8.4f\n',mean(b_hat_all(1,:),2));
fprintf('Standard errors (standard deviation of estimates at Monte Carlo repetitions)\n');
fprintf('  Beta:%8.4f',standard_dev2);
fprintf('\n  Number of observations we keep:%8.4f',number_to_keep);

