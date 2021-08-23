close all;
clear all;
clc;

% true parameters

alpha = 1;
beta  = 1.5;
sigma = sqrt(0.5);
b_true = [alpha;beta;sigma];

T = 500; % number of observations
T_sub = 0.6 * T;
reps = 1000; % number of Monte Carlo repetitions

% explanatory variable
x = [1:T]';

s = RandStream('mlfg6331_64'); 
x = datasample(s,x,T_sub,'Replace',false);

T = T_sub;

% error terms

randn('seed',202101);
eps = normrnd(0,sigma, [T,reps]);  %generate (T x reps) matrix of normally distributed i.i.d. errors,
    %with mean 0 and variance sigma^2
    
x_mean = mean(x);
x_std = std(x);

x_standard = (x-x_mean)/(x_std*sqrt(2));

eps_endog = eps + 10 * x_standard;

% Make epsilon endogenous
%eps = eps_endog;

% random shuffle
x_and_eps = [x eps];

s = RandStream('mlfg6331_64');
random_x_and_eps = x_and_eps(randperm(s,size(x_and_eps, 1)), :);

x = random_x_and_eps(:,1);
eps = random_x_and_eps(:,2:reps+1);

% dependent variables, in each of the repetitions

y = alpha+beta*x+eps;  % (T x reps) matrix of dependent variables

% sort
x_y_eps = [x y eps];

x_y_eps = sortrows(x_y_eps,1);

x = x_y_eps(:,1);
y = x_y_eps(:,2:reps+1);
eps = x_y_eps(:,reps+2:end);

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
%fprintf('Se(a):%8.4f',var_true(1,1)^0.5);
%fprintf('  Se(b):%7.4f\n',var_true(2,2)^0.5);
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
% PAIRWISE ESTIMATION %
%%%%%%%%%%%%%%

b_hat_all = zeros(2,reps);  % store estimated betahats, r-th repetition in r-th column

r = 1;
while r < reps+0.5
    number_of_betas = T * (T-1) /2;
    pairwise_betas=zeros(2,number_of_betas);
    delta_y = zeros(1,number_of_betas);
    delta_x = zeros(1,number_of_betas);
    counter=1;

    % iterate over all pairs
    for i=(1:1:T-1)
        for j=(2:1:T)
            if i<j
                % calculate y difference
                delta_y(1, counter) = y(j,r) - y(i,r);
                delta_x(1, counter) = x(j,1) - x(i,1);                
                % calculate betahat
                x_avg     = (x(i,1)+x(j,1))/2;
                y_avg     = (y(i,r)+y(j,r))/2;          
                numerator = y(j,r) - y(i,r);
                denominator = x(j,1) - x(i,1);
                b_hat_i     = numerator/denominator;
                alpha_hat_i = y_avg - b_hat_i*x_avg;
                pairwise_betas(1,counter)=alpha_hat_i;
                pairwise_betas(2,counter)=b_hat_i;
                counter   = counter+1;
            end;
        end;
    end;
    
    %delta_y = 1./delta_y;
    %delta_y = abs(delta_y);
    %delta_x = 1./delta_x;
    %delta_x = abs(delta_x);
    weighting_delta = delta_y;
    sum_weighting_delta = sum(weighting_delta);
    weighted_parwise_betas = pairwise_betas*weighting_delta';
    weighted_average_parwise_betas = weighted_parwise_betas./sum_weighting_delta;
   
    
    %average_parwise_betas = mean(pairwise_betas,2);
    
    %b_hat_all(1,r)        = average_parwise_betas(1,:);
    %b_hat_all(2,r)        = average_parwise_betas(2,:);

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
fprintf('\nFULL-PAIRWISE PAIRWISE ESTIMATION\n');
fprintf('Estimated parameters (mean of Monte Carlo repetitions)\n');
fprintf('Alpha:%8.4f',mean(b_hat_all(1,:),2));
fprintf('  Beta:%8.4f\n',mean(b_hat_all(2,:),2));
fprintf('Standard errors (standard deviation of estimates at Monte Carlo repetitions)\n');
fprintf('Alpha:%8.4f',standard_dev1);
fprintf('  Beta:%8.4f',standard_dev2);
