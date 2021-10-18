close all;
clear all;
clc;

% true parameters

alpha = 0;
beta  = 1.5;
sigma = sqrt(10);
a = 0;
b = 1;
integers_to_keep = [19 20];
b_true = [alpha;beta;sigma];

T = 5000; % number of observations
reps = 1000; % number of Monte Carlo repetitions

% explanatory variable
rng('default');
i = unidrnd(20,T,1);

x = a + b * i;

% error terms

randn('seed',202101);
eps = normrnd(0,sigma, [T,reps]);  %generate (T x reps) matrix of normally distributed i.i.d. errors,
    %with mean 0 and variance sigma^2
    
eps_endog = x + eps;

% Make epsilon endogenous
eps = eps_endog;

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

counter = 0;

for i=(1:1:T)
    if any(x(i,1) == integers_to_keep)
        counter = counter + 1;
    end
end

x_kept = zeros(counter, 1);
y_kept = zeros(counter, reps);

holder = 1;

for i=(1:1:T)
    if any(x(i,1) == integers_to_keep)
        x_kept(holder, 1) = x(i,1);
        y_kept(holder, :) = y(i, :);
        holder = holder +1;
    end
end

x = x_kept;
y = y_kept;

%x_matr = [ones(T,1) x];  % this is matrix X in betahat = (X'X)^(-1)*(X'y)
%var_true = (sigma^2)*inv(x_matr'*x_matr);  % true variance-covariance matrix

b_hat_all = zeros(1,reps);  % store estimated betahats, r-th repetition in r-th column

r = 1;
while r < reps+0.5
    numerator = x(:, 1)' * y(:,r);
    denominator = x(:, 1)' * x(:, 1);
    
    b_hat              = numerator/denominator;
    b_hat_all(1,r)     = b_hat;

    r = r + 1;
end

standard_dev1=std(b_hat_all(1,:));

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
fprintf('  Beta:%8.4f\n',mean(b_hat_all(1,:),2));
fprintf('Standard errors (standard deviation of estimates at Monte Carlo repetitions)\n');
fprintf('  Beta:%8.4f\n',standard_dev1);
fprintf('Number of observations we keep\n');
fprintf('  N:%8.4f',counter);