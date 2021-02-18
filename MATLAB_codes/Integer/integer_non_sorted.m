close all;
clear all;
clc;

% true parameters

alpha = 1;
beta  = 0.5;
sigma = 1;

b_true = [alpha;beta;sigma];

T = 500; % number of observations

%%%%%%%%%%%%%%%%%%%
% DATA GENERATION %
%%%%%%%%%%%%%%%%%%%

% explanatory variable
% rand('seed',202101);
% generate x: (Tx1) vector of integers
x = [1:T]';
x_sorted=x;

% error terms
% generate i.i.d. error terms, each normally distributed with 0
%    expected value and sigma^2 variance
randn('seed',202102);
eps = normrnd(0,sigma, [T,1]);

% generate the dependent variable
y = alpha+beta*x_sorted+eps;

%%%%%%%%%%%%%%
% OLS ESTIMATION %
%%%%%%%%%%%%%%

% calculate betahat
x_avg     = mean(x);
y_avg     = mean(y);
numerator = 0;
denominator = 0;
for i=(1:1:T)
    x_dev = x(i,1)-x_avg;
    y_dev = y(i,1)-y_avg;
    numerator = numerator + x_dev*y_dev;
    denominator = denominator + x_dev*x_dev;
end;
b_hat     = numerator/denominator;
alpha_hat = mean(y) - b_hat*mean(x);

%%%%%%%%%%%%
% PRINTING %
%%%%%%%%%%%%

fprintf('True parameters\n');
fprintf('Alpha:%8.4f',b_true(1));
fprintf('  Beta:%8.4f',b_true(2));
fprintf('  Sigma:%8.4f\n',b_true(3));
fprintf('\n');

fprintf('OLS Estimation\n');
fprintf('Estimated parameters\n');
fprintf('Alpha:%8.4f',alpha_hat);
fprintf('  Beta:%8.4f',b_hat);
fprintf('\n');

%%%%%%%%%%%%%%
% NON-SORTED PAIRWISE ESTIMATION %
%%%%%%%%%%%%%%
number_of_betas = T * (T-1) /2;
pairwise_betas_2=zeros(2,number_of_betas);
counter=1;
% iterate over all pairs
for i=(1:1:T-1)
    for j=(2:1:T)
        if i<j
            % calculate betahat
            x_avg     = (x(i,1)+x(j,1))/2;
            y_avg     = (y(i,1)+y(j,1))/2;
            x_dev_i   = x(i,1)-x_avg;
            y_dev_i   = y(i,1)-y_avg;
            x_dev_j   = x(j,1)-x_avg;
            y_dev_j   = y(j,1)-y_avg;           
            numerator = x_dev_i*y_dev_i + x_dev_j*y_dev_j;
            denominator = x_dev_i*x_dev_i + x_dev_j*x_dev_j;
            b_hat_i     = numerator/denominator;
            alpha_hat_i = y_avg - b_hat*y_avg;
            pairwise_betas_2(1,counter)=alpha_hat_i;
            pairwise_betas_2(2,counter)=b_hat_i;
            counter   = counter+1;
        end;
    end;
end;

%display average of betas
average_parwise_betas_2=mean(pairwise_betas_2,2);
fprintf('\n');
fprintf('NON-SORTED INTEGER PAIRWISE ESTIMATION \n');
fprintf('Alpha:%8.4f',average_parwise_betas_2(1,:));
fprintf('  Beta:%8.4f\n',average_parwise_betas_2(2,:));
