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
rand('seed',202101);
% generate x: (Tx1) vector of uniformly distributed random
%    variables on the interval (-1;+1) 
x = rand(T,1)*2-1;
x_sorted = sort(x);

% error terms
% generate i.i.d. error terms, each normally distributed with 0
%    expected value and sigma^2 variance
randn('seed',202102);
eps = normrnd(0,sigma, [T,1]);

% generate the dependent variable
y = alpha+beta*x_sorted+eps;
x = x_sorted;

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
% SORTED PAIRWISE ESTIMATION (WITHOUT CONNECTING FIRST AND LAST)  %
%%%%%%%%%%%%%%
pairwise_betas_1=zeros(2,T-1);
% iterate over all pairs
for i=(1:1:T-1)
    % calculate betahat
    x_avg     = mean(x(i:i+1));
    y_avg     = mean(y(i:i+1));
    numerator = y(i+1,1) - y(i,1);
    denominator = x(i+1,1) - x(i,1);
    b_hat     = numerator/denominator;
    alpha_hat = mean(y(i:i+1)) - b_hat*mean(x(i:i+1));
    pairwise_betas_1(1,i)=alpha_hat;
    pairwise_betas_1(2,i)=b_hat;
end;

%display average of betas
average_parwise_betas_1=mean(pairwise_betas_1,2);
fprintf('\n');
fprintf('SORTED PAIRWISE ESTIMATION (WITHOUT CONNECTING FIRST AND LAST)\n');
fprintf('Alpha:%8.4f',average_parwise_betas_1(1,:));
fprintf('  Beta:%8.4f\n',average_parwise_betas_1(2,:));


%%%%%%%%%%%%%%
% SORTED PAIRWISE ESTIMATION (WITH CONNECTING FIRST AND LAST)  %
%%%%%%%%%%%%%%
pairwise_betas_2=zeros(2,T);
% iterate over all pairs
for i=(1:1:T)
    if i~=T
        % calculate betahat
        x_avg     = mean(x(i:i+1));
        y_avg     = mean(y(i:i+1));
        numerator = y(i+1,1) - y(i,1);
        denominator = x(i+1,1) - x(i,1);
        b_hat     = numerator/denominator;
        alpha_hat = mean(y(i:i+1)) - b_hat*mean(x(i:i+1));
        pairwise_betas_2(1,i)=alpha_hat;
        pairwise_betas_2(2,i)=b_hat;
    else
        x_first_and_last = [x(1); x(T)];
        y_first_and_last = [y(1); y(T)];
        % calculate betahat
        x_avg     = mean(x_first_and_last);
        y_avg     = mean(y_first_and_last);
        numerator = y(i,1) - y(1,1);
        denominator = x(i,1) - x(1,1);
        b_hat     = numerator/denominator;
        alpha_hat = mean(y_first_and_last) - b_hat*mean(x_first_and_last);
        pairwise_betas_2(1,i)=alpha_hat;
        pairwise_betas_2(2,i)=b_hat;
    end;
    % calculate (Tx1) vector of residuals
    %resid     = y-x_matr*b_hat;
    % estimate sigma: (root of) sum of squared residuals divided by degrees of freedom
    %sigma_hat = sqrt(sum(resid'*resid)/(T-2));
    %b_hat     = [b_hat; sigma_hat];
    % calculate estimated standard errors: this is (sigma^2*(X'X)^(-1))^0.5
    %stderr    = (sigma_hat*sigma_hat*inv(x_matr'*x_matr)).^0.5;
end;

%display average of betas
average_parwise_betas_2=mean(pairwise_betas_2,2);
fprintf('\n');
fprintf('SORTED PAIRWISE ESTIMATION (WITH CONNECTING FIRST AND LAST)\n');
fprintf('Alpha:%8.4f',average_parwise_betas_2(1,:));
fprintf('  Beta:%8.4f\n',average_parwise_betas_2(2,:));
