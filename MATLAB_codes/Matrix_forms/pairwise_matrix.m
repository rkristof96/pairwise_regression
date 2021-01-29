close all;
clear all;
clc;

% true parameters

alpha = 1;
beta  = 0.5;
sigma = 1;

b_true = [alpha;beta;sigma];

T = 100; % number of observations

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

% x_matr is the X in betahat = (X'X)^(-1)(X'y) expression
x_matr = [ones(T,1), x];

% calculate betahat
b_hat     = inv(x_matr'*x_matr)*x_matr'*y;
% calculate (Tx1) vector of residuals
%resid     = y-x_matr*b_hat;
% estimate sigma: (root of) sum of squared residuals divided by degrees of freedom
%sigma_hat = sqrt(sum(resid'*resid)/(T-2));
%b_hat     = [b_hat; sigma_hat];
% calculate estimated standard errors: this is (sigma^2*(X'X)^(-1))^0.5
%stderr    = (sigma_hat*sigma_hat*inv(x_matr'*x_matr)).^0.5;

%%%%%%%%%%%%
% PRINTING %
%%%%%%%%%%%%

fprintf('True parameters\n');
fprintf('Alpha:%8.4f',b_true(1));
fprintf('  Beta:%8.4f',b_true(2));
fprintf('  Sigma:%8.4f\n',b_true(3));
fprintf('\n',b_true(3));

fprintf('OLS Estimation\n');
fprintf('Estimated parameters\n');
fprintf('Alpha:%8.4f',b_hat(1));
fprintf('  Beta:%8.4f',b_hat(2));
%fprintf('  Sigma:%8.4f\n',b_hat(3));
%fprintf('Standard errors\n');
%fprintf('Alpha:%8.4f',stderr(1,1));
%fprintf('  Beta:%8.4f\n',stderr(2,2));

%%%%%%%%%%%%%%
% SORTED PAIRWISE ESTIMATION (WITHOUT CONNECTING FIRST AND LAST)  %
%%%%%%%%%%%%%%
pairwise_betas_1=zeros(2,T-1);
% iterate over all pairs
for i=(1:1:T-1)
    % x_matr is the X in betahat = (X'X)^(-1)(X'y) expression
    x_matr    = [ones(2,1), x(i:i+1)];
    % calculate betahat
    b_hat     = inv(x_matr'*x_matr)*x_matr'*y(i:i+1);
    pairwise_betas_1(:,i)=b_hat;
    % calculate (Tx1) vector of residuals
    %resid     = y-x_matr*b_hat;
    % estimate sigma: (root of) sum of squared residuals divided by degrees of freedom
    %sigma_hat = sqrt(sum(resid'*resid)/(T-2));
    %b_hat     = [b_hat; sigma_hat];
    % calculate estimated standard errors: this is (sigma^2*(X'X)^(-1))^0.5
    %stderr    = (sigma_hat*sigma_hat*inv(x_matr'*x_matr)).^0.5;
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
        % x_matr is the X in betahat = (X'X)^(-1)(X'y) expression
        x_matr    = [ones(2,1), x(i:i+1)];
        % calculate betahat
        b_hat     = inv(x_matr'*x_matr)*x_matr'*y(i:i+1);
        pairwise_betas_2(:,i)=b_hat;
    else
        x_first_and_last = [x(1); x(T)];
        y_first_and_last = [y(1); y(T)];
        x_matr    = [ones(2,1), x_first_and_last];
        b_hat     = inv(x_matr'*x_matr)*x_matr'*y_first_and_last;
        pairwise_betas_2(:,i)=b_hat;
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
