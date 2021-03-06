close all;
clear all;
clc;

% true parameters

alpha = 1;
beta  = 0.5;
sigma = 1;

b_true = [alpha;beta;sigma];

T = 5000; % number of observations
reps = 1000; % number of Monte Carlo repetitions


%%%%%%%%%%%%%%%%%%%
% DATA GENERATION %
%%%%%%%%%%%%%%%%%%%

% explanatory variable
rand('seed',202101);
% generate x: (Tx1) vector of uniformly distributed random
%    variables on the interval (-1;+1) 
%x = rand(T,1)*2-1;
x = normrnd(50,25, [T,1]);


% error terms
% generate i.i.d. error terms, each normally distributed with 0
%    expected value and sigma^2 variance
randn('seed',202102);
eps = normrnd(0,sigma, [T,reps]);

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
% ADJACENT PAIRWISE ESTIMATION (WITHOUT CONNECTING FIRST AND LAST)  %
%%%%%%%%%%%%%%
b_hat_all = zeros(2,reps);  % store estimated betahats, r-th repetition in r-th column

r = 1;
while r < reps+0.5
    pairwise_betas=zeros(2,T);
    
    for i=(1:1:T-1)
        % calculate betahat
        x_avg     = mean(x(i:i+1));
        y_avg     = mean(y(i:i+1, r));
        numerator = y(i+1,r) - y(i,r);
        denominator = x(i+1,1) - x(i,1);
        b_hat     = numerator/denominator;
        alpha_hat = mean(y(i:i+1,r)) - b_hat*mean(x(i:i+1));
        pairwise_betas(1,i)=alpha_hat;
        pairwise_betas(2,i)=b_hat;
    end;
    
    average_parwise_betas = mean(pairwise_betas,2);
    b_hat_all(1,r)        = average_parwise_betas(1,:);
    b_hat_all(2,r)        = average_parwise_betas(2,:);

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


%%%%%%%%%%%%%%
% SORTED PAIRWISE ESTIMATION (WITH CONNECTING FIRST AND LAST)  %
%%%%%%%%%%%%%%
%pairwise_betas_2=zeros(2,T);
% iterate over all pairs
%for i=(1:1:T)
%    if i~=T
%        % calculate betahat
%        x_avg     = mean(x(i:i+1));
%        y_avg     = mean(y(i:i+1));
%        numerator = y(i+1,1) - y(i,1);
%        denominator = x(i+1,1) - x(i,1);
%        b_hat     = numerator/denominator;
%        alpha_hat = mean(y(i:i+1)) - b_hat*mean(x(i:i+1));
%        pairwise_betas_2(1,i)=alpha_hat;
%        pairwise_betas_2(2,i)=b_hat;
%    else
%        x_first_and_last = [x(1); x(T)];
%        y_first_and_last = [y(1); y(T)];
%        % calculate betahat
%        x_avg     = mean(x_first_and_last);
%        y_avg     = mean(y_first_and_last);
%        numerator = y(i,1) - y(1,1);
%        denominator = x(i,1) - x(1,1);
%        b_hat     = numerator/denominator;
%        alpha_hat = mean(y_first_and_last) - b_hat*mean(x_first_and_last);
%        pairwise_betas_2(1,i)=alpha_hat;
%        pairwise_betas_2(2,i)=b_hat;
%    end;
%    % calculate (Tx1) vector of residuals
%    %resid     = y-x_matr*b_hat;
%    % estimate sigma: (root of) sum of squared residuals divided by degrees of freedom
%    %sigma_hat = sqrt(sum(resid'*resid)/(T-2));
%    %b_hat     = [b_hat; sigma_hat];
%    % calculate estimated standard errors: this is (sigma^2*(X'X)^(-1))^0.5
%    %stderr    = (sigma_hat*sigma_hat*inv(x_matr'*x_matr)).^0.5;
%end;

%display average of betas
%average_parwise_betas_2=mean(pairwise_betas_2,2);
%fprintf('\n');
%fprintf('SORTED PAIRWISE ESTIMATION (WITH CONNECTING FIRST AND LAST)\n');
%fprintf('Alpha:%8.4f',average_parwise_betas_2(1,:));
%fprintf('  Beta:%8.4f\n',average_parwise_betas_2(2,:));
