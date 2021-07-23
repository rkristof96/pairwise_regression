close all;
clear all;
clc;

% true parameters

alpha = 1;
beta  = 0;
sigma = 1;

b_true = [alpha;beta;sigma];

T = 50; % number of observations
reps = 1000; % number of Monte Carlo repetitions


%%%%%%%%%%%%%%%%%%%
% DATA GENERATION %
%%%%%%%%%%%%%%%%%%%

% explanatory variable
rand('seed',202101);
% generate x: (Tx1) vector of uniformly distributed random
%    variables  
x = rand(T,1)*500;
x_sorted = sort(x);
x = x_sorted;

% error terms
% generate i.i.d. error terms, each normally distributed with 0
%    expected value and sigma^2 variance
randn('seed',202102);
eps = normrnd(0,sigma, [T,reps]);

% generate the dependent variable
y = alpha+beta*x+eps;

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
    beta_hat              = numerator/denominator;
    b_hat_all(2,r)     = beta_hat;
    b_hat_all(1,r)     = y_avg_r - beta_hat*x_avg;

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
% AVERAGING Y %
%%%%%%%%%%%%%%

x_floored = unique(floor(x));
T = length(x_floored);
y_averaged = zeros(T,reps);

for r=(1:reps)
    place_holder=1;
    for i=(1:T)
        if i==1
            current_sum=y(i,r);
            current_counter=1;
        elseif x(i,1)~=x(i-1,1)
            current_y=current_sum/current_counter;
            y_averaged(place_holder,r)=current_y;
            current_sum=y(i,r);
            current_counter=1;        
            place_holder=place_holder+1;
            if i==T
                y_averaged(place_holder,r)=current_sum;
            end
        else
            current_sum=current_sum+y(i,r);
            current_counter=current_counter+1;
            if i==T
                current_y=current_sum/current_counter;
                y_averaged(place_holder,r)=current_y;
            end            
        end
    end
end

x = x_floored;
y = y_averaged;
T = length(x);


%%%%%%%%%%%%%%
% AVERAGED FLOORED SORTED PAIRWISE ESTIMATION (WITH CONNECTING FIRST AND LAST) %
%%%%%%%%%%%%%%
b_hat_all = zeros(2,reps);  % store estimated betahats, r-th repetition in r-th column

r = 1;
while r < reps+0.5
    pairwise_betas=zeros(2,T);
    for i=(1:1:T)
        if i~=T
            % calculate betahat
            x_avg     = mean(x(i:i+1));
            y_avg     = mean(y(i:i+1, r));
            numerator = y(i+1,r) - y(i,r);
            denominator = x(i+1,1) - x(i,1);
            b_hat_i     = numerator/denominator;
            alpha_hat_i = mean(y(i:i+1,r)) - b_hat_i*mean(x(i:i+1));
            pairwise_betas(1,i)=alpha_hat_i;
            pairwise_betas(2,i)=b_hat_i;
        else
            x_first_and_last = [x(1); x(T)];
            y_first_and_last = [y(1,r); y(T,r)];
            x_avg     = mean(x_first_and_last);
            y_avg     = mean(y_first_and_last);
            numerator = y(i,r) - y(1,r);
            denominator = x(i,1) - x(1,1);
            b_hat_i     = numerator/denominator;
            alpha_hat_i = mean(y_first_and_last) - b_hat_i*mean(x_first_and_last);
            pairwise_betas(1,i)=alpha_hat_i;
            pairwise_betas(2,i)=b_hat_i;
        end;
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
fprintf('\n SORTED PAIRWISE ESTIMATION (WITH CONNECTING FIRST AND LAST)\n');
fprintf('Estimated parameters (mean of Monte Carlo repetitions)\n');
fprintf('Alpha:%8.4f',mean(b_hat_all(1,:),2));
fprintf('  Beta:%8.4f\n',mean(b_hat_all(2,:),2));
fprintf('Standard errors (standard deviation of estimates at Monte Carlo repetitions)\n');
fprintf('Alpha:%8.4f',standard_dev1);
fprintf('  Beta:%8.4f',standard_dev2);
