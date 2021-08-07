close all;
clear all;
clc;

% true parameters

alpha = 1;
beta  = 1.5;
sigma = 1;
epsilon = 0.01;
d = 0.5;

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
x = rand(T,1)*2-1;
%x = normrnd(50,25, [T,1]);

% error terms

randn('seed',202101);
eps = normrnd(0,sigma, [T,reps]);  %generate (T x reps) matrix of normally distributed i.i.d. errors,
    %with mean 0 and variance sigma^2

% dependent variables, in each of the repetitions

y = alpha+beta*x+eps;  % (T x reps) matrix of dependent variables

% sort
xy = [x y];

xy = sortrows(xy,1);

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
% Full PAIRWISE ESTIMATION %
%%%%%%%%%%%%%%

b_hat_all = zeros(1,reps);  % store estimated betahats, r-th repetition in r-th column

% Keep only observations that meet our criteria
no_to_keep = 0;
for i=(1:1:T)
    for j=(0:1:20)
        absolute_deviation = abs(abs(x(i, 1))-j*d);
        if absolute_deviation<epsilon
            no_to_keep = no_to_keep + 1;
        end
    end         
end

x_kept = zeros(no_to_keep,1);
y_kept = zeros(no_to_keep,reps);

counter = 1;

for i=(1:1:T)
    for j=(0:1:20)
        absolute_deviation = abs(abs(x(i, 1))-j*d);
        if absolute_deviation<epsilon
            x_kept(counter,1) = x(i,1);
            y_kept(counter,:) = y(i,:);
            counter = counter +1;
        end
    end         
end

x = x_kept;
y = y_kept;
T = no_to_keep;

r = 1;
while r < reps+0.5
    number_of_pairs = T * (T-1) /2;
    delta_x = zeros(1,number_of_pairs);
    delta_y = zeros(1,number_of_pairs);
    counter = 1;
    
    x1_aux = zeros(1,number_of_pairs);
    x2_aux = zeros(1,number_of_pairs);
    y1_aux = zeros(1,number_of_pairs);
    y2_aux = zeros(1,number_of_pairs);    
    
    % iterate over all pairs
    for i=(1:1:T-1)
        for j=(2:1:T)
            if i<j
                % calculate x difference
                delta_x(1, counter) = x(j,1) - x(i,1);
                delta_y(1, counter) = y(j,r) - y(i,r);
                x1_aux(1, counter) = x(j,1);
                x2_aux(1, counter) = x(i,1);
                y1_aux(1, counter) = y(j,r);
                y2_aux(1, counter) = y(i,r);
                counter = counter + 1;
            end
        end
    end
    
    sum_delta_y = 0;
    N = 0;
    pairwise_beta_1_sum = 0;
    pairwise_beta_0_sum = 0;
    for i=(1:1:number_of_pairs)
        absolute_deviation = abs(abs(delta_x(1, i))-d);
        if absolute_deviation<2*epsilon
            sum_delta_y = sum_delta_y + delta_y(1, i);
            x_avg     = mean([x1_aux(1, i) x2_aux(1, i)]);
            y_avg     = mean([y1_aux(1, i) y2_aux(1, i)]);
            beta_hat_i = delta_y(1,i)/delta_x(1,i);
            alpha_hat_i = y_avg - beta_hat_i*x_avg;
            pairwise_beta_1_sum = pairwise_beta_1_sum + beta_hat_i;
            pairwise_beta_0_sum =pairwise_beta_0_sum + alpha_hat_i;              
            N = N+1;
        end
    end
    
    %estimate beta
    %beta_hat = sum_delta_y/(N*d);
    beta_hat = pairwise_beta_1_sum/N;
    alpha_hat = pairwise_beta_0_sum/N;
    
    b_hat_all(1,r)        = alpha_hat;
    b_hat_all(2,r)        = beta_hat;
    
    r = r + 1;
end

standard_dev1=std(b_hat_all(1,:));
standard_dev2=std(b_hat_all(2,:));

%%%%%%%%%%%%
% PRINTING %
%%%%%%%%%%%%

fprintf('\n');
fprintf('\n FULL PAIRWISE ESTIMATION\n');
fprintf('Estimated parameters (mean of Monte Carlo repetitions)\n');
fprintf('Alpha:%8.4f',mean(b_hat_all(1,:),2));
fprintf('  Beta:%8.4f\n',mean(b_hat_all(2,:),2));
fprintf('Standard errors (standard deviation of estimates at Monte Carlo repetitions)\n');
fprintf('Alpha:%8.4f',standard_dev1);
fprintf('  Beta:%8.4f',standard_dev2);
fprintf('\n  Number of observations we keep:%8.4f',N);
fprintf('\n  d is equal to:%8.4f',d);
