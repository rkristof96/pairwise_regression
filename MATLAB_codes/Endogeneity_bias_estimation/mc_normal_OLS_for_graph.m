close all;
clear all;
clc;

% true parameters

alpha = 0;
beta  = 4.2;
x_sigma = sqrt(100);
x_uniform_limit = 100;
v_sigma = sqrt(10);
%a = 0;
%b = 1;
integers_to_keep = [1];
%d = 1;
epsilon = 0.1;
b_true = [alpha;beta;x_sigma];

T = 500; % number of observations
reps = 1000; % number of Monte Carlo repetitions

for theta=(0:0.2:20)
    % explanatory variable
    rng('default');
    x = unifrnd(-x_uniform_limit,x_uniform_limit, [T,1]);
    %x = normrnd(0,x_sigma, [T,1]);

    % error terms

    randn('seed',202101);
    eps = normrnd(0,v_sigma, [T,reps]);  %generate (T x reps) matrix of normally distributed i.i.d. errors,
        %with mean 0 and variance sigma^2

    eps_endog = theta * x + eps;

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
        for j=1:length(integers_to_keep)
            if (integers_to_keep(j)-epsilon <= x(i,1)) && (x(i,1) <= integers_to_keep(j)+epsilon)
                counter = counter + 1;
            end
        end
    end

    x_kept = zeros(counter, 1);
    y_kept = zeros(counter, reps);
    eps_kept = zeros(counter, reps);

    holder = 1;

    for i=(1:1:T)
        for j=1:length(integers_to_keep)
            if (integers_to_keep(j)-epsilon <= x(i,1)) && (x(i,1) <= integers_to_keep(j)+epsilon)
                x_kept(holder, 1) = integers_to_keep(j);
                y_kept(holder, :) = y(i, :);
                eps_kept(holder, :) = eps(i,:);
                holder = holder +1;
            end
        end
    end

    %x = x_kept;
    %y = y_kept;
    %eps = eps_kept;

    b_hat_all = zeros(1,reps);  % store estimated betahats, r-th repetition in r-th column
    bias_hat_all = zeros(1,reps);
    bias_hat_all_perc = zeros(1,reps);

    x_eps_corr = zeros(1,reps);

    x_sigma = sqrt(var(x));

    r = 1;
    while r < reps+0.5
        numerator = x(:, 1)' * y(:,r);
        denominator = x(:, 1)' * x(:, 1);

        b_hat              = numerator/denominator;
        %bias_hat           = b_hat-beta;

        corr_matrix       = corrcoef(x,eps(:,r));
        x_eps_corr(1,r)   = corr_matrix(2,1);

        bias_hat           = (corr_matrix(2,1) / sqrt(1-corr_matrix(2,1)^2)) * v_sigma / x_sigma;

        b_hat_all(1,r)     = b_hat;
        bias_hat_all(1,r)  = bias_hat;
        
        bias_hat_all_perc(1,r) = bias_hat/b_hat*100;

        r = r + 1;
    end

    standard_dev1=std(b_hat_all(1,:));

    standard_dev2=std(bias_hat_all(1,:));
    
    standard_dev3=std(bias_hat_all_perc(1,:));

    x_eps_corr_mean = mean(x_eps_corr,2);

    standard_dev4 = std(x_eps_corr);

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
    fprintf('Estimated parameters (mean of Monte Carlo repetitions)\n');
    fprintf('  Bias:%8.4f\n',mean(bias_hat_all(1,:),2));
    fprintf('Standard errors (standard deviation of estimates at Monte Carlo repetitions)\n');
    fprintf('  Bias:%8.4f\n',standard_dev2);
    fprintf('Estimated parameters (mean of Monte Carlo repetitions)\n');
    fprintf('  Bias percentage:%8.4f\n',mean(bias_hat_all_perc(1,:),2));
    fprintf('Standard errors (standard deviation of estimates at Monte Carlo repetitions)\n');
    fprintf('  Bias percentage:%8.4f\n',standard_dev3);
    %fprintf('Number of observations we keep\n');
    %fprintf('  N:%8.4f\n',counter);
    fprintf('  Correlation coefficient:%8.4f\n',x_eps_corr_mean);
    fprintf('  St dev of correlation coefficient:%8.4f\n',standard_dev4);
    fprintf('  Estimated variance of x:%8.4f',x_sigma^2);

    to_save = [theta, mean(bias_hat_all_perc(1,:),2), standard_dev3, x_eps_corr_mean];

    fid = fopen('hello.txt', 'a+');
    fprintf(fid, '%6.4f %6.4f %6.4f %6.4f\n', to_save);
    fclose(fid);
end