close all;
clear all;
clc;

% true parameters

alpha = 0;
beta  = 4.2;
x_sigma = sqrt(1);
z_sigma = sqrt(1);
v_sigma = sqrt(1);
epsilon_0 = 0.8;
b_true = [alpha;beta;x_sigma];

T = 500; % number of observations
reps = 1000; % number of Monte Carlo repetitions

for theta=(0:0.1:0.7)
    %%%%%%%%%%%%%%
    % IV ESTIMATION %
    %%%%%%%%%%%%%%
    
    b_hat_all = zeros(1,reps);  % store estimated betahats, r-th repetition in r-th column
    bias_hat_all = zeros(1,reps);
    bias_hat_all_perc = zeros(1,reps);

    x_u_corr = zeros(1,reps);

    r = 1;
    while r < reps+0.5
        % explanatory variable, instrumental variable and v ->
        % multivariate standard normal
           
        mu = zeros(3, 1);
        
        var_covar_matrix = zeros(3, 3);

        var_covar_matrix(1,1) = z_sigma^2;
        var_covar_matrix(2,2) = x_sigma^2;
        var_covar_matrix(3,3) = v_sigma^2;
        
        var_covar_matrix(1,2) = epsilon_0 * z_sigma * x_sigma;
        var_covar_matrix(2,1) = epsilon_0 * z_sigma * x_sigma;
        
        var_covar_matrix(1,3) = - theta * epsilon_0 * z_sigma * x_sigma;
        var_covar_matrix(3,1) = - theta * epsilon_0 * z_sigma * x_sigma;

        rng('default');
        randn('seed',r);
        R = mvnrnd(mu,var_covar_matrix,T);
        
        z = R(:,1);
        x = R(:,2);
        v = R(:,3);
        
        % endogenous u
        
        u = theta * x + v;
        
        % dependent variables, in each of the repetitions

        y = alpha+beta*x+u;  % (T x reps) matrix of dependent variables
        
        % IV formula
        
        numerator = z' * y;
        denominator = z' * x;

        b_hat              = numerator/denominator;
        
        bias_hat           = b_hat-beta;
        
        b_hat_all(1,r)     = b_hat;
        bias_hat_all(1,r)  = bias_hat;
        
        bias_hat_all_perc(1,r) = bias_hat/b_hat*100;
        
        corr_matrix       = corrcoef(x,u);
        x_u_corr(1,r)   = corr_matrix(2,1);
        
        r = r + 1;
    end

    standard_dev1=std(b_hat_all(1,:));

    standard_dev2=std(bias_hat_all(1,:));
    
    standard_dev3=std(bias_hat_all_perc(1,:));

    x_u_corr_mean = mean(x_u_corr,2);

    standard_dev4 = std(x_u_corr);
    
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
    fprintf('  Correlation coefficient:%8.4f\n',x_u_corr_mean);
    fprintf('  St dev of correlation coefficient:%8.4f\n',standard_dev4);
    
    to_save = [theta, mean(bias_hat_all(1,:),2), standard_dev2, x_u_corr_mean];

    fid = fopen('hello.txt', 'a+');
    fprintf(fid, '%6.4f %6.4f %6.4f %6.4f\n', to_save);
    fclose(fid);
end