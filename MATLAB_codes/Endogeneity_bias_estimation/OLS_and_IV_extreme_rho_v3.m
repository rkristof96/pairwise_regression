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
Psi_0 = 0.4;
b_true = [alpha;beta;x_sigma];
extreme_value = 0.99;

T = 500; % number of observations
reps = 1000; % number of Monte Carlo repetitions

max_theta = 1.7;

% explanatory variable
rng('default');
x = normrnd(0,x_sigma, [T,1]);
    
% error terms
randn('seed',202101);
v = normrnd(0,v_sigma, [T,reps]);  %generate (T x reps) matrix of normally distributed i.i.d. errors,
        %with mean 0 and variance sigma^2

counter = 1;
        
for theta=(0:0.1:max_theta)
    
    %%%%%%%%%%%%%%
    % OLS ESTIMATION %
    %%%%%%%%%%%%%%

    u = theta * x + v;

    % dependent variables, in each of the repetitions

    y = alpha+beta*x+u;  % (T x reps) matrix of dependent variables
    
    % Calculate rho_zero using formula
    
    rho_zero = theta * x_sigma / sqrt(theta^2 * x_sigma^2 + v_sigma^2);
    
    b_hat_all_ols = zeros(1,reps);  % store estimated betahats, r-th repetition in r-th column
    bias_as_percentage_all_ols = zeros(1,reps);  % store estimated betahats, r-th repetition in r-th column
    
    r = 1;
    while r < reps+0.5
        numerator = x(:, 1)' * y(:,r);
        denominator = x(:, 1)' * x(:, 1);

        b_hat_ols              = numerator/denominator;

        bias_as_percentage_ols = (b_hat_ols - beta)/beta;
        
        b_hat_all_ols(1,r)     = b_hat_ols;
        bias_as_percentage_all_ols(1,r) = bias_as_percentage_ols;
        
        r = r + 1;
    end

    b_hat_all_ols_mean = mean(b_hat_all_ols,2);
    
    standard_dev1_ols=std(b_hat_all_ols(1,:));
        
    bias_as_percentage_all_ols_mean = mean(bias_as_percentage_all_ols,2);
    
    standard_dev2_ols=std(bias_as_percentage_all_ols(1,:));
    
    %%%%%%%%%%%%%%
    % IV ESTIMATION %
    %%%%%%%%%%%%%%
    
    b_hat_all_iv = zeros(1,reps);  % store estimated betahats, r-th repetition in r-th column
    bias_as_percentage_all_iv = zeros(1,reps);

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
        
        var_covar_matrix(1,3) = Psi_0 * sqrt(theta^2*x_sigma^2 + v_sigma^2) ...
                               * z_sigma - theta * epsilon_0 * z_sigma * x_sigma;
        var_covar_matrix(3,1) = Psi_0 * sqrt(theta^2*x_sigma^2 + v_sigma^2) ...
                               * z_sigma - theta * epsilon_0 * z_sigma * x_sigma;

        rng('default');
        randn('seed',r);
        R = mvnrnd(mu,var_covar_matrix,T);
        
        z_iv = R(:,1);
        x_iv = R(:,2);
        v_iv = R(:,3);
        
        % endogenous u
        
        u_iv = theta * x_iv + v_iv;
        
        % dependent variables

        y_iv = alpha+beta*x_iv+u_iv;  
        
        % IV formula
        
        numerator = z_iv' * y_iv;
        denominator = z_iv' * x_iv;

        b_hat_iv              = numerator/denominator;
        bias_as_percentage_iv = (b_hat_iv - beta)/beta;
        
        b_hat_all_iv(1,r)     = b_hat_iv;
        bias_as_percentage_all_iv(1,r) = bias_as_percentage_iv;
        
        r = r + 1;
    end

    b_hat_all_iv_mean = mean(b_hat_all_iv,2);
    
    standard_dev1_iv=std(b_hat_all_iv(1,:));
    
    bias_as_percentage_all_iv_mean = mean(bias_as_percentage_all_iv,2);
    
    standard_dev2_iv = std(bias_as_percentage_all_iv(1,:));
    
    %%%%%%%%%%%%
    % PRINTING %
    %%%%%%%%%%%%

    to_save = [theta, rho_zero, bias_as_percentage_all_ols_mean, standard_dev2_ols, ...
               bias_as_percentage_all_iv_mean, standard_dev2_iv];

    fid = fopen('hello.txt', 'a+');
    fprintf(fid, '%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n', to_save);
    fclose(fid);
    
    fprintf('Theta:%8.4f\n',theta)
    
    counter = counter + 1;
end