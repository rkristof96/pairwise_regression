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
Psi_0 = 0.6;
b_true = [alpha;beta;x_sigma];
extreme_value = 0.99;

T = 500; % number of observations
reps = 1000; % number of Monte Carlo repetitions

max_theta = 1;
        
for theta=(0:0.1:max_theta)
    
    % Calculate rho_zero using formula
    
    rho_zero = theta * x_sigma / sqrt(theta^2 * x_sigma^2 + v_sigma^2);
    
    b_hat_all_ols = zeros(1,reps);  % store estimated betahats, r-th repetition in r-th column
    bias_as_percentage_all_ols = zeros(1,reps);  % store estimated betahats, r-th repetition in r-th column
    
    b_hat_all_iv = zeros(1,reps);  % store estimated betahats, r-th repetition in r-th column
    bias_as_percentage_all_iv = zeros(1,reps);
    
    r = 1;
    while r < reps+0.5
        
        %%%%%%%%%%%%%%
        % GENERATE DATA %
        %%%%%%%%%%%%%%
        
        % explanatory variable, instrumental variable and v ->
        % multivariate standard normal -> transfor x and v to uniform
           
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
        
        R_normcdf = normcdf(R);
        
        R_transformed = [norminv(R_normcdf(:,1)) R_normcdf(:,2) norminv(R_normcdf(:,3))];
        
        z = R_transformed(:,1);
        x = R_transformed(:,2);
        v = R_transformed(:,3);
        
        % endogenous u
        
        u = theta * x + v;
        
        % dependent variables

        y = alpha+beta*x+u; 
        
        %%%%%%%%%%%%%%
        % OLS ESTIMATION %
        %%%%%%%%%%%%%%
        numerator = x' * y;
        denominator = x' * x;

        b_hat_ols              = numerator/denominator;
        bias_as_percentage_ols = (b_hat_ols - beta)/beta;
        
        b_hat_all_ols(1,r)     = b_hat_ols;
        bias_as_percentage_all_ols(1,r) = bias_as_percentage_ols;
        
        r = r + 1;
    
        %%%%%%%%%%%%%%
        % IV ESTIMATION %
        %%%%%%%%%%%%%%
        
        numerator = z' * y;
        denominator = z' * x;

        b_hat_iv              = numerator/denominator;
        bias_as_percentage_iv = (b_hat_iv - beta)/beta;
        
        b_hat_all_iv(1,r)     = b_hat_iv;
        bias_as_percentage_all_iv(1,r) = bias_as_percentage_iv;
        
        r = r + 1;
    end
    
    b_hat_all_ols_mean = mean(b_hat_all_ols,2);
    
    standard_dev1_ols=std(b_hat_all_ols(1,:));
        
    bias_as_percentage_all_ols_mean = mean(bias_as_percentage_all_ols,2);
    
    standard_dev2_ols=std(bias_as_percentage_all_ols(1,:));

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
end

figure
scatterhist(z,x,'Direction','out')