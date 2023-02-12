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
extreme_value = 0.9792;

T = 500; % number of observations
reps = 1000; % number of Monte Carlo repetitions

for theta=(0:0.1:4.8)
    %%%%%%%%%%%%%%
    % OLS ESTIMATION %
    %%%%%%%%%%%%%%
    
    % explanatory variable
    rng('default');
    x = normrnd(0,x_sigma, [T,1]);
    
    % error terms

    randn('seed',202101);
    v = normrnd(0,v_sigma, [T,reps]);  %generate (T x reps) matrix of normally distributed i.i.d. errors,
        %with mean 0 and variance sigma^2

    u = theta * x + v;

    % dependent variables, in each of the repetitions

    y = alpha+beta*x+u;  % (T x reps) matrix of dependent variables
    
    b_hat_all_ols = zeros(1,reps);  % store estimated betahats, r-th repetition in r-th column
    bias_corrected_b_hat_lower_ols = zeros(1,reps);
    bias_corrected_b_hat_upper_ols = zeros(1,reps);

    x_u_corr_ols = zeros(1,reps);

    x_sigma_estimated = sqrt(var(x));

    r = 1;
    while r < reps+0.5
        numerator = x(:, 1)' * y(:,r);
        denominator = x(:, 1)' * x(:, 1);

        b_hat_ols              = numerator/denominator;

        corr_matrix       = corrcoef(x,u(:,r));
        x_u_corr_ols(1,r)   = corr_matrix(2,1);
        
        estimated_error = y(:,r) - b_hat_ols * x(:, 1);
        
        v_sigma_estimated = std(estimated_error);

        bias_hat_upper_ols           = (extreme_value / sqrt(1-extreme_value^2)) ...
                              * v_sigma_estimated / x_sigma_estimated;

        bias_hat_lower_ols           = (-extreme_value / sqrt(1-extreme_value^2)) ...
                              * v_sigma_estimated / x_sigma_estimated;
                          
        bias_corrected_beta_hat_lower = b_hat_ols - bias_hat_upper_ols;
        bias_corrected_beta_hat_upper = b_hat_ols - bias_hat_lower_ols;
                              
        b_hat_all_ols(1,r)     = b_hat_ols;
        bias_corrected_b_hat_lower_ols(1,r) = bias_corrected_beta_hat_lower;
        bias_corrected_b_hat_upper_ols(1,r) = bias_corrected_beta_hat_upper;
        
        r = r + 1;
    end

    b_hat_all_ols_mean = mean(b_hat_all_ols,2);
    
    standard_dev1_ols=std(b_hat_all_ols(1,:));
        
    bias_corrected_b_hat_lower_ols_mean = mean(bias_corrected_b_hat_lower_ols,2);
    
    standard_dev2_ols = std(bias_corrected_b_hat_lower_ols(1,:));
    
    bias_corrected_b_hat_upper_ols_mean = mean(bias_corrected_b_hat_upper_ols,2);
    
    standard_dev3_ols = std(bias_corrected_b_hat_upper_ols(1,:));
    
    x_u_corr_ols_mean = mean(x_u_corr_ols,2);
    
    %%%%%%%%%%%%%%
    % IV ESTIMATION %
    %%%%%%%%%%%%%%
    
    b_hat_all_iv = zeros(1,reps);  % store estimated betahats, r-th repetition in r-th column
    bias_hat_all_iv = zeros(1,reps);

    x_u_corr_iv = zeros(1,reps);

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

        b_hat_iv              = numerator/denominator;
        
        bias_hat_iv           = b_hat_iv-beta;
        
        b_hat_all_iv(1,r)     = b_hat_iv;
        bias_hat_all_iv(1,r)  = bias_hat_iv;
                
        corr_matrix       = corrcoef(x,u);
        x_u_corr_iv(1,r)   = corr_matrix(2,1);
        
        r = r + 1;
    end

    b_hat_all_iv_mean = mean(b_hat_all_iv,2);
    
    standard_dev1_iv=std(b_hat_all_iv(1,:));
    
    x_u_corr_iv_mean = mean(x_u_corr_iv,2);
    
    %%%%%%%%%%%%
    % PRINTING %
    %%%%%%%%%%%%

    to_save = [theta, b_hat_all_ols_mean, standard_dev1_ols, bias_corrected_b_hat_lower_ols_mean, ...
               standard_dev2_ols, bias_corrected_b_hat_upper_ols_mean, standard_dev3_ols, ...
               x_u_corr_ols_mean, b_hat_all_iv_mean, standard_dev1_iv, x_u_corr_iv_mean];

    fid = fopen('hello.txt', 'a+');
    fprintf(fid, '%6.4f %6.4f %6.10f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n', to_save);
    fclose(fid);
    
    fprintf('Theta:%8.4f\n',theta)
end