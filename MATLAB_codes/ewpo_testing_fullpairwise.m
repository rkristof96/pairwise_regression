% close all;
% clear all;
% clc;

% if you run the parallel tool you can set the number of core you use
% not necessary, bu you can speed up the simulations
%clust = parcluster;
%clust.NumWorkers=4;

% define date_when the run starts
date_var = datestr(datetime('today'));


% parameter sets
 T_list = [10,50,100,500,1000];%[10,50,100,500,1000];
 % correlation between x and u
 sigma_list = [0,0.2,0.5,0.8];
 % correlations between residuals used to generate the data
 if strcmp(version,'C')
    cov_list = [0,0.0391,0.2493,0.6397];
 else
    cov_list = [0,0,0,0];
 end
 
 reps = 1000;


%corr_matrix = zeros(length(sigma_list)*length(T_list)*2,2);
%beta_matrix = zeros(length(sigma_list)*length(T_list)*2,2);

%abs_dummy = 'True'; % or Abs
%dist = 'norm';
%sorted = 'False';
 
for sample_size_ind = 1:length(T_list)
%avg_bias_matrix = zeros(length(sigma_list)*2,2);
corr_matrix = zeros(length(sigma_list)*2,2);
beta_matrix = zeros(length(sigma_list)*2,4);
diff_alpha = zeros(length(sigma_list), reps);
diff_alpha2 = zeros(length(sigma_list), reps);
diff_beta = zeros(length(sigma_list), reps);
diff_beta2 = zeros(length(sigma_list), reps);
diff_matrix = zeros(length(sigma_list)*2,8);


T = T_list(sample_size_ind);

for sigma_ind = 1:length(sigma_list)

% true parameters

alpha = 1;
beta  = 0.5;
sigma = sigma_list(sigma_ind); % assign correlation between x and u
%rand('seed',202101);
rand('seed',20230318);
% generate x: (Tx1) u: (Txreps) vector of bi-variate normal distributed
% random variables

mu = zeros(1,reps+1);
cov_matrix = zeros(reps+1);
% assign correlation values between residuals
% necessary for the data generation
cov_matrix(:,:) = cov_list(sigma_ind); 
% first row and column from the covariance matrix includes the correlation
cov_matrix(1,:) = sigma;
cov_matrix(:,1) = sigma;

% diagonal of the covariance matrix includes ones
for diag_i=(1:1:reps+1)
    cov_matrix(diag_i,diag_i) = 1;
end

% select case uniform or normal DGP
if strcmp(dist,'unif')
%uniform case
    rng('default')  % For reproducibility
    R = mvnrnd(mu,cov_matrix,T); % generate data
    x = R(:,1);
    x = normcdf(x);
    x = unifinv(x,-3,2); %transform data to uniformly distributed data
end
if strcmp(dist,'norm')
    % normal case
    cov_matrix(1,1) = 1;
    rng('default')  % For reproducibility
    R = mvnrnd(mu,cov_matrix,T);
    x = R(:,1); % select x
end

if strcmp(dist,'student-t')
%uniform case
    rng('default')  % For reproducibility
    R = mvnrnd(mu,cov_matrix,T); % generate data
    x = R(:,1);
    x = normcdf(x);
    x = tinv(x,5); %transform data to t-distribution data
end

% in the first column are the x values all other are the error terms for
% the repetitions (u values)
eps = R(:,2:end);

if strcmp(version,'A')
	u = sigma*eps.^2.*(x-(-0.5));
end
if strcmp(version,'B')
	u = sigma*(x-(-0.5))+eps;
    %u = sigma*(x-(-0.5));
end

if strcmp(version,'C')
    % eps is not necessarily uncorrelated with x, if sigma is non-zero
    % v is uncorrelated with x
    v = normrnd(0,1,[T, reps]);
	u = eps + v;
end
% regression equation
y = alpha + beta*x + u;

% combine dependent and independent variables to a dataset
xy = [x y];
% sort values
if strcmp(sorted,'True')
    xy = sortrows(xy,1);
end
% reaasign values
x = xy(:,1);
y = xy(:,2:reps+1);


%%%%%%%%%%%%%%%%%%
% OLS ESTIMATION %
%%%%%%%%%%%%%%%%%%

x_matr = [ones(T,1) x];  % this is matrix X in betahat = (X'X)^(-1)*(X'y)
var_true = (sigma^2)*inv(x_matr'*x_matr);  % true variance-covariance matrix

% coefficients 
b_hat_all_ols = zeros(2,reps);  % sFtore estimated betahats, r-th repetition in r-th column
% correlation values between x and u are stored here
x_u_corr_all = zeros(1,reps);

%repetetion starts
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
    
   % betas beta_hat = (x-x_avg)'*(y-y_avg)./sum((x-x_avg).^2);
   % alpha_hat = y_avg - beta_hat*x_avg;
    b_hat              = numerator/denominator;
    b_hat_all_ols(2,r)     = b_hat;
    b_hat_all_ols(1,r)     = y_avg_r - b_hat*x_avg;
    
    u_hat              = y(:,r) - b_hat_all_ols(1,r) - b_hat_all_ols(2,r)*x;
    x_u_corr_all(1,r)  = corr(x, u_hat);
    
    r = r + 1;
end

standard_dev1=std(b_hat_all_ols(1,:));
standard_dev2=std(b_hat_all_ols(2,:));
%ols_alpha = mean(b_hat_all_ols(1,:));
%ols_beta = mean(b_hat_all_ols(2,:));

mean_x_u_corr = mean(x_u_corr_all);
stdev_x_u_corr = std(x_u_corr_all);

%%%%%%%%%%%%
% PRINTING %
%%%%%%%%%%%%

% fprintf('\nTrue parameters\n');
% fprintf('Alpha:%8.4f',b_true(1));
% fprintf('  Beta:%8.4f',b_true(2));
% fprintf('  Sigma:%8.4f\n',b_true(3));
% fprintf('Se(a):%8.4f',var_true(1,1)^0.5);
% fprintf('  Se(b):%7.4f\n',var_true(2,2)^0.5);
% print your results: means across Monte-Carlo repetitions
% fprintf('\n');
% fprintf('OLS Estimation\n');
% fprintf('Estimated parameters (mean of Monte Carlo repetitions)\n');
% fprintf('Alpha:%8.4f',mean(b_hat_all_ols(1,:),2));
% fprintf('  Beta:%8.4f\n',mean(b_hat_all_ols(2,:),2));
% fprintf('Standard errors (standard deviation of estimates at Monte Carlo repetitions)\n');
% fprintf('Alpha:%8.4f',standard_dev1);
% fprintf('  Beta:%8.4f\n',standard_dev2);
% fprintf('Mean correlation of x and u hat.\n');
% fprintf('Correlation:%8.4f\n',mean_x_u_corr);
% fprintf('Standard error of correlation of x and u hat.\n');
% fprintf('Standard error:%8.4f',stdev_x_u_corr);


%corr_matrix((sigma_ind-1)*2+1,1) = mean_x_u_corr;
%corr_matrix(sigma_ind*2,1) = stdev_x_u_corr;

%ols_avg_bias_alpha = mean(b_hat_all_ols(1,:)-alpha);
ols_avg_bias_beta = mean(b_hat_all_ols(2,:)-beta);

%ols_std_bias_alpha = std(b_hat_all_ols(1,:)-alpha);
ols_std_bias_beta = std(b_hat_all_ols(2,:)-beta);

% bias
%avg_bias_matrix((sigma_ind-1)*2+1,3) = ols_avg_bias_alpha;
%avg_bias_matrix(sigma_ind*2,3) = ols_std_bias_alpha;

avg_bias_matrix((sigma_ind-1)*2+1,2) = ols_avg_bias_beta;
avg_bias_matrix(sigma_ind*2,2) = ols_std_bias_beta;

% store coefficients
beta_matrix((sigma_ind-1)*2+1,3) = mean(b_hat_all_ols(1,:),2);
beta_matrix(sigma_ind*2,3) = standard_dev1;

beta_matrix((sigma_ind-1)*2+1,4) = mean(b_hat_all_ols(2,:),2);
beta_matrix(sigma_ind*2,4) = standard_dev2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FULL PAIRWISE ESTIMATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b_hat_all = zeros(2,reps);  % store estimated betahats, r-th repetition in r-th column
x_u_corr_all = zeros(1,reps);

r = 1;
while r < reps+0.5
    number_of_betas = T * (T-1) /2;
    pairwise_betas = zeros(2,number_of_betas);
    delta_y = zeros(1,number_of_betas);
    delta_x = zeros(1,number_of_betas);
    counter=1;

    % iterate over all pairs
    for i=(1:1:T-1)
        for j=((i+1):1:T)
           % if i<j
                % calculate x difference

                if strcmp(abs_dummy, 'True')
                    delta_x(1, counter) = abs(x(j,1) - x(i,1));
                else
                    delta_x(1, counter) = x(j,1) - x(i,1);
                end
                % calculate betahat
                x_avg     = (x(i,1)+x(j,1))/2;
                y_avg     = (y(i,r)+y(j,r))/2;
                numerator = y(j,r) - y(i,r);
                delta_y(1, counter) = numerator;
                denominator = x(j,1) - x(i,1);
                b_hat_i     = numerator/denominator;
                alpha_hat_i = y_avg - b_hat_i*x_avg;
                pairwise_betas(1,counter)=alpha_hat_i;
                pairwise_betas(2,counter)=b_hat_i;
                counter   = counter + 1;
          %  end
        end
    end

    % Obtain the delta-x weighted average of pairwise betas
            
    sum_delta_x = nansum(delta_x);
    inf_index =(pairwise_betas==Inf);
    minusinf_index = (pairwise_betas==-Inf);
    sum_noninf = max(max(inf_index(1,:),inf_index(2,:)),max(minusinf_index(1,:),minusinf_index(2,:)));
    %sum((inf_index+minusinf_index),1);
    
    noninf_index = (1-sum_noninf)==1; 
    
    weighted_parwise_betas = pairwise_betas(:,noninf_index)*delta_x( noninf_index)';

    weighted_average_parwise_betas = weighted_parwise_betas./sum_delta_x;
    
    b_hat_all(1,r)        = weighted_average_parwise_betas(1);
    b_hat_all(2,r)        = weighted_average_parwise_betas(2);
    
    u_hat              = y(:,r) - b_hat_all(1,r) - b_hat_all(2,r)*x;
    x_u_corr_all(1,r)  = corr(x, u_hat);
    
    % test stat
    delta_u = delta_y - b_hat_all(2,r)*delta_x;
    %test_stat(sigma_ind, r) = (sum(delta_x.*(delta_x*(beta-b_hat_all(2,r))+delta_u))/(T_list(sample_size_ind)*(T_list(sample_size_ind)-1))*2);
    
    r = r + 1;
end

standard_dev1=std(b_hat_all(1,:));
standard_dev2=std(b_hat_all(2,:));

%avg_bias_alpha = mean(b_hat_all(1,:)-alpha,2);
avg_bias_beta = mean(b_hat_all(2,:)-beta,2);

%std_bias_alpha = std(b_hat_all(1,:)-alpha);
std_bias_beta = std(b_hat_all(2,:)-beta);

mean_x_u_corr = mean(x_u_corr_all);
stdev_x_u_corr = std(x_u_corr_all);

%%%%%%%%%%%%
% PRINTING %
%%%%%%%%%%%%

% fprintf('\n');
% fprintf('\n FULL PAIRWISE ESTIMATION\n');
% fprintf('Estimated parameters (mean of Monte Carlo repetitions)\n');
% fprintf('Alpha:%8.4f',mean(b_hat_all(1,:),2));
% fprintf('  Beta:%8.4f\n',mean(b_hat_all(2,:),2));
% fprintf('Standard errors (standard deviation of estimates at Monte Carlo repetitions)\n');
% fprintf('Alpha:%8.4f',standard_dev1);
% fprintf('  Beta:%8.4f\n',standard_dev2);
% fprintf('Mean correlation of x and u hat.\n');
% fprintf('Correlation:%8.4f\n',mean_x_u_corr);
% fprintf('Standard error of correlation of x and u hat.\n');
% fprintf('Standard error:%8.8f',stdev_x_u_corr);

% results
%corr_matrix((sigma_ind-1)*2+1,2) = mean_x_u_corr;
%corr_matrix(sigma_ind*2,2) = stdev_x_u_corr;

%bias
%col 1 - alpha, col 2 - beta
%avg_bias_matrix((sigma_ind-1)*2+1,1) = avg_bias_alpha;
%avg_bias_matrix(sigma_ind*2,1) = std_bias_alpha;

avg_bias_matrix((sigma_ind-1)*2+1,1) = avg_bias_beta;
avg_bias_matrix(sigma_ind*2,1) = std_bias_beta;

beta_matrix((sigma_ind-1)*2+1,1) = nanmean(b_hat_all(1,:),2);
beta_matrix((sigma_ind-1)*2+1,2) = nanmean(b_hat_all(2,:),2);

nan_err =sum(isnan(b_hat_all(2,:)));
if nan_err ~= 0
   disp(num2str(nan_err))
end
beta_matrix(sigma_ind*2,1) = standard_dev1;
beta_matrix(sigma_ind*2,2) = standard_dev2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% calculate test statistics %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diff_alpha(sigma_ind,:) = (b_hat_all_ols(1,:)-b_hat_all(1,:));
diff_alpha2(sigma_ind,:) = (b_hat_all_ols(1,:)-b_hat_all(1,:)).^2;
diff_beta(sigma_ind,:) = (b_hat_all_ols(2,:)-b_hat_all(2,:));
diff_beta2(sigma_ind,:) = (b_hat_all_ols(2,:)-b_hat_all(2,:)).^2;
% mean
diff_matrix(sigma_ind*2-1,1) = mean((b_hat_all_ols(1,:)-b_hat_all(1,:)));
diff_matrix(sigma_ind*2-1,2) = mean((b_hat_all_ols(1,:)-b_hat_all(1,:)).^2);
diff_matrix(sigma_ind*2-1,3) = mean((b_hat_all_ols(2,:)-b_hat_all(2,:)));
diff_matrix(sigma_ind*2-1,4) = mean((b_hat_all_ols(2,:)-b_hat_all(2,:)).^2);
%std
diff_matrix(sigma_ind*2,1) = std((b_hat_all_ols(1,:)-b_hat_all(1,:)));
diff_matrix(sigma_ind*2,2) = std((b_hat_all_ols(1,:)-b_hat_all(1,:)).^2);
diff_matrix(sigma_ind*2,3) = std((b_hat_all_ols(2,:)-b_hat_all(2,:)));
diff_matrix(sigma_ind*2,4) = std((b_hat_all_ols(2,:)-b_hat_all(2,:)).^2);

% sum
diff_matrix(sigma_ind,5) = sum((b_hat_all_ols(1,:)-b_hat_all(1,:)));
diff_matrix(sigma_ind,6) = sum((b_hat_all_ols(1,:)-b_hat_all(1,:)).^2);
diff_matrix(sigma_ind,7) = sum((b_hat_all_ols(2,:)-b_hat_all(2,:)));
diff_matrix(sigma_ind,8) = sum((b_hat_all_ols(2,:)-b_hat_all(2,:)).^2);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(sorted,'True') 
    sorted_text = '_sorted';
else
    sorted_text = '';
end

if strcmp(abs_dummy,'True')==1
    if strcmp(dist,'norm')
        beta_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Norm_beta_version_',version,'_Abs_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        %corr_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Norm_corr_Abs_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        %tstat_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Norm_tstat_Abs_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        diff_alpha_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Norm_diff_alpha_version_',version,'_Abs_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        diff_alpha2_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Norm_diff_alpha2_version_',version,'_Abs_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        diff_beta_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Norm_diff_beta_version_',version,'_Abs_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        diff_beta2_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Norm_diff_beta2_version_',version,'_Abs_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        diff_matrix_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Norm_diff_matrix_version_',version,'_Abs_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        %bias_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Norm_bias_Abs_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
    else
        beta_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Unif_beta_version_',version,'_Abs_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        %corr_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Unif_corr_Abs_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        %tstat_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Unif_tstat_Abs_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        diff_alpha_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Unif_diff_alpha_version_',version,'_Abs_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        diff_alpha2_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Unif_diff_alpha2_version_',version,'_Abs_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        diff_beta_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Unif_diff_beta_version_',version,'_Abs_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        diff_beta2_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Unif_diff_beta2_version_',version,'_Abs_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        diff_matrix_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Unif_diff_matrix_version_',version,'_Abs_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat'); 
        %bias_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Unif_bias_Abs_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');  
    end 
else
    if strcmp(dist,'norm')
        beta_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Norm_beta_version_',version,'_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        %corr_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Norm_corr_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        %tstat_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Norm_tstat_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        diff_alpha_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Norm_diff_alpha_version_',version,'_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        diff_alpha2_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Norm_diff_alpha2_version_',version,'_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        diff_beta_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Norm_diff_beta_version_',version,'_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        diff_beta2_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Norm_diff_beta2_version_',version,'_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        diff_matrix_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Norm_diff_matrix_version_',version,'_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        %bias_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Norm_bias_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
    else
        beta_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Unif_beta_version_',version,'_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        %corr_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Unif_corr_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        %tstat_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Unif_tstat_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        diff_alpha_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Unif_diff_alpha_version_',version,'_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        diff_alpha2_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Unif_diff_alpha2_version_',version,'_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        diff_beta_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Unif_diff_beta_version_',version,'_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        diff_beta2_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Unif_diff_beta2_version_',version,'_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        diff_matrix_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Unif_diff_matrix_version_',version,'_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat'); 
        %bias_title = strcat('Results/ewpo_testing_fullpairwise',sorted_text,'_Unif_bias_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');  
    end
end
%Results/
parsave(beta_title, beta_matrix)
parsave(diff_alpha_title, diff_alpha)
parsave(diff_alpha2_title, diff_alpha2)
parsave(diff_beta_title, diff_beta)
parsave(diff_beta2_title, diff_beta2)
parsave(diff_matrix_title, diff_matrix)
%parsave(corr_title, corr_matrix)
%parsave(tstat_title, test_stat)
%parsave(bias_title, avg_bias_matrix)

end




