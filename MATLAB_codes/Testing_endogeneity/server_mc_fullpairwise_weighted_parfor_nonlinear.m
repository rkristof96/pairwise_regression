% close all;
% clear all;
% clc;

% if you run the parallel tool you can set the number of core you use
% not necessary, bu you can speed up the simulations
%clust = parcluster;
%clust.NumWorkers=4;

% define date_when the run starts
date_var = '20072023';%datestr(datetime('today'));


% parameter sets
 T_list = [500];%[50, 500, 1000];%[50];, 5000
 % correlation between x and u
 sigma_list = [0,0.2,0.5,0.8];
 % correlations between residuals used to generate the data
 cov_list = [0,0.0391,0.2493,0.6397];
 
 reps = 1000;


%corr_matrix = zeros(length(sigma_list)*length(T_list)*2,2);
%beta_matrix = zeros(length(sigma_list)*length(T_list)*2,2);

abs_dummy = 'False'; % or Abs
dist = 'uniform';
 
for sample_size_ind = 1:length(T_list)
corr_matrix = zeros(length(sigma_list)*2,2);
beta_matrix = zeros(length(sigma_list)*2,2);
test_stat = zeros(length(sigma_list), reps);

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

v = randn([T, reps]);

%mu = zeros(1,reps+1);
%cov_matrix = zeros(reps+1);
% assign correlation values between residuals
% necessary for the data generation
%cov_matrix(:,:) = cov_list(sigma_ind); 
% first row and column from the covariance matrix includes the correlation
% value between x and u called sigma
%cov_matrix(1,:) = sigma;
%cov_matrix(:,1) = sigma;

% diagonal of the covariance matrix includes ones
%for diag_i=(1:1:reps+1)
%    cov_matrix(diag_i,diag_i) = 1;
%end

% select case uniform or normal DGP
if strcmp(dist,'uniform')
%uniform case
    %rng('default')  % For reproducibility
    %R = mvnrnd(mu,cov_matrix,T); % generate data
    %x = R(:,1);
    %x = normcdf(x);
    %x = unifinv(x,0,50); %transform data to uniformly distributed data
    x = unifrnd(0,50,[T,1]);
    mu = mean(x);
else
    % normal case
    %cov_matrix(1,1) = 5;
    %rng('default')  % For reproducibility
    %R = mvnrnd(mu,cov_matrix,T);
    %x = exp(R(:,1)); % select x
end

u = zeros(T,reps);
y = zeros(T,reps);

for i =1:reps
u(:,i) = (v(:,i) + sigma).*(x-mu);
y(:,i) = alpha+beta*x+u(:,i);  
end
%u = (x-mu).*(v + sigma);


% in the first column are the x values all other are the error terms for
% the repetitions (u values)
%eps = R(:,2:end);

% dependent variables, in each of the repetitions
%y = alpha+beta*x+u;  % (T x reps) matrix of dependent variables

var_y = std(y);

omega = 1/T*sum(exp(y).*logistic(y).^2.*(1+2*exp(-y).*logistic(y)),1);



% combine dependent and independent variables to a dataset
%xy = [x y];
% sort values
%xy = sortrows(xy,1);
% reaasign values
%x = xy(:,1);
%y = xy(:,2:reps+1);


%%%%%%%%%%%%%%
% OLS ESTIMATION %
%%%%%%%%%%%%%%

x_matr = [ones(T,1) x];  % this is matrix X in betahat = (X'X)^(-1)*(X'y)
var_true = (sigma^2)*inv(x_matr'*x_matr);  % true variance-covariance matrix

% coefficients 
b_hat_all = zeros(2,reps);  % sFtore estimated betahats, r-th repetition in r-th column
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
    b_hat              = numerator/denominator;
    b_hat_all(2,r)     = b_hat;
    b_hat_all(1,r)     = y_avg_r - b_hat*x_avg;
    
    u_hat              = y(:,r) - b_hat_all(1,r) - b_hat_all(2,r)*x;
    x_u_corr_all(1,r)  = corr(x, u_hat);
    
    r = r + 1;
end

standard_dev1=std(b_hat_all(1,:));
standard_dev2=std(b_hat_all(2,:));

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
% fprintf('Alpha:%8.4f',mean(b_hat_all(1,:),2));
% fprintf('  Beta:%8.4f\n',mean(b_hat_all(2,:),2));
% fprintf('Standard errors (standard deviation of estimates at Monte Carlo repetitions)\n');
% fprintf('Alpha:%8.4f',standard_dev1);
% fprintf('  Beta:%8.4f\n',standard_dev2);
% fprintf('Mean correlation of x and u hat.\n');
% fprintf('Correlation:%8.4f\n',mean_x_u_corr);
% fprintf('Standard error of correlation of x and u hat.\n');
% fprintf('Standard error:%8.4f',stdev_x_u_corr);


corr_matrix((sigma_ind-1)*2+1,1) = mean_x_u_corr;
corr_matrix(sigma_ind*2,1) = stdev_x_u_corr;


beta_matrix((sigma_ind-1)*2+1,1) = mean(b_hat_all(2,:),2);
beta_matrix(sigma_ind*2,1) = standard_dev2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FULL PAIRWISE ESTIMATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%
% 1st step %
%%%%%%%%%%%%
% Define variables
% number of betas
number_of_betas = T * (T-1) /2;
% define delta vectors
dz0 = zeros(number_of_betas,reps);
dz_star = zeros(number_of_betas,reps);
dz_star_wave = zeros(number_of_betas,reps);
% define the Mz projection of the delta vectors
M_dz0 =  zeros(number_of_betas,reps);
M_dz_star_wave = zeros(number_of_betas,reps);
%Mz = cell(number_of_betas,number_of_betas,reps);%zeros(number_of_betas,number_of_betas,reps);


pairwise_betas = zeros(2,number_of_betas);
b_hat_all = zeros(2,reps);  % store estimated betahats, r-th repetition in r-th column
x_u_corr_all = zeros(1,reps);
u_hat_mat = zeros(number_of_betas,reps);
mu_z = zeros(1,reps);

x_avg = zeros(number_of_betas,1);
y_avg = zeros(number_of_betas,reps);

r = 1;
while r < reps+0.5
    counter=1;
    % define z vectors
        % iterate over all pairs
    for i=(1:1:T-1)
        for j=((i+1):1:T)
           % if i<j
                % calculate z values
                r_j = (y(i,r) + y(j,r))/2;
                dz0(counter,r) = logistic(y(j,r)) - logistic(y(i,r)) + omega(1,r)*var_y(1,r);
                dz_star(counter,r) = exp(-r_j)*logistic(r_j).^2 *(x(j,1) - x(i,1));
                dz_star_wave(counter,r) = -logistic_prime_prime(r_j)*(x(j,1) - x(i,1));
                %
                x_avg(counter,1)     = (x(i,1)+x(j,1))/2;
                y_avg(counter,r)     = (y(i,r)+y(j,r))/2;
                % step counter
                counter = counter + 1;
        end
    end
    % projection matrix
    %Mz(:,:,r) = eye(number_of_betas) - dz_star(:,r)*inv(dz_star(:,r)'*dz_star(:,r))*dz_star(:,r)';
    Mz = eye(number_of_betas) - dz_star(:,r)*inv(dz_star(:,r)'*dz_star(:,r))*dz_star(:,r)';
    
    Mz_dz0 = Mz*dz0(:,r);
    Mz_dz_star_wave = Mz*dz_star_wave(:,r);
    mu_z(1,r) = mean(Mz_dz_star_wave);
    

                if strcmp(abs_dummy, 'True')
                    delta_x = abs(Mz_dz_star_wave);
                else
                    delta_x = Mz_dz_star_wave;
                end
                % calculate betahat
                
                x_avg     = Mz_dz_star_wave;
                y_avg     = Mz_dz0;
                
                
                numerator = Mz_dz0;
                %delta_y(1, counter) = numerator;
                denominator = Mz_dz_star_wave;
                b_hat_i     = numerator./denominator;
                alpha_hat_i = y_avg - b_hat_i.*x_avg;
                pairwise_betas(1,:) = alpha_hat_i;
                pairwise_betas(2,:) = b_hat_i;



    % Obtain the delta-x weighted average of pairwise betas
            
    sum_delta_x = nansum(delta_x);
    inf_index =(pairwise_betas==Inf);
    minusinf_index = (pairwise_betas==-Inf);
    sum_noninf = max(max(inf_index(1,:),inf_index(2,:)),max(minusinf_index(1,:),minusinf_index(2,:)));
    %sum((inf_index+minusinf_index),1);
    
    noninf_index = (1-sum_noninf)==1; 
    
    weighted_parwise_betas = pairwise_betas(:,noninf_index)*delta_x( noninf_index);

    weighted_average_parwise_betas = weighted_parwise_betas./sum_delta_x;
    
    b_hat_all(1,r)        = weighted_average_parwise_betas(1);
    b_hat_all(2,r)        = weighted_average_parwise_betas(2);
    
    u_hat              = Mz_dz0 - b_hat_all(1,r) - b_hat_all(2,r)*Mz_dz_star_wave;
    u_hat_mat(:,r) = u_hat;
    x_u_corr_all(1,r)  = corr(Mz_dz_star_wave, u_hat);
    
    %mu_z(r,1) = mean(Mz_dz_star_wave);

    r = r + 1;
end
    
 
    
    %%%%%%%%%%%%
    % 2nd step %
    %%%%%%%%%%%%
    
    % test stat
    % delta_u = delta_y - b_hat_all(2,:).*delta_x;
    %test_stat(sigma_ind, r) = (sum(delta_x.*(delta_x*(beta-b_hat_all(2,r))+delta_u))/(T_list(sample_size_ind)*(T_list(sample_size_ind)-1))*2);
    
    % bias correction
    
    d_eta = sum(u_hat_mat,1);
    correction = (1/number_of_betas*mu_z).*d_eta;
    delta2_hat = b_hat_all(2,:) + correction;
    beta1_hat = (delta2_hat).^0.5;
    
    
    % sign of beta 2
    d_beta1_plus =  1/number_of_betas*sum(sum((dz0 - beta1_hat.*dz_star -  dz_star_wave).^2));
    d_beta1_minus =  1/number_of_betas*sum(sum((dz0 + beta1_hat.*dz_star -  dz_star_wave).^2));
    
    if d_beta1_plus > d_beta1_minus 
        beta1_hat = (-1).*beta1_hat;
    end
        
    
    
 
%end

standard_dev1=std(beta1_hat);
standard_dev2=std(delta2_hat);
%standard_dev1=std(b_hat_all(1,:));
%standard_dev2=std(b_hat_all(2,:));

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
corr_matrix((sigma_ind-1)*2+1,2) = mean_x_u_corr;
corr_matrix(sigma_ind*2,2) = stdev_x_u_corr;

%beta_matrix((sigma_ind-1)*2+1,2) = nanmean(b_hat_all(2,:),2);
beta_matrix((sigma_ind-1)*2+1,2) = nanmean(delta2_hat);

%nan_err =sum(isnan(b_hat_all(2,:)));
nan_err =sum(isnan(delta2_hat));
if nan_err ~= 0
   disp(num2str(nan_err))
end
beta_matrix(sigma_ind*2,2) = standard_dev2;

end



if strcmp(abs_dummy,'True')
    if strcmp(dist,'norm')
        beta_title = strcat('Results/Nonlinear_Norm_beta_Abs_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        corr_title = strcat('Results/Nonlinear_Norm_corr_Abs_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        tstat_title = strcat('Results/Nonlinear_Norm_tstat_Abs_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
    else
        beta_title = strcat('Results/Nonlinear_Unif_beta_Abs_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        corr_title = strcat('Results/Nonlinear_Unif_corr_Abs_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        tstat_title = strcat('Results/Nonlinear_Unif_tstat_Abs_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');  
    end 
else
    if strcmp(dist,'norm')
        beta_title = strcat('Results/Nonlinear_Norm_beta_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        corr_title = strcat('Results/Nonlinear_Norm_corr_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        tstat_title = strcat('Results/Nonlinear_Norm_tstat_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
    else
        beta_title = strcat('Results/Nonlinear_Unif_beta_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        corr_title = strcat('Results/Nonlinear_Unif_corr_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        tstat_title = strcat('Results/Nonlinear_Unif_tstat_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');  
    end
end
%Results/
parsave(beta_title, beta_matrix)
parsave(corr_title, corr_matrix)
parsave(tstat_title, test_stat)

end





