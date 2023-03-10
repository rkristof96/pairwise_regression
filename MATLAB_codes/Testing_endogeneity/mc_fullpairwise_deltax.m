close all;
clear all;
clc;

% parameter sets
 T_list = [50];%, 500, 5000];
 sigma_list = [0,0.2,0.5,0.8];
 cov_list = [0,0.0071,0.0491,0.1272];
% 0.0071 with sigma = 0.2
% 0.0491 with sigma = 0.5
% 0.1272 with sigma = 0.8

corr_matrix = zeros(length(sigma_list)*length(T_list)*2,2);
beta_matrix = zeros(length(sigma_list)*length(T_list)*2,2);

abs_dummy = 'True';
 
for sample_size_ind = 1:length(T_list)
% corr_matrix = zeros(length(sigma_list)*2,2);
% beta_matrix = zeros(length(sigma_list)*2,2);
for sigma_ind = 1:length(sigma_list)

% true parameters

alpha = 1;
beta  = 0.5;
sigma = sigma_list(sigma_ind);%0.8;

b_true = [alpha;beta;sigma];

T = T_list(sample_size_ind); % number of observations
reps = 1000; % number of Monte Carlo repetitions

% explanatory variable
rand('seed',202101);
% generate x: (Tx1) u: (Txreps) vector of bi-variate normal distributed
% random variables

mu = zeros(1,reps+1);
cov_matrix = zeros(reps+1);
cov_matrix(:,:) = cov_list(sigma_ind); %0.1272;
% 0.0071 with sigma = 0.2
% 0.0491 with sigma = 0.5
% 0.1272 with sigma = 0.8

cov_matrix(1,:) = sigma;
cov_matrix(:,1) = sigma;

for diag_i=(1:1:reps+1)
    cov_matrix(diag_i,diag_i) = 1; %ones([(reps+1) 1]); %1;
end

cov_matrix(1,1) = 5;

rng('default')  % For reproducibility

R = mvnrnd(mu,cov_matrix,T);

x = R(:,1);
eps = R(:,2:end);
    
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
x_u_corr_all = zeros(1,reps);

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

corr_matrix((sample_size_ind-1)*length(T_list)+(sigma_ind-1)*2+1,1) = mean_x_u_corr;
corr_matrix((sample_size_ind-1)*length(T_list)+sigma_ind*2,1) = stdev_x_u_corr;


beta_matrix((sample_size_ind-1)*length(T_list)+(sigma_ind-1)*2+1,1) = mean(b_hat_all(2,:),2);
beta_matrix((sample_size_ind-1)*length(T_list)+sigma_ind*2,1) = standard_dev2;

%%%%%%%%%%%%%%
% FULL PAIRWISE ESTIMATION %
%%%%%%%%%%%%%%

b_hat_all = zeros(2,reps);  % store estimated betahats, r-th repetition in r-th column
x_u_corr_all = zeros(1,reps);

r = 1;
while r < reps+0.5
    number_of_betas = T * (T-1) /2;
    pairwise_betas = zeros(2,number_of_betas);
    delta_u = zeros(1,number_of_betas);
    delta_x = zeros(1,number_of_betas);
    counter=1;

    % iterate over all pairs
    for i=(1:1:T-1)
        for j=((i+1):1:T)
           % if i<j
                % calculate x difference
                if abs_dummy == 'True'
                    delta_x(1, counter) = abs(x(j,1) - x(i,1));
                else
                    delta_x(1, counter) = x(j,1) - x(i,1);
                end
                % calculate betahat
                x_avg     = (x(i,1)+x(j,1))/2;
                y_avg     = (y(i,r)+y(j,r))/2;
                numerator = y(j,r) - y(i,r);
                denominator = x(j,1) - x(i,1);
                b_hat_i     = numerator/denominator;
                alpha_hat_i = y_avg - b_hat_i*x_avg;
                pairwise_betas(1,counter)=alpha_hat_i;
                pairwise_betas(2,counter)=b_hat_i;
                delta_u = numerator;
                counter   = counter+1;
          %  end
        end
    end

    % Obtain the delta-x weighted average of pairwise betas
            
    sum_delta_x = sum(delta_x);
    weighted_parwise_betas = pairwise_betas*delta_x';
    weighted_average_parwise_betas = weighted_parwise_betas./sum_delta_x;
    
    b_hat_all(1,r)        = weighted_average_parwise_betas(1);
    b_hat_all(2,r)        = weighted_average_parwise_betas(2);
    
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


corr_matrix((sample_size_ind-1)*length(T_list)+(sigma_ind-1)*2+1,2) = mean_x_u_corr;
corr_matrix((sample_size_ind-1)*length(T_list)+sigma_ind*2,2) = stdev_x_u_corr;


beta_matrix((sample_size_ind-1)*length(T_list)+(sigma_ind-1)*2+1,2) = mean(b_hat_all(2,:),2);
beta_matrix((sample_size_ind-1)*length(T_list)+sigma_ind*2,2) = standard_dev2;

end
% if abs_dummy == 'True'
%     beta_title = strcat('Results/Norm_beta_abs_deltax_sigma_',num2str(sigma_list(sigma_ind)),'_',today(string('datetime')),'.xlsx');
%     corr_title = strcat('Results/Norm_beta_abs_deltax_sigma_',num2str(sigma_list(sigma_ind)),'_',today(string('datetime')),'.xlsx')
% else
%     strcat('Results/Norm_beta_deltax_sigma_',num2str(sigma_list(sigma_ind)),'_',today(string('datetime')),today(string('datetime')),'.xlsx')
%     strcat('Results/Norm_beta_deltax_sigma_',num2str(sigma_list(sigma_ind)),'_',today(string('datetime')),'.xlsx')
% end 
% % write output
% writetable(beta_matrix,beta_title)
% writetable(corr_matrix,corr_title)
end

if abs_dummy == 'True'
    beta_title = 'Norm_beta_abs_deltax.xlsx';
    corr_title = 'Norm_corr_abs_deltax.xlsx';
else
    beta_title = 'Norm_beta_deltax.xlsx';
    corr_title = 'Norm_corr_deltax.xlsx';
end 


data_table_beta = table(beta_matrix);
data_table_corr = table(corr_matrix);

%data_table_beta.Properties.VariableNames = {'OLS' ,'Pair_wise'};
%data_table_corr.Properties.VariableNames = {'OLS' ,'Pair_wise'};


% write output
writetable(data_table_beta,beta_title)
writetable(data_table_corr,corr_title)