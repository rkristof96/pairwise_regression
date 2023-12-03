close all;
clear all;
clc;

%clust = parcluster;
%clust.NumWorkers=4;
clust = parcluster;
clust.NumWorkers=4;
digits(64)

% define date_when the run starts
date_var = datestr(datetime('today'));

significance_level = 0.05;

% parameter sets
 T_list = [50, 500,1000, 5000];
 sigma_list = [0,0.2,0.5,0.8];
 cov_list = [0,0.0391,0.2493,0.6397];
 
 %sigma_list =[0];
 %cov_list = [0];
 
 %0.8 - 0.6397
 %0.5 - 0.2493
 %0.2 - 0.0391


%beta_matrix = zeros(length(sigma_list)*length(T_list)*2,2);

%critical_value_matrix = zeros(length(sigma_list)*length(T_list),2);

abs_dummy = 'True'; % or Abs
dist = 'norm';
 
parfor sample_size_ind = 1:length(T_list)
critical_value_matrix = zeros(2,length(sigma_list));
beta_value_matrix = zeros(1,length(sigma_list));


for sigma_ind = 1:length(sigma_list)

% true parameters

alpha = 1;
beta  = 0.5;
sigma = sigma_list(sigma_ind);%0.8;

b_true = [alpha;beta;sigma];

T = T_list(sample_size_ind); % number of observations
reps = 1;%1000; % number of Monte Carlo repetitions

critical_values = zeros(2,reps);
beta_values = zeros(1,reps);


% explanatory variable
%rand('seed',202101);
rand('seed',20230318);
% generate x: (Tx1) u: (Txreps) vector of bi-variate normal distributed
% random variables

mu = zeros(1,reps+1);
cov_matrix = zeros(reps+1);
cov_matrix(:,:) = cov_list(sigma_ind);

cov_matrix(1,:) = sigma;
cov_matrix(:,1) = sigma;

for diag_i=(1:1:reps+1)
    cov_matrix(diag_i,diag_i) = 1; %ones([(reps+1) 1]); %1;
end

%cov_matrix(1,1) = 5;

%rng('default')  % For reproducibility

%R = mvnrnd(mu,cov_matrix,T);

%x = R(:,1);

if strcmp(dist,'uniform')
    %cov_matrix(1,1) = 1;
    %rng('default')  % For reproducibility
    R = mvnrnd(mu,cov_matrix,T);
    x = R(:,1);
    x = normcdf(x);
    x = unifinv(x,-5,5);
else
    cov_matrix(1,1) = 5;
    %rng('default')  % For reproducibility
    R = mvnrnd(mu,cov_matrix,T);
    x = R(:,1);
end

eps = R(:,2:end);
    


% dependent variables, in each of the repetitions

y = alpha+beta*x+eps;  % (T x reps) matrix of dependent variables

% sort
xy = [x y];

xy = sortrows(xy,1);

x = xy(:,1);
y = xy(:,2:reps+1);


%%%%%%%%%%%%%%
% FULL PAIRWISE ESTIMATION %
%%%%%%%%%%%%%%

% Jack-knife
% d - number of observations t drop

d = fix((sqrt(T)+T)/2);
if d < fix(sqrt(T))
    d = ceil(sqrt(T)); % if d < n^0.5
end

Reps_max =  nchoosek(T,(T-d));
% reps is usually a very large number
% maximize it in 10,000
reps_min = min(Reps_max, 10000);

r = 1;

while r < reps +0.5

b_hat_all = zeros(2,reps_min);  % store estimated betahats, r-th repetition in r-th column

%critical_values = zeros(2,reps);

s = 1;
while s < reps_min+0.5
   % while s < reps_min + 0.5

    
    randsample_index = randsample(T,(T-d));
    
    x_jackknife = x(randsample_index);
    y_jackknife = y(randsample_index,r);
    
    jack_length = length(y_jackknife);
    
    % define objects to store results
    number_of_betas = jack_length * (jack_length-1) /2;
    pairwise_betas = zeros(2,number_of_betas);
    delta_y = zeros(1,number_of_betas);
    delta_x = zeros(1,number_of_betas);
    counter=1;

    % iterate over all pairs
    for i=(1:1:jack_length-1)
        for j=((i+1):1:jack_length)
           % if i<j
                % calculate x difference

                if strcmp(abs_dummy, 'True')
                    delta_x(1, counter) = abs(x_jackknife(j,1) - x_jackknife(i,1));
                else
                    delta_x(1, counter) = x_jackknife(j,1) - x_jackknife(i,1);
                end
                % calculate betahat
                x_avg     = (x_jackknife(i,1)+x_jackknife(j,1))/2;
                %y_avg     = (y_jackknife(i,r)+y_jackknife(j,r))/2;
                y_avg     = (y_jackknife(i,1)+y_jackknife(j,1))/2;
                %numerator = y_jackknife(j,r) - y_jackknife(i,r);
                numerator = y_jackknife(j,1) - y_jackknife(i,1);
                delta_y(1, counter) = numerator;
                denominator = x_jackknife(j,1) - x_jackknife(i,1);
                b_hat_i     = numerator/denominator;
                alpha_hat_i = y_avg - b_hat_i*x_avg;
                pairwise_betas(1,counter)=alpha_hat_i;
                pairwise_betas(2,counter)=b_hat_i;
                counter   = counter+1;
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
    
    b_hat_all(1,s)        = weighted_average_parwise_betas(1);
    b_hat_all(2,s)        = weighted_average_parwise_betas(2);
    
    u_hat              = y_jackknife(:,1) - b_hat_all(1,s) - b_hat_all(2,s)*x_jackknife;
   
    
    % test stat
    delta_u = delta_y - b_hat_all(2,s)*delta_x;
    
    s = s + 1;
    %end
end

b_hat2_sorted = sort(b_hat_all(2,:));

lower_bound = fix(significance_level/2*reps_min);
upper_bound = fix((1-significance_level/2)*reps_min);

%lower limit
critical_values(1,r) = b_hat2_sorted(lower_bound);
critical_values(2,r) = b_hat2_sorted(upper_bound);
beta_values(1,r) = mean(b_hat2_sorted)

   r = r + 1;
    %end
    
    % beta_matrix((sigma_ind-1)*2+1,2) = nanmean(b_hat_all(2,:),2);

    nan_err =sum(isnan(b_hat_all(2,:)));
    if nan_err ~= 0
       disp(num2str(nan_err))
    end
    
end


critical_value_matrix(:,sigma_ind) = mean(critical_values,2);
beta_value_matrix(:,sigma_ind) = mean(beta_values);




end



if strcmp(abs_dummy,'True')
    if strcmp(dist,'norm')
        cv_title = strcat('Results/Jack_knife/Norm_CriticalValues_Abs_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        coef_title = strcat('Results/Jack_knife/Norm_Coef_Abs_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');

    else
        cv_title = strcat('Results/Jack_knife/Unif_CriticalValues_Abs_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        coef_title = strcat('Results/Jack_knife/Unif_Coef_Abs_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');    end 
else
    if strcmp(dist,'norm')
        cv_title = strcat('Results/Jack_knife/Norm_CriticalValues_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        coef_title = strcat('Results/Jack_knife/Norm_Coef_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
    else
        cv_title = strcat('Results/Jack_knife/Unif_CriticalValues_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
        coef_title = strcat('Results/Jack_knife/Unif_Coef_deltax_T_',num2str(T_list(sample_size_ind)),'_',date_var,'.mat');
    end
end
%Results/
parsave(cv_title, critical_value_matrix)
parsave(coef_title, beta_value_matrix)

end





