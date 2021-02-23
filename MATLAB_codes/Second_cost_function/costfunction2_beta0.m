function retval = costfunction2_beta0(beta)

%using global variables w and x;
global assigned_weight;
global pairwise_beta0;

weighted_product = pairwise_beta0 .* assigned_weight';
beta_diff = weighted_product - beta;
beta_diff_squared = beta_diff.^2;
sum_of_beta_diff_squared=sum(beta_diff_squared);

retval = sum_of_beta_diff_squared; %define the return value of the function