function retval = lossfunction2_beta1(beta)

%using global variables w and x;
global assigned_weight;
global pairwise_beta1;

weighted_product = pairwise_beta1 .* assigned_weight';
beta_diff = weighted_product - beta;
beta_diff_squared = beta_diff.^2;
sum_of_beta_diff_squared=sum(beta_diff_squared);

retval = sum_of_beta_diff_squared; %define the return value of the function