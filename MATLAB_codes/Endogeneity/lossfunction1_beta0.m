function retval = lossfunction1_beta0(beta)

%using global variables w and x;
global assigned_weight;
global pairwise_beta0;

beta_diff = pairwise_beta0 - beta;
weighted_product = beta_diff .* assigned_weight';
weighted_product_squared = weighted_product.^2;
sum_of_weighted_product_squared=sum(weighted_product_squared);

retval = sum_of_weighted_product_squared; %define the return value of the function