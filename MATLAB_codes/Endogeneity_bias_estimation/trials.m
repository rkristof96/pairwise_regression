n = 1000; rho = .7;
Z = mvnrnd([0 0],[1 rho; rho 1],n);
U = normcdf(Z);
X = [U(:,1) norminv(U(:,2))];

% draw the scatter plot of data with histograms 
figure
scatterhist(X(:,1),X(:,2),'Direction','out')

x = X(:,1);
u = X(:,2);

z1 = Z(:,1);
z2 = Z(:,2);

corr_x_u = corrcoef(x,u);
corr_z1_z2 = corrcoef(z1,z2);