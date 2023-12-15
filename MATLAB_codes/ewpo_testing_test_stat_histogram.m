%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test-statistic histograms %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This code generates histograms about the test statistics for different
% DGPs and different correlation between x and u.

clear all
close all
clc

%%
% Uniform DGP

req_date = '14-Dec-2023';
sorted =  '_sorted';%''
version = 'A';
abs = ''; %'_Abs';
%%
% (alpha_OLS-alpha_EWPO)
% sample sizes
sample_size_list = [10, 50, 100, 500, 1000];
for index = 1:length(sample_size_list)
filename = strcat('ewpo_testing_fullpairwise',sorted,'_Unif_diff_alpha_version_',version,abs,'_deltax_T_',num2str(sample_size_list(index)),'_',req_date);
%filename=sprintf('Unif_tstat_deltax_T_%d_09-Apr-2023',sample_size_list(index));
% load uniform data
st = load(['Results/',filename,'.mat']);

close all
% define figure
f = figure();
% title
if strcmp(sorted,'_sorted')
    sorted_text = 'sorted,';
else
    sorted_text = '';
end

if strcmp(abs,'')
    weight_text = '\Delta x';
else
    weight_text = '|\Delta x|';
end

title_name = strcat('(\Alpha_OLS-\Alpha_EWPO) distribution, DGP = U(-2,3), estimator: full pairwise, ',sorted_text,' version: ',version,'  ', weight_text);
title(title_name);
%title(' distribution, DGP = U(-2,3), estimator: full pairwise, weights = \Delta x')
% figure parameetrs
x0=10;
y0=10;
width=750;
height = 500;
set(f,'position',[x0,y0,width,height])

% corr(x,u) = 0
f(1)=subplot(2,2,1);
% plot histogram
hist(st.x(1,:))
hold on;
% add mean value
line([mean(st.x(1,:)), mean(st.x(1,:))], [0 300], 'LineWidth', 2, 'Color', 'r');
% add title and labels
title(f(1),'\rho = 0')
xlabel(f(1),'Value')
ylabel(f(1),'Frequency')
hold off;

% corr(x,u) = 0.2
f(2)=subplot(2,2,2);
% plot histogram
hist(st.x(2,:))
hold on;
% add mean value
line([mean(st.x(2,:)), mean(st.x(2,:))], [0 300], 'LineWidth', 2, 'Color', 'r');
% add title and labels
title(f(2),'\rho = 0.2')
xlabel(f(2),'Value')
ylabel(f(2),'Frequency')
hold off;

% corr(x,u) = 0.5
f(3)=subplot(2,2,3);
hist(st.x(3,:))
hold on;
line([mean(st.x(3,:)), mean(st.x(3,:))], [0 300], 'LineWidth', 2, 'Color', 'r');
title(f(3),'\rho = 0.5')
xlabel(f(3),'Value')
ylabel(f(3),'Frequency')
hold off;

% corr(x,u) = 0.8
f(4)=subplot(2,2,4);
% plot histogram
hist(st.x(4,:))
hold on;
% add mean value
line([mean(st.x(4,:)), mean(st.x(4,:))], [0 300], 'LineWidth', 2, 'Color', 'r');
% add title and labels
title(f(4),'\rho = 0.8')
xlabel(f(4),'Value')
ylabel(f(4),'Frequency')
hold off;

hist_name = strcat('Plots/Hist_',filename,'.png');

  try
     saveas(gcf,hist_name);
   catch ME
     fprintf('There was a problem saving the file.\n');
     fprintf('List of open files is:\n')
     fopen('all')
     fprintf('End of list of open files\n');
     rethrow(ME)
   end


end
%%
% (alpha_OLS-alpha_EWPO)^2
for index = 1:length(sample_size_list)
filename = strcat('ewpo_testing_fullpairwise',sorted,'_Unif_diff_alpha2_version_',version,abs,'_deltax_T_',num2str(sample_size_list(index)),'_',req_date);
%filename=sprintf('Unif_tstat_deltax_T_%d_09-Apr-2023',sample_size_list(index));
% load uniform data
st = load(['Results/',filename,'.mat']);

close all
% define figure
f = figure();
% title
if strcmp(sorted,'_sorted')
    sorted_text = 'sorted,';
else
    sorted_text = '';
end

if strcmp(abs,'')
    weight_text = '\Delta x';
else
    weight_text = '|\Delta x|';
end

title_name = strcat('(\Alpha_OLS-\Alpha_EWPO)^2 distribution, DGP = U(-2,3), estimator: full pairwise, ',sorted_text,' version ',version,'  ', weight_text);
title(title_name);
%title(' distribution, DGP = U(-2,3), estimator: full pairwise, weights = \Delta x')
% figure parameetrs
x0=10;
y0=10;
width=750;
height = 500;
set(f,'position',[x0,y0,width,height])

% corr(x,u) = 0
f(1)=subplot(2,2,1);
% plot histogram
hist(st.x(1,:))
hold on;
% add mean value
line([mean(st.x(1,:)), mean(st.x(1,:))], [0 300], 'LineWidth', 2, 'Color', 'r');
% add title and labels
title(f(1),'\rho = 0')
xlabel(f(1),'Value')
ylabel(f(1),'Frequency')
hold off;

% corr(x,u) = 0.2
f(2)=subplot(2,2,2);
% plot histogram
hist(st.x(2,:))
hold on;
% add mean value
line([mean(st.x(2,:)), mean(st.x(2,:))], [0 300], 'LineWidth', 2, 'Color', 'r');
% add title and labels
title(f(2),'\rho = 0.2')
xlabel(f(2),'Value')
ylabel(f(2),'Frequency')
hold off;

% corr(x,u) = 0.5
f(3)=subplot(2,2,3);
hist(st.x(3,:))
hold on;
line([mean(st.x(3,:)), mean(st.x(3,:))], [0 300], 'LineWidth', 2, 'Color', 'r');
title(f(3),'\rho = 0.5')
xlabel(f(3),'Value')
ylabel(f(3),'Frequency')
hold off;

% corr(x,u) = 0.8
f(4)=subplot(2,2,4);
% plot histogram
hist(st.x(4,:))
hold on;
% add mean value
line([mean(st.x(4,:)), mean(st.x(4,:))], [0 300], 'LineWidth', 2, 'Color', 'r');
% add title and labels
title(f(4),'\rho = 0.8')
xlabel(f(4),'Value')
ylabel(f(4),'Frequency')
hold off;

hist_name = strcat('Plots/Hist_',filename,'.png');

  try
     saveas(gcf,hist_name);
   catch ME
     fprintf('There was a problem saving the file.\n');
     fprintf('List of open files is:\n')
     fopen('all')
     fprintf('End of list of open files\n');
     rethrow(ME)
   end


end
%%
% (beta_OLS-beta_EWPO)
for index = 1:length(sample_size_list)
filename = strcat('ewpo_testing_fullpairwise',sorted,'_Unif_diff_beta_version_',version,abs,'_deltax_T_',num2str(sample_size_list(index)),'_',req_date);
%filename=sprintf('Unif_tstat_deltax_T_%d_09-Apr-2023',sample_size_list(index));
% load uniform data
st = load(['Results/',filename,'.mat']);

close all
% define figure
f = figure();
% title
if strcmp(sorted,'_sorted')
    sorted_text = 'sorted,';
else
    sorted_text = '';
end

if strcmp(abs,'')
    weight_text = '\Delta x';
else
    weight_text = '|\Delta x|';
end

title_name = strcat('(\Alpha_OLS-\Alpha_EWPO)^2 distribution, DGP = U(-2,3), estimator: full pairwise, ',sorted_text,' version ',version,'  ', weight_text);
title(title_name);
%title(' distribution, DGP = U(-2,3), estimator: full pairwise, weights = \Delta x')
% figure parameetrs
x0=10;
y0=10;
width=750;
height = 500;
set(f,'position',[x0,y0,width,height])

% corr(x,u) = 0
f(1)=subplot(2,2,1);
% plot histogram
hist(st.x(1,:))
hold on;
% add mean value
line([mean(st.x(1,:)), mean(st.x(1,:))], [0 300], 'LineWidth', 2, 'Color', 'r');
% add title and labels
title(f(1),'\rho = 0')
xlabel(f(1),'Value')
ylabel(f(1),'Frequency')
hold off;

% corr(x,u) = 0.2
f(2)=subplot(2,2,2);
% plot histogram
hist(st.x(2,:))
hold on;
% add mean value
line([mean(st.x(2,:)), mean(st.x(2,:))], [0 300], 'LineWidth', 2, 'Color', 'r');
% add title and labels
title(f(2),'\rho = 0.2')
xlabel(f(2),'Value')
ylabel(f(2),'Frequency')
hold off;

% corr(x,u) = 0.5
f(3)=subplot(2,2,3);
hist(st.x(3,:))
hold on;
line([mean(st.x(3,:)), mean(st.x(3,:))], [0 300], 'LineWidth', 2, 'Color', 'r');
title(f(3),'\rho = 0.5')
xlabel(f(3),'Value')
ylabel(f(3),'Frequency')
hold off;

% corr(x,u) = 0.8
f(4)=subplot(2,2,4);
% plot histogram
hist(st.x(4,:))
hold on;
% add mean value
line([mean(st.x(4,:)), mean(st.x(4,:))], [0 300], 'LineWidth', 2, 'Color', 'r');
% add title and labels
title(f(4),'\rho = 0.8')
xlabel(f(4),'Value')
ylabel(f(4),'Frequency')
hold off;

hist_name = strcat('Plots/Hist_',filename,'.png');

  try
     saveas(gcf,hist_name);
   catch ME
     fprintf('There was a problem saving the file.\n');
     fprintf('List of open files is:\n')
     fopen('all')
     fprintf('End of list of open files\n');
     rethrow(ME)
   end


end
%%
% (beta_OLS-beta_EWPO)^2
for index = 1:length(sample_size_list)
filename = strcat('ewpo_testing_fullpairwise',sorted,'_Unif_diff_beta2_version_',version,abs,'_deltax_T_',num2str(sample_size_list(index)),'_',req_date);
%filename=sprintf('Unif_tstat_deltax_T_%d_09-Apr-2023',sample_size_list(index));
% load uniform data
st = load(['Results/',filename,'.mat']);

close all
% define figure
f = figure();
% title
if strcmp(sorted,'_sorted')
    sorted_text = 'sorted,';
else
    sorted_text = '';
end

if strcmp(abs,'')
    weight_text = '\Delta x';
else
    weight_text = '|\Delta x|';
end

title_name = strcat('(\Alpha_OLS-\Alpha_EWPO)^2 distribution, DGP = U(-2,3), estimator: full pairwise, ',sorted_text,' version ',version,'  ', weight_text);
title(title_name);
%title(' distribution, DGP = U(-2,3), estimator: full pairwise, weights = \Delta x')
% figure parameetrs
x0=10;
y0=10;
width=750;
height = 500;
set(f,'position',[x0,y0,width,height])

% corr(x,u) = 0
f(1)=subplot(2,2,1);
% plot histogram
hist(st.x(1,:))
hold on;
% add mean value
line([mean(st.x(1,:)), mean(st.x(1,:))], [0 300], 'LineWidth', 2, 'Color', 'r');
% add title and labels
title(f(1),'\rho = 0')
xlabel(f(1),'Value')
ylabel(f(1),'Frequency')
hold off;

% corr(x,u) = 0.2
f(2)=subplot(2,2,2);
% plot histogram
hist(st.x(2,:))
hold on;
% add mean value
line([mean(st.x(2,:)), mean(st.x(2,:))], [0 300], 'LineWidth', 2, 'Color', 'r');
% add title and labels
title(f(2),'\rho = 0.2')
xlabel(f(2),'Value')
ylabel(f(2),'Frequency')
hold off;

% corr(x,u) = 0.5
f(3)=subplot(2,2,3);
hist(st.x(3,:))
hold on;
line([mean(st.x(3,:)), mean(st.x(3,:))], [0 300], 'LineWidth', 2, 'Color', 'r');
title(f(3),'\rho = 0.5')
xlabel(f(3),'Value')
ylabel(f(3),'Frequency')
hold off;

% corr(x,u) = 0.8
f(4)=subplot(2,2,4);
% plot histogram
hist(st.x(4,:))
hold on;
% add mean value
line([mean(st.x(4,:)), mean(st.x(4,:))], [0 300], 'LineWidth', 2, 'Color', 'r');
% add title and labels
title(f(4),'\rho = 0.8')
xlabel(f(4),'Value')
ylabel(f(4),'Frequency')
hold off;

hist_name = strcat('Plots/Hist_',filename,'.png');

  try
     saveas(gcf,hist_name);
   catch ME
     fprintf('There was a problem saving the file.\n');
     fprintf('List of open files is:\n')
     fopen('all')
     fprintf('End of list of open files\n');
     rethrow(ME)
   end


end
close all