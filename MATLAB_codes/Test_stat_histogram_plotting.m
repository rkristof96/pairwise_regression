%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test-statistic histograms %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This code generates histograms about the test statistics for different
% DGPs and different correlation between x and u.

clear all
close all
clc
%%
% Normal DGP
% sample size
sample_size_list = [50, 500, 5000];

% iterate over on dataset
for index = 1:length(sample_size_list)
filename=sprintf('Norm_tstat_deltax_T_%d_09-Apr-2023',sample_size_list(index));
% read datasets
st = load(['Results/',filename,'.mat']);
% close figures
close all
f = figure();
% add title
title('Test-statistics distribution, DGP = Normal(0,5), estimator: full pairwise, weights = \Delta x')
% picture parameters
x0=10;
y0=10;
width=750;
height=500;
set(f,'position',[x0,y0,width,height])

% no correlation
f(1)=subplot(2,2,1);
% plot histogram
hist(st.x(1,:))
hold on;
% add mean value
line([mean(st.x(1,:)), mean(st.x(1,:))], [0 300], 'LineWidth', 2, 'Color', 'r');
% title and axis labels
title(f(1),'corr(x,u) = 0')
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
% add title
title(f(2),'corr(x,u) = 0.2')
xlabel(f(2),'Value')
ylabel(f(2),'Frequency')
hold off;

% corr(x,u) = 0.5
f(3)=subplot(2,2,3);
% plot histogram
hist(st.x(3,:))
hold on;
% add mean value
line([mean(st.x(3,:)), mean(st.x(3,:))], [0 300], 'LineWidth', 2, 'Color', 'r');
title(f(3),'corr(x,u) = 0.5')
% labels and title
xlabel(f(3),'Value')
ylabel(f(3),'Frequency')
hold off;

% corr(x,u) = 0.8
f(4)=subplot(2,2,4);
hist(st.x(4,:))
hold on;
% add mean value
line([mean(st.x(4,:)), mean(st.x(4,:))], [0 300], 'LineWidth', 2, 'Color', 'r');
% title and labels
title(f(4),'corr(x,u) = 0.8')
xlabel(f(4),'Value')
ylabel(f(4),'Frequency')
hold off;
% histogram names
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


%%
% Uniform DGP
% sample sizes
sample_size_list = [50, 500, 5000];
for index = 1:length(sample_size_list)
filename=sprintf('Unif_tstat_deltax_T_%d_09-Apr-2023',sample_size_list(index));
% load uniform data
st = load(['Results/',filename,'.mat']);

close all
% define figure
f = figure();
% title
title('Test-statistics distribution, DGP = U(-5,5), estimator: full pairwise, weights = \Delta x')
% figure parameetrs
x0=10;
y0=10;
width=750;
height=500;
set(f,'position',[x0,y0,width,height])

% corr(x,u) = 0
f(1)=subplot(2,2,1);
% plot histogram
hist(st.x(1,:))
hold on;
% add mean value
line([mean(st.x(1,:)), mean(st.x(1,:))], [0 300], 'LineWidth', 2, 'Color', 'r');
% add title and labels
title(f(1),'corr(x,u) = 0')
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
title(f(2),'corr(x,u) = 0.2')
xlabel(f(2),'Value')
ylabel(f(2),'Frequency')
hold off;

% corr(x,u) = 0.5
f(3)=subplot(2,2,3);
hist(st.x(3,:))
hold on;
line([mean(st.x(3,:)), mean(st.x(3,:))], [0 300], 'LineWidth', 2, 'Color', 'r');
title(f(3),'corr(x,u) = 0.5')
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
title(f(4),'corr(x,u) = 0.8')
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