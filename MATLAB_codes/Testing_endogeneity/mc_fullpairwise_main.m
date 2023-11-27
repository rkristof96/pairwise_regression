%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%prepare WD
clear all
close all
clc

%% Normal DGP

clear all
clust = parcluster;
clust.NumWorkers=4;
digits(64)
% define parameter
% absolute weighted difference
abs_dummy = 'True'; 
% DGP
dist = 'norm';
run('mc_fullpairwise_deltax_parfor.m')

% define parameter
% simple weighted difference
abs_dummy = 'False'; 
% DGP
dist = 'norm';
run('mc_fullpairwise_deltax_parfor.m')

%% Uniform DGP
clear all
clust = parcluster;
clust.NumWorkers=4;
clc
digits(64)
% define parameter
% absolute weighted difference
abs_dummy = 'True'; 
% DGP
dist = 'uniform';
run('mc_fullpairwise_deltax_parfor.m')

% define parameter
% simple weighted difference
abs_dummy = 'False'; 
% DGP
dist = 'uniform';
run('mc_fullpairwise_deltax_parfor.m')
%%
%Normal
 T_list = [50, 500, 5000];%[50, 500, 5000];
%reload and combine data
req_date ='09-Apr-2023';
abs_beta_out =[];
abs_corr_out =[];
abs_tstat_out = [];

noabs_beta_out =[];
noabs_corr_out =[];
noabs_tstat_out =[];
for ss =1:length(T_list)
    
    abs_beta_add =load(strcat('Results/Norm_beta_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    abs_corr_add = load(strcat('Results/Norm_corr_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    abs_tstat_add = load(strcat('Results/Norm_tstat_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    
    abs_beta_out = [abs_beta_out,abs_beta_add.x'];
    abs_corr_out = [abs_corr_out,abs_corr_add.x'];
    abs_tstat_out = [abs_tstat_out,abs_tstat_add.x'];
    
    noabs_beta_add = load(strcat('Results/Norm_beta_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    noabs_corr_add = load(strcat('Results/Norm_corr_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    noabs_tstat_add = load(strcat('Results/Norm_tstat_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));

        
    noabs_beta_out = [noabs_beta_out,noabs_beta_add.x'];
    noabs_corr_out = [noabs_corr_out,noabs_corr_add.x'];
    noabs_tstat_out = [noabs_tstat_out,noabs_tstat_add.x'];

end


abs_data_table_beta = table(abs_beta_out');
abs_data_table_corr = table(abs_corr_out');
abs_data_table_tstat = table(mean(abs_tstat_out,1)');
abs_data_table_tstat_std = table(std(abs_tstat_out,1)');
abs_data_table_tstat_skewness = table(skewness(abs_tstat_out,1)');
abs_data_table_tstat_kurtosis = table(kurtosis(abs_tstat_out,1)');


abs_data_table_whole_tstat = table(abs_tstat_out);

noabs_data_table_beta = table(noabs_beta_out');
noabs_data_table_corr = table(noabs_corr_out');
noabs_data_table_tstat = table(mean(noabs_tstat_out,1)');
noabs_data_table_tstat_std = table(std(noabs_tstat_out,1)');
noabs_data_table_tstat_skewness = table(skewness(noabs_tstat_out,1)');
noabs_data_table_tstat_kurtosis = table(kurtosis(noabs_tstat_out,1)');
noabs_data_table_whole_tstat = table(noabs_tstat_out);

% % write output
writetable(abs_data_table_beta,strcat('Results/Excels/Norm_abs_beta_', req_date,'.xlsx'))
writetable(abs_data_table_corr,strcat('Results/Excels/Norm_abs_corr_', req_date,'.xlsx'))
writetable(abs_data_table_tstat,strcat('Results/Excels/Norm_abs_tstat_', req_date,'.xlsx'))
writetable(abs_data_table_tstat_std,strcat('Results/Excels/Norm_abs_tstat_std_', req_date,'.xlsx'))
writetable(abs_data_table_tstat_skewness,strcat('Results/Excels/Norm_abs_tstat_skewness_', req_date,'.xlsx'))
writetable(abs_data_table_tstat_kurtosis,strcat('Results/Excels/Norm_abs_tstat_kurtosis_', req_date,'.xlsx'))
writetable(abs_data_table_whole_tstat,strcat('Results/Excels/Norm_abs_whole_tstat_', req_date,'.xlsx'))


writetable(noabs_data_table_beta,strcat('Results/Excels/Norm_noabs_beta_', req_date,'.xlsx'))
writetable(noabs_data_table_corr,strcat('Results/Excels/Norm_noabs_corr_', req_date,'.xlsx'))
writetable(noabs_data_table_tstat,strcat('Results/Excels/Norm_noabs_tstat_', req_date,'.xlsx'))
writetable(noabs_data_table_tstat_std,strcat('Results/Excels/Norm_noabs_tstat_std_', req_date,'.xlsx'))
writetable(noabs_data_table_tstat_skewness,strcat('Results/Excels/Norm_noabs_tstat_skewness_', req_date,'.xlsx'))
writetable(noabs_data_table_tstat_kurtosis,strcat('Results/Excels/Norm_noabs_tstat_kurtosis_', req_date,'.xlsx'))
writetable(noabs_data_table_whole_tstat,strcat('Results/Excels/Norm_noabs_whole_tstat_', req_date,'.xlsx'))

%%
%%%%%%%%%%%%%%%
% Uniform DGP %
%%%%%%%%%%%%%%%

%reload and combine data
T_list = [50, 500, 5000];%[50, 500, 5000];
req_date ='09-Apr-2023';
abs_beta_out =[];
abs_corr_out =[];
abs_tstat_out = [];

noabs_beta_out =[];
noabs_corr_out =[];
noabs_tstat_out =[];
for ss =1:length(T_list)
   %%{ 
    abs_beta_add =load(strcat('Results/Unif_beta_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    abs_corr_add = load(strcat('Results/Unif_corr_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    abs_tstat_add = load(strcat('Results/Unif_tstat_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    
    
    abs_beta_out = [abs_beta_out,abs_beta_add.x'];
    abs_corr_out = [abs_beta_out,abs_corr_add.x'];
    abs_tstat_out = [abs_tstat_out,abs_tstat_add.x'];
    
    %}
    %%{
    noabs_beta_add = load(strcat('Results/Unif_beta_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    noabs_corr_add = load(strcat('Results/Unif_corr_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    noabs_tstat_add = load(strcat('Results/Unif_tstat_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));

        
    noabs_beta_out = [noabs_beta_out,noabs_beta_add.x'];
    noabs_corr_out = [noabs_corr_out,noabs_corr_add.x'];
    noabs_tstat_out = [noabs_tstat_out,noabs_tstat_add.x'];
%%}
end


abs_data_table_beta = table(abs_beta_out');
abs_data_table_corr = table(abs_corr_out');
abs_data_table_tstat = table(mean(abs_tstat_out,1)');

abs_data_table_tstat_std = table(std(abs_tstat_out,1)');
abs_data_table_tstat_skewness = table(skewness(abs_tstat_out,1)');
abs_data_table_tstat_kurtosis = table(kurtosis(abs_tstat_out,1)');

abs_data_table_whole_tstat = table(abs_tstat_out);

noabs_data_table_beta = table(noabs_beta_out');
noabs_data_table_corr = table(noabs_corr_out');

noabs_data_table_tstat = table(mean(noabs_tstat_out,1)');
noabs_data_table_tstat_std = table(std(noabs_tstat_out,1)');
noabs_data_table_tstat_skewness = table(skewness(noabs_tstat_out,1)');
noabs_data_table_tstat_kurtosis = table(kurtosis(noabs_tstat_out,1)');

noabs_data_table_whole_tstat = table(noabs_tstat_out);

% % write output
writetable(abs_data_table_beta,strcat('Results/Excels/Unif_abs_beta_', req_date,'.xlsx'))
writetable(abs_data_table_corr,strcat('Results/Excels/Unif_abs_corr_', req_date,'.xlsx'))
writetable(abs_data_table_tstat,strcat('Results/Excels/Unif_abs_tstat_', req_date,'.xlsx'))
writetable(abs_data_table_tstat_std,strcat('Results/Excels/Unif_abs_tstat_std_', req_date,'.xlsx'))
writetable(abs_data_table_tstat_skewness,strcat('Results/Excels/Unif_abs_tstat_skewness_', req_date,'.xlsx'))
writetable(abs_data_table_tstat_kurtosis,strcat('Results/Excels/Unif_abs_tstat_kurtosis_', req_date,'.xlsx'))
writetable(abs_data_table_whole_tstat,strcat('Results/Excels/Unif_abs_whole_tstat_', req_date,'.xlsx'))

writetable(noabs_data_table_beta,strcat('Results/Excels/Unif_noabs_beta_', req_date,'.xlsx'))
writetable(noabs_data_table_corr,strcat('Results/Excels/Unif_noabs_corr_', req_date,'.xlsx'))
writetable(noabs_data_table_tstat,strcat('Results/Excels/Unif_noabs_tstat_', req_date,'.xlsx'))
writetable(noabs_data_table_tstat_std,strcat('Results/Excels/Unif_noabs_tstat_std_', req_date,'.xlsx'))
writetable(noabs_data_table_tstat_skewness,strcat('Results/Excels/Unif_noabs_tstat_skewness_', req_date,'.xlsx'))
writetable(noabs_data_table_tstat_kurtosis,strcat('Results/Excels/Unif_noabs_tstat_kurtosis_', req_date,'.xlsx'))
writetable(noabs_data_table_whole_tstat,strcat('Results/Excels/Unif_noabs_whole_tstat_', req_date,'.xlsx'))


%------------------------------------------------------------------------------------------------------------------------------
% Usual suspects
% Exercise 1

%%%%%%%%%%%%%%%%%%%%%%%%
%%% MC Full-pairwise %%%
%%%%%%%%%%%%%%%%%%%%%%%%
% case 1
abs_dummy = 'True'; 
dist = 'norm';
sorted = 'False';
run('usual_suspects_mc_fullpairwise_weighted_parfor')
% case 2
abs_dummy = 'True'; 
dist = 'norm';
sorted = 'True';
run('usual_suspects_mc_fullpairwise_weighted_parfor')
% case 3
abs_dummy = 'False'; 
dist = 'norm';
sorted = 'False';
run('usual_suspects_mc_fullpairwise_weighted_parfor')
% case 4
abs_dummy = 'False'; 
dist = 'norm';
sorted = 'True';
run('usual_suspects_mc_fullpairwise_weighted_parfor')
%%%
%Unif
% case 5
abs_dummy = 'False'; 
dist = 'unif';
sorted = 'True';
run('usual_suspects_mc_fullpairwise_weighted_parfor')
% case 6
abs_dummy = 'False'; 
dist = 'unif';
sorted = 'False';
run('usual_suspects_mc_fullpairwise_weighted_parfor')
% case 7
abs_dummy = 'True'; 
dist = 'unif';
sorted = 'False';
run('usual_suspects_mc_fullpairwise_weighted_parfor')
% case 8
abs_dummy = 'True'; 
dist = 'unif';
sorted = 'True';
run('usual_suspects_mc_fullpairwise_weighted_parfor')


%%
%sorted = '';
sorted = '_sorted' ;%'_sorted';

%%%%%%%%
%Normal%
%%%%%%%%

T_list = [10, 50, 100, 500, 1000];%[10, 50, 100, 500, 1000]

%reload and combine data
req_date ='26-Nov-2023';
abs_beta_out =[];
abs_bias_out =[];


noabs_beta_out =[];
noabs_bias_out =[];

for ss =1:length(T_list)
    
    abs_beta_add = load(strcat('Results/usual_suspects_fullpairwise',sorted,'_Norm_beta_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    abs_bias_add = load(strcat('Results/usual_suspects_fullpairwise',sorted,'_Norm_bias_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    
    abs_beta_out = [abs_beta_out,abs_beta_add.x'];
    abs_bias_out = [abs_bias_out,abs_bias_add.x'];
    
    noabs_beta_add = load(strcat('Results/usual_suspects_fullpairwise',sorted,'_Norm_beta_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    noabs_bias_add = load(strcat('Results/usual_suspects_fullpairwise',sorted,'_Norm_bias_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    
    noabs_beta_out = [abs_beta_out,abs_beta_add.x'];
    noabs_bias_out = [abs_bias_out,abs_bias_add.x'];

end


abs_data_table_beta = table(abs_beta_out');
abs_data_table_bias = table(abs_bias_out');


noabs_data_table_beta = table(noabs_beta_out');
noabs_data_table_bias = table(noabs_bias_out');


% % write output
writetable(abs_data_table_beta,strcat('Results/Excels/usual_suspects_fullpairwise',sorted,'_Norm_beta_Abs_deltax_', req_date,'.xlsx'))
writetable(abs_data_table_bias,strcat('Results/Excels/usual_suspects_fullpairwise',sorted,'_Norm_bias_Abs_deltax_', req_date,'.xlsx'))

writetable(noabs_data_table_beta,strcat('Results/Excels/usual_suspects_fullpairwise',sorted,'_Norm_beta_deltax_', req_date,'.xlsx'))
writetable(noabs_data_table_bias,strcat('Results/Excels/usual_suspects_fullpairwise',sorted,'_Norm_bias_deltax_', req_date,'.xlsx'))

%%%%%%%%%
%Uniform%
%%%%%%%%%

%Uniform
T_list = [10, 50, 100, 500, 1000];%[10, 50, 100, 500, 1000]
%reload and combine data
req_date ='26-Nov-2023';
abs_beta_out =[];
abs_bias_out =[];


noabs_beta_out =[];
noabs_bias_out =[];

for ss =1:length(T_list)
    
    abs_beta_add = load(strcat('Results/usual_suspects_fullpairwise',sorted,'_Unif_beta_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    abs_bias_add = load(strcat('Results/usual_suspects_fullpairwise',sorted,'_Unif_bias_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    
    abs_beta_out = [abs_beta_out,abs_beta_add.x'];
    abs_bias_out = [abs_bias_out,abs_bias_add.x'];
    
    noabs_beta_add = load(strcat('Results/usual_suspects_fullpairwise',sorted,'_Unif_beta_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    noabs_bias_add = load(strcat('Results/usual_suspects_fullpairwise',sorted,'_Unif_bias_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    
    noabs_beta_out = [abs_beta_out,abs_beta_add.x'];
    noabs_bias_out = [abs_bias_out,abs_bias_add.x'];

end


abs_data_table_beta = table(abs_beta_out');
abs_data_table_bias = table(abs_bias_out');


noabs_data_table_beta = table(noabs_beta_out');
noabs_data_table_bias = table(noabs_bias_out');


% % write output
writetable(abs_data_table_beta,strcat('Results/Excels/usual_suspects_fullpairwise',sorted,'_Unif_beta_Abs_deltax_', req_date,'.xlsx'))
writetable(abs_data_table_bias,strcat('Results/Excels/usual_suspects_fullpairwise',sorted,'_Unif_bias_Abs_deltax_', req_date,'.xlsx'))

writetable(noabs_data_table_beta,strcat('Results/Excels/usual_suspects_fullpairwise',sorted,'_Unif_beta_deltax_', req_date,'.xlsx'))
writetable(noabs_data_table_bias,strcat('Results/Excels/usual_suspects_fullpairwise',sorted,'_Unif_bias_deltax_', req_date,'.xlsx'))
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Non-sorted adjacent %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
abs_dummy = 'True';
sorted = 'False';
dist = 'norm';
run('usual_suspects_mc_non_sorted_adjacent_parfor')
dist = 'unif';
run('usual_suspects_mc_non_sorted_adjacent_parfor')

%%%%%%%%
%Normal%
%%%%%%%%

sorted = '';
T_list = [10, 50, 100, 500, 1000];%[10, 50, 100, 500, 1000]

%reload and combine data
req_date ='26-Nov-2023';
abs_beta_out =[];
abs_bias_out =[];



for ss =1:length(T_list)
    
    abs_beta_add = load(strcat('Results/usual_suspects_nonsorted_adjacent_Norm_beta_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    abs_bias_add = load(strcat('Results/usual_suspects_nonsorted_adjacent_Norm_bias_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    
    abs_beta_out = [abs_beta_out,abs_beta_add.x'];
    abs_bias_out = [abs_bias_out,abs_bias_add.x'];

end


abs_data_table_beta = table(abs_beta_out');
abs_data_table_bias = table(abs_bias_out');



% % write output
writetable(abs_data_table_beta,strcat('Results/Excels/usual_suspects_nonsorted_adjacent_Norm_beta_Abs_deltax_', req_date,'.xlsx'))
writetable(abs_data_table_bias,strcat('Results/Excels/usual_suspects_nonsorted_adjacent_Norm_bias_Abs_deltax_', req_date,'.xlsx'))


%%%%%%%%%
%Uniform%
%%%%%%%%%

%Uniform
T_list = [10, 50, 100, 500, 1000];%[10, 50, 100, 500, 1000]
%reload and combine data
req_date ='26-Nov-2023';
abs_beta_out =[];
abs_bias_out =[];


for ss =1:length(T_list)
    
    abs_beta_add = load(strcat('Results/usual_suspects_nonsorted_adjacent_Unif_beta_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    abs_bias_add = load(strcat('Results/usual_suspects_nonsorted_adjacent_Unif_bias_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    
    abs_beta_out = [abs_beta_out,abs_beta_add.x'];
    abs_bias_out = [abs_bias_out,abs_bias_add.x'];
    
end


abs_data_table_beta = table(abs_beta_out');
abs_data_table_bias = table(abs_bias_out');



% % write output
writetable(abs_data_table_beta,strcat('Results/Excels/usual_suspects_nonsorted_adjacent_Unif_beta_Abs_deltax_', req_date,'.xlsx'))
writetable(abs_data_table_bias,strcat('Results/Excels/usual_suspects_nonsorted_adjacent_Unif_bias_Abs_deltax_', req_date,'.xlsx'))

