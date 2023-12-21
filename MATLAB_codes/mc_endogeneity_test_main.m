%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%prepare WD
clear all
close all
clc

%%
clear all
%clust = parcluster;
%clust.NumWorkers=4;
digits(64)

%%% version = 'A';
% case 1
abs_dummy = 'False'; 
dist = 'unif';
sorted = 'True';
version = 'A';
run('ewpo_testing_fullpairwise')
% case 2
abs_dummy = 'False'; 
dist = 'unif';
sorted = 'False';
version = 'A';
run('ewpo_testing_fullpairwise')
% case 3
abs_dummy = 'True'; 
dist = 'unif';
sorted = 'False';
version = 'A';
run('ewpo_testing_fullpairwise')
% case 4
abs_dummy = 'True'; 
dist = 'unif';
sorted = 'True';
version = 'A';
run('ewpo_testing_fullpairwise')
%%% version = 'B';
% case 1
abs_dummy = 'False'; 
dist = 'unif';
sorted = 'True';
version = 'B';
run('ewpo_testing_fullpairwise')
% case 2
abs_dummy = 'False'; 
dist = 'unif';
sorted = 'False';
version = 'B';
run('ewpo_testing_fullpairwise')
% case 3
abs_dummy = 'True'; 
dist = 'unif';
sorted = 'False';
version = 'B';
run('ewpo_testing_fullpairwise')
% case 4
abs_dummy = 'True'; 
dist = 'unif';
sorted = 'True';
version = 'B';
run('ewpo_testing_fullpairwise')


%sorted = '';
sorted = '_sorted';
%version = 'A';
version = 'B';


%%%%%%%%%
%Uniform%
%%%%%%%%%

%Uniform
T_list = [10, 50, 100, 500, 1000];
%reload and combine data
req_date ='21-Dec-2023';
abs_beta_out =[];
abs_diff_out =[];


noabs_beta_out =[];
noabs_diff_out =[];

for ss =1:length(T_list)
%{
    abs_beta_add = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Unif_beta_version_',version,'_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    abs_diff_add = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Unif_diff_matrix_version_',version,'_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    
    abs_beta_out = [abs_beta_out,abs_beta_add.x'];
    abs_diff_out = [abs_diff_out,abs_diff_add.x'];
    %}
    noabs_beta_add = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Unif_beta_version_',version,'_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    noabs_diff_add = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Unif_diff_matrix_version_',version,'_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    
    noabs_beta_out = [noabs_beta_out,noabs_beta_add.x'];
    noabs_diff_out = [noabs_diff_out,noabs_diff_add.x'];

end


abs_data_table_beta = table(abs_beta_out');
abs_data_table_diff = table(abs_diff_out');


noabs_data_table_beta = table(noabs_beta_out');
noabs_data_table_diff = table(noabs_diff_out');


% % write output
writetable(abs_data_table_beta,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Unif_beta_version',version,'_Abs_deltax_', req_date,'.xlsx'))
writetable(abs_data_table_diff,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Unif_diff_version',version,'_Abs_deltax_', req_date,'.xlsx'))

writetable(noabs_data_table_beta,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Unif_beta_version',version,'_deltax_', req_date,'.xlsx'))
writetable(noabs_data_table_diff,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Unif_diff_version',version,'_deltax_', req_date,'.xlsx'))
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Non-sorted adjacent %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
abs_dummy = 'True';
sorted = 'False';
dist = 'norm';
run('usual_suspects_mc_non_sorted_adjacent_different_alphas')
dist = 'unif';
run('usual_suspects_mc_non_sorted_adjacent_different_alphas')


%%%%%%%%%
%Uniform%
%%%%%%%%%

%Uniform
T_list = [10, 50, 100, 500, 1000];
%reload and combine data
req_date ='14-Dec-2023';
abs_beta_out =[];
abs_diff_out =[];


for ss =1:length(T_list)
    
    abs_beta_add = load(strcat('Results/usual_suspects_nonsorted_adjacent_different_alphas_Unif_beta_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    abs_diff_add = load(strcat('Results/usual_suspects_nonsorted_adjacent_different_alphas_Unif_diff_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    
    abs_beta_out = [abs_beta_out,abs_beta_add.x'];
    abs_diff_out = [abs_diff_out,abs_diff_add.x'];
    
end


abs_data_table_beta = table(abs_beta_out');
abs_data_table_diff = table(abs_diff_out');



% % write output
writetable(abs_data_table_beta,strcat('Results/Excels/usual_suspects_nonsorted_adjacent_different_alphas_Unif_beta_Abs_deltax_', req_date,'.xlsx'))
writetable(abs_data_table_diff,strcat('Results/Excels/usual_suspects_nonsorted_adjacent_different_alphas_Unif_diff_Abs_deltax_', req_date,'.xlsx'))
