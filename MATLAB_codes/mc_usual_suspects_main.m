%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%prepare WD
clear all
close all
clc
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
req_date ='29-Nov-2023';
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
req_date ='29-Nov-2023';
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
req_date ='29-Nov-2023';
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
req_date ='29-Nov-2023';
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

% -----------------------------------------------------
% Exercise 2 and version B

%%%%%%%%%%%%%%%%%%%%%%%%
%%% MC Full-pairwise %%%
%%%%%%%%%%%%%%%%%%%%%%%%
% case 1
abs_dummy = 'True'; 
dist = 'norm';
sorted = 'False';
run('usual_suspects_mc_fullpairwise_different_alphas')
% case 2
abs_dummy = 'True'; 
dist = 'norm';
sorted = 'True';
run('usual_suspects_mc_fullpairwise_different_alphas')
% case 3
abs_dummy = 'False'; 
dist = 'norm';
sorted = 'False';
run('usual_suspects_mc_fullpairwise_different_alphas')
% case 4
abs_dummy = 'False'; 
dist = 'norm';
sorted = 'True';
run('usual_suspects_mc_fullpairwise_different_alphas')
%%%
%Unif
% case 5
abs_dummy = 'False'; 
dist = 'unif';
sorted = 'True';
run('usual_suspects_mc_fullpairwise_different_alphas')
% case 6
abs_dummy = 'False'; 
dist = 'unif';
sorted = 'False';
run('usual_suspects_mc_fullpairwise_different_alphas')
% case 7
abs_dummy = 'True'; 
dist = 'unif';
sorted = 'False';
run('usual_suspects_mc_fullpairwise_different_alphas')
% case 8
abs_dummy = 'True'; 
dist = 'unif';
sorted = 'True';
run('usual_suspects_mc_fullpairwise_different_alphas')


%%
%sorted = '';
sorted = '_sorted' ;%'_sorted';

%%%%%%%%
%Normal%
%%%%%%%%

T_list = [10, 50, 100, 500, 1000];%[10, 50, 100, 500, 1000]

%reload and combine data
req_date ='29-Nov-2023';
abs_beta_out =[];
abs_bias_out =[];


noabs_beta_out =[];
noabs_bias_out =[];

for ss =1:length(T_list)
    
    abs_beta_add = load(strcat('Results/usual_suspects_fullpairwise_different_alphas',sorted,'_Norm_beta_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    abs_bias_add = load(strcat('Results/usual_suspects_fullpairwise_different_alphas',sorted,'_Norm_bias_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    
    abs_beta_out = [abs_beta_out,abs_beta_add.x'];
    abs_bias_out = [abs_bias_out,abs_bias_add.x'];
    
    noabs_beta_add = load(strcat('Results/usual_suspects_fullpairwise_different_alphas',sorted,'_Norm_beta_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    noabs_bias_add = load(strcat('Results/usual_suspects_fullpairwise_different_alphas',sorted,'_Norm_bias_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    
    noabs_beta_out = [abs_beta_out,abs_beta_add.x'];
    noabs_bias_out = [abs_bias_out,abs_bias_add.x'];

end


abs_data_table_beta = table(abs_beta_out');
abs_data_table_bias = table(abs_bias_out');


noabs_data_table_beta = table(noabs_beta_out');
noabs_data_table_bias = table(noabs_bias_out');


% % write output
writetable(abs_data_table_beta,strcat('Results/Excels/usual_suspects_fullpairwise_different_alphas',sorted,'_Norm_beta_Abs_deltax_', req_date,'.xlsx'))
writetable(abs_data_table_bias,strcat('Results/Excels/usual_suspects_fullpairwise_different_alphas',sorted,'_Norm_bias_Abs_deltax_', req_date,'.xlsx'))

writetable(noabs_data_table_beta,strcat('Results/Excels/usual_suspects_fullpairwise_different_alphas',sorted,'_Norm_beta_deltax_', req_date,'.xlsx'))
writetable(noabs_data_table_bias,strcat('Results/Excels/usual_suspects_fullpairwise_different_alphas',sorted,'_Norm_bias_deltax_', req_date,'.xlsx'))

%%%%%%%%%
%Uniform%
%%%%%%%%%

%Uniform
T_list = [10, 50, 100, 500, 1000];%[10, 50, 100, 500, 1000]
%reload and combine data
req_date ='29-Nov-2023';
abs_beta_out =[];
abs_bias_out =[];


noabs_beta_out =[];
noabs_bias_out =[];

for ss =1:length(T_list)
    
    abs_beta_add = load(strcat('Results/usual_suspects_fullpairwise_different_alphas',sorted,'_Unif_beta_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    abs_bias_add = load(strcat('Results/usual_suspects_fullpairwise_different_alphas',sorted,'_Unif_bias_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    
    abs_beta_out = [abs_beta_out,abs_beta_add.x'];
    abs_bias_out = [abs_bias_out,abs_bias_add.x'];
    
    noabs_beta_add = load(strcat('Results/usual_suspects_fullpairwise_different_alphas',sorted,'_Unif_beta_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    noabs_bias_add = load(strcat('Results/usual_suspects_fullpairwise_different_alphas',sorted,'_Unif_bias_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    
    noabs_beta_out = [abs_beta_out,abs_beta_add.x'];
    noabs_bias_out = [abs_bias_out,abs_bias_add.x'];

end


abs_data_table_beta = table(abs_beta_out');
abs_data_table_bias = table(abs_bias_out');


noabs_data_table_beta = table(noabs_beta_out');
noabs_data_table_bias = table(noabs_bias_out');


% % write output
writetable(abs_data_table_beta,strcat('Results/Excels/usual_suspects_fullpairwise_different_alphas',sorted,'_Unif_beta_Abs_deltax_', req_date,'.xlsx'))
writetable(abs_data_table_bias,strcat('Results/Excels/usual_suspects_fullpairwise_different_alphas',sorted,'_Unif_bias_Abs_deltax_', req_date,'.xlsx'))

writetable(noabs_data_table_beta,strcat('Results/Excels/usual_suspects_fullpairwise_different_alphas',sorted,'_Unif_beta_deltax_', req_date,'.xlsx'))
writetable(noabs_data_table_bias,strcat('Results/Excels/usual_suspects_fullpairwise_different_alphas',sorted,'_Unif_bias_deltax_', req_date,'.xlsx'))
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

%%%%%%%%
%Normal%
%%%%%%%%

sorted = '';
T_list = [10, 50, 100, 500, 1000];%[10, 50, 100, 500, 1000]

%reload and combine data
req_date ='29-Nov-2023';
abs_beta_out =[];
abs_bias_out =[];



for ss =1:length(T_list)
    
    abs_beta_add = load(strcat('Results/usual_suspects_nonsorted_adjacent_different_alphas_Norm_beta_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    abs_bias_add = load(strcat('Results/usual_suspects_nonsorted_adjacent_different_alphas_Norm_bias_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    
    abs_beta_out = [abs_beta_out,abs_beta_add.x'];
    abs_bias_out = [abs_bias_out,abs_bias_add.x'];

end


abs_data_table_beta = table(abs_beta_out');
abs_data_table_bias = table(abs_bias_out');



% % write output
writetable(abs_data_table_beta,strcat('Results/Excels/usual_suspects_nonsorted_adjacent_different_alphas_Norm_beta_Abs_deltax_', req_date,'.xlsx'))
writetable(abs_data_table_bias,strcat('Results/Excels/usual_suspects_nonsorted_adjacent_different_alphas_Norm_bias_Abs_deltax_', req_date,'.xlsx'))


%%%%%%%%%
%Uniform%
%%%%%%%%%

%Uniform
T_list = [10, 50, 100, 500, 1000];%[10, 50, 100, 500, 1000]
%reload and combine data
req_date ='29-Nov-2023';
abs_beta_out =[];
abs_bias_out =[];


for ss =1:length(T_list)
    
    abs_beta_add = load(strcat('Results/usual_suspects_nonsorted_adjacent_different_alphas_Unif_beta_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    abs_bias_add = load(strcat('Results/usual_suspects_nonsorted_adjacent_different_alphas_Unif_bias_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    
    abs_beta_out = [abs_beta_out,abs_beta_add.x'];
    abs_bias_out = [abs_bias_out,abs_bias_add.x'];
    
end


abs_data_table_beta = table(abs_beta_out');
abs_data_table_bias = table(abs_bias_out');



% % write output
writetable(abs_data_table_beta,strcat('Results/Excels/usual_suspects_nonsorted_adjacent_different_alphas_Unif_beta_Abs_deltax_', req_date,'.xlsx'))
writetable(abs_data_table_bias,strcat('Results/Excels/usual_suspects_nonsorted_adjacent_different_alphas_Unif_bias_Abs_deltax_', req_date,'.xlsx'))
