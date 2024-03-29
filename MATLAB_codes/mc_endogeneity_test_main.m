C:\Users\Kristof\Desktop\Endogeneity\MATLAB_codes%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%% version = 'C';
% case 1
abs_dummy = 'False'; 
dist = 'student-t';
sorted = 'True';
version = 'C';
run('ewpo_testing_fullpairwise')
% case 2
abs_dummy = 'False'; 
dist = 'student-t';
sorted = 'False';
version = 'C';
run('ewpo_testing_fullpairwise')
% case 3
abs_dummy = 'True'; 
dist = 'student-t';
sorted = 'False';
version = 'C';
run('ewpo_testing_fullpairwise')
% case 4
abs_dummy = 'True'; 
dist = 'student-t';
sorted = 'True';
version = 'C';
run('ewpo_testing_fullpairwise')


sorted = '';
%sorted = '_sorted';
%version = 'A';
%version = 'B';
version = 'C';


%%%%%%%%%
%Uniform%
%%%%%%%%%

%Uniform
T_list = [10, 50, 100, 500, 1000];
%reload and combine data
req_date ='08-Jan-2024';
abs_beta_out =[];
abs_diff_out =[];
abs_diff_beta_out =[];


noabs_beta_out =[];
noabs_diff_out =[];
noabs_diff_beta_out =[];

for ss =1:length(T_list)

    abs_beta_add = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Unif_beta_version_',version,'_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    abs_diff_add = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Unif_diff_matrix_version_',version,'_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    abs_diff_beta_add = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Unif_diff_beta_version_',version,'_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    
    abs_beta_out = [abs_beta_out,abs_beta_add.x'];
    abs_diff_out = [abs_diff_out,abs_diff_add.x'];
    abs_diff_beta_out = [abs_diff_beta_out,abs_diff_beta_add.x'];

   
    noabs_beta_add = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Unif_beta_version_',version,'_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    noabs_diff_add = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Unif_diff_matrix_version_',version,'_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    noabs_diff_beta_add = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Unif_diff_beta_version_',version,'_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));

    noabs_beta_out = [noabs_beta_out,noabs_beta_add.x'];
    noabs_diff_out = [noabs_diff_out,noabs_diff_add.x'];
    noabs_diff_beta_out = [noabs_diff_beta_out,noabs_diff_beta_add.x'];

end


abs_data_table_beta = table(abs_beta_out');
abs_data_table_diff = table(abs_diff_out');
abs_data_table_diff_beta = table(abs_diff_beta_out');


noabs_data_table_beta = table(noabs_beta_out');
noabs_data_table_diff = table(noabs_diff_out');
noabs_data_table_diff_beta = table(abs_diff_beta_out');



% % write output
writetable(abs_data_table_beta,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Unif_beta_version',version,'_Abs_deltax_', req_date,'.xlsx'))
writetable(abs_data_table_diff,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Unif_diff_version',version,'_Abs_deltax_', req_date,'.xlsx'))
writetable(abs_data_table_diff_beta,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Unif_diff_beta_version',version,'_Abs_deltax_', req_date,'.xlsx'))


writetable(noabs_data_table_beta,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Unif_beta_version',version,'_deltax_', req_date,'.xlsx'))
writetable(noabs_data_table_diff,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Unif_diff_version',version,'_deltax_', req_date,'.xlsx'))
writetable(noabs_data_table_diff_beta,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Unif_diff_beta_version',version,'_deltax_', req_date,'.xlsx'))


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% EwPO testing 2 estimators %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% case 1
abs_dummy = 'True'; 
dist = 'norm';
sorted = 'True';
run('ewpo_testing_2estimators_fullpairwise')
% case 2
abs_dummy = 'True'; 
dist = 'norm';
sorted = 'False';
run('ewpo_testing_2estimators_fullpairwise')
% case 3
abs_dummy = 'True'; 
dist = 'student-t';
sorted = 'True';
run('ewpo_testing_2estimators_fullpairwise')
% case 4
abs_dummy = 'True'; 
dist = 'student-t';
sorted = 'False';
run('ewpo_testing_2estimators_fullpairwise')

%%%%%%%%%%
% Normal %
%%%%%%%%%%

sorted = '';
%sorted = '_sorted';

%Normal
T_list = [10, 50, 100, 500, 1000];
%reload and combine data
req_date ='18-Jan-2024';
abs_beta_out =[];
abs_beta_corrected_out =[];
abs_beta_ols_out = [];
abs_beta_ewpo_out = [];
abs_beta_ols_corrected_out = [];
abs_beta_ewpo_corrected_out = [];
abs_ols_correction_out = [];
abs_ewpo_correction_out = [];


%noabs_beta_out =[];
%noabs_beta_corrected_out =[];
%noabs_beta_ols_out = [];
%noabs_beta_ewpo_out = [];
%noabs_beta_ols_corrected_out = [];
%noabs_beta_ewpo_corrected_out = [];
%noabs_ols_correction_out = [];
%noabs_ewpo_correction_out = [];

for ss =1:length(T_list)

    abs_beta_add = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Norm_beta_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    abs_beta_corrected_add = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Norm_beta_corrected_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    abs_beta_ols_add = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Norm_beta_ols_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    abs_beta_ewpo_add = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Norm_beta_ewpo_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    abs_beta_ols_corrected_add = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Norm_beta_ols_corrected_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    abs_beta_ewpo_corrected_add = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Norm_beta_ewpo_corrected_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    abs_ols_correction_add = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Norm_ols_correction_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    abs_ewpo_correction_add = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Norm_ewpo_correction_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    
    abs_beta_out = [abs_beta_out,abs_beta_add.x'];
    abs_beta_corrected_out = [abs_beta_corrected_out,abs_beta_corrected_add.x'];
    abs_beta_ols_out = [abs_beta_ols_out,abs_beta_ols_add.x'];
    abs_beta_ewpo_out = [abs_beta_ewpo_out,abs_beta_ewpo_add.x'];
    abs_beta_ols_corrected_out = [abs_beta_ols_corrected_out,abs_beta_ols_corrected_add.x'];
    abs_beta_ewpo_corrected_out = [abs_beta_ewpo_corrected_out,abs_beta_ewpo_corrected_add.x'];
    abs_ols_correction_out = [abs_ols_correction_out,abs_ols_correction_add.x'];
    abs_ewpo_correction_out = [abs_ewpo_correction_out,abs_ewpo_correction_add.x'];

    %{
    noabs_beta_add = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Norm_beta_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    noabs_beta_corrected_add = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Norm_beta_corrected_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    noabs_beta_ols_out = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Norm_beta_ols_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    noabs_beta_ewpo_out = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Norm_beta_ewpo_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    noabs_beta_ols_corrected_out = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Norm_beta_ols_corrected_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    noabs_beta_ewpo_corrected_out = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Norm_beta_ewpo_corrected_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    noabs_ols_correction_out = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Norm_ols_correction_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    noabs_ewpo_correction_out = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Norm_ewpo_correction_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
   
    
    noabs_beta_out = [noabs_beta_out,noabs_beta_add.x'];
    noabs_beta_corrected_out = [noabs_beta_corrected_out,noabs_beta_corrected_add.x'];
    noabs_beta_ols_out = [noabs_beta_ols_out,noabs_betaols_out_add.x'];
    noabs_beta_ewpo_out = [noabs_beta_ewpo_out,noabs_beta_ewpo_out_add.x'];
    noabs_beta_ols_corrected_out = [noabs_beta_ols_corrected_out,noabs_beta_ols_corrected_out_add.x'];
    noabs_beta_ewpo_corrected_out = [noabs_beta_ewpo_corrected_out,noabs_beta_ewpo_corrected_out_add.x'];
    noabs_ols_correction_out = [noabs_ols_correction_out,noabs_ols_correction_out_add.x'];
    noabs_ewpo_correction_out = [noabs_ewpo_correction_out,noabs_ewpo_correction_out_add.x'];
%}
end


abs_data_table_beta = table(abs_beta_out');
abs_data_table_beta_corrected = table(abs_beta_corrected_out');
abs_data_table_beta_ols = table(abs_beta_ols_out');
abs_data_table_beta_ewpo = table(abs_beta_ewpo_out');
abs_data_table_beta_ols_corrected = table(abs_beta_ols_corrected_out');
abs_data_table_beta_ewpo_corrected = table(abs_beta_ewpo_corrected_out');
abs_data_table_ols_correction = table(abs_ols_correction_out');
abs_data_table_ewpo_correction = table(abs_ewpo_correction_out');

%{
noabs_data_table_beta = table(noabs_beta_out');
noabs_data_table_beta_corrected = table(noabs_beta_corrected_out');
noabs_data_table_beta_ols = table(noabs_beta_ols_out');
noabs_data_table_beta_ewpo = table(noabs_beta_ewpo_out');
noabs_data_table_beta_ols_corrected = table(noabs_beta_ols_corrected_out');
noabs_data_table_beta_ewpo_corrected = table(noabs_beta_ewpo_corrected_out');
noabs_data_table_ols_correction = table(noabs_ols_correction_out');
noabs_data_table_ewpo_correction = table(noabs_ewpo_correction_out');
%}

% % write output
writetable(abs_data_table_beta,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Norm_beta_Abs_deltax_', req_date,'.xlsx'))
writetable(abs_data_table_beta_corrected,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Norm_beta_corrected_Abs_deltax_', req_date,'.xlsx'))
writetable(abs_data_table_beta_ols,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Norm_beta_ols_Abs_deltax_', req_date,'.xlsx'))
writetable(abs_data_table_beta_ewpo,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Norm_beta_ewpo_Abs_deltax_', req_date,'.xlsx'))
writetable(abs_data_table_beta_ols_corrected,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Norm_beta_ols_corrected_Abs_deltax_', req_date,'.xlsx'))
writetable(abs_data_table_beta_ewpo_corrected,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Norm_beta_ewpo_corrected_Abs_deltax_', req_date,'.xlsx'))
writetable(abs_data_table_ols_correction,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Norm_beta_ols_correction_Abs_deltax_', req_date,'.xlsx'))
writetable(abs_data_table_ewpo_correction,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Norm_beta_ewpo_correction_Abs_deltax_', req_date,'.xlsx'))

%{
writetable(noabs_data_table_beta,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Norm_beta_deltax_', req_date,'.xlsx'))
writetable(noabs_data_table_beta_corrected,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Norm_beta_corrected_deltax_', req_date,'.xlsx'))
writetable(noabs_data_table_beta_ols,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Norm_beta_ols_deltax_', req_date,'.xlsx'))
writetable(noabs_data_table_beta_ewpo,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Norm_beta_ewpo_deltax_', req_date,'.xlsx'))
writetable(noabs_data_table_beta_ols_corrected,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Norm_beta_ols_corrected_deltax_', req_date,'.xlsx'))
writetable(noabs_data_table_beta_ewpo_corrected,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Norm_beta_ewpo_corrected_deltax_', req_date,'.xlsx'))
writetable(noabs_data_table_ols_correction,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Norm_beta_ols_correction_deltax_', req_date,'.xlsx'))
writetable(noabs_data_table_ewpo_correction,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Norm_beta_ewpo_correction_deltax_', req_date,'.xlsx'))
%}

%%%%%%%%%%%%%
% Student-t %
%%%%%%%%%%%%%

sorted = '_sorted';

%Student-t
T_list = [10, 50, 100, 500, 1000];
%reload and combine data
req_date ='18-Jan-2024';
abs_beta_out =[];
abs_beta_corrected_out =[];
abs_beta_ols_out = [];
abs_beta_ewpo_out = [];
abs_beta_ols_corrected_out = [];
abs_beta_ewpo_corrected_out = [];
abs_ols_correction_out = [];
abs_ewpo_correction_out = [];


%noabs_beta_out =[];
%noabs_beta_corrected_out =[];
%noabs_beta_ols_out = [];
%noabs_beta_ewpo_out = [];
%noabs_beta_ols_corrected_out = [];
%noabs_beta_ewpo_corrected_out = [];
%noabs_ols_correction_out = [];
%noabs_ewpo_correction_out = [];

for ss =1:length(T_list)

    abs_beta_add = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Student-t_beta_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    abs_beta_corrected_add = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Student-t_beta_corrected_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    abs_beta_ols_add = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Student-t_beta_ols_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    abs_beta_ewpo_add = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Student-t_beta_ewpo_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    abs_beta_ols_corrected_add = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Student-t_beta_ols_corrected_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    abs_beta_ewpo_corrected_add = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Student-t_beta_ewpo_corrected_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    abs_ols_correction_add = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Student-t_ols_correction_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    abs_ewpo_correction_add = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Student-t_ewpo_correction_Abs_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    
    abs_beta_out = [abs_beta_out,abs_beta_add.x'];
    abs_beta_corrected_out = [abs_beta_corrected_out,abs_beta_corrected_add.x'];
    abs_beta_ols_out = [abs_beta_ols_out,abs_beta_ols_add.x'];
    abs_beta_ewpo_out = [abs_beta_ewpo_out,abs_beta_ewpo_add.x'];
    abs_beta_ols_corrected_out = [abs_beta_ols_corrected_out,abs_beta_ols_corrected_add.x'];
    abs_beta_ewpo_corrected_out = [abs_beta_ewpo_corrected_out,abs_beta_ewpo_corrected_add.x'];
    abs_ols_correction_out = [abs_ols_correction_out,abs_ols_correction_add.x'];
    abs_ewpo_correction_out = [abs_ewpo_correction_out,abs_ewpo_correction_add.x'];

    %{
    noabs_beta_add = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Student-t_beta_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    noabs_beta_corrected_add = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Student-t_beta_corrected_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    noabs_beta_ols_out = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Student-t_beta_ols_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    noabs_beta_ewpo_out = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Student-t_beta_ewpo_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    noabs_beta_ols_corrected_out = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Student-t_beta_ols_corrected_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    noabs_beta_ewpo_corrected_out = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Student-t_beta_ewpo_corrected_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    noabs_ols_correction_out = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Student-t_ols_correction_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
    noabs_ewpo_correction_out = load(strcat('Results/ewpo_testing_fullpairwise',sorted,'_Student-t_ewpo_correction_deltax_T_',num2str(T_list(ss)),'_',req_date,'.mat'));
   
    
    noabs_beta_out = [noabs_beta_out,noabs_beta_add.x'];
    noabs_beta_corrected_out = [noabs_beta_corrected_out,noabs_beta_corrected_add.x'];
    noabs_beta_ols_out = [noabs_beta_ols_out,noabs_betaols_out_add.x'];
    noabs_beta_ewpo_out = [noabs_beta_ewpo_out,noabs_beta_ewpo_out_add.x'];
    noabs_beta_ols_corrected_out = [noabs_beta_ols_corrected_out,noabs_beta_ols_corrected_out_add.x'];
    noabs_beta_ewpo_corrected_out = [noabs_beta_ewpo_corrected_out,noabs_beta_ewpo_corrected_out_add.x'];
    noabs_ols_correction_out = [noabs_ols_correction_out,noabs_ols_correction_out_add.x'];
    noabs_ewpo_correction_out = [noabs_ewpo_correction_out,noabs_ewpo_correction_out_add.x'];
%}
end


abs_data_table_beta = table(abs_beta_out');
abs_data_table_beta_corrected = table(abs_beta_corrected_out');
abs_data_table_beta_ols = table(abs_beta_ols_out');
abs_data_table_beta_ewpo = table(abs_beta_ewpo_out');
abs_data_table_beta_ols_corrected = table(abs_beta_ols_corrected_out');
abs_data_table_beta_ewpo_corrected = table(abs_beta_ewpo_corrected_out');
abs_data_table_ols_correction = table(abs_ols_correction_out');
abs_data_table_ewpo_correction = table(abs_ewpo_correction_out');

%{
noabs_data_table_beta = table(noabs_beta_out');
noabs_data_table_beta_corrected = table(noabs_beta_corrected_out');
noabs_data_table_beta_ols = table(noabs_beta_ols_out');
noabs_data_table_beta_ewpo = table(noabs_beta_ewpo_out');
noabs_data_table_beta_ols_corrected = table(noabs_beta_ols_corrected_out');
noabs_data_table_beta_ewpo_corrected = table(noabs_beta_ewpo_corrected_out');
noabs_data_table_ols_correction = table(noabs_ols_correction_out');
noabs_data_table_ewpo_correction = table(noabs_ewpo_correction_out');
%}

% % write output
writetable(abs_data_table_beta,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Student-t_beta_Abs_deltax_', req_date,'.xlsx'))
writetable(abs_data_table_beta_corrected,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Student-t_beta_corrected_Abs_deltax_', req_date,'.xlsx'))
writetable(abs_data_table_beta_ols,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Student-t_beta_ols_Abs_deltax_', req_date,'.xlsx'))
writetable(abs_data_table_beta_ewpo,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Student-t_beta_ewpo_Abs_deltax_', req_date,'.xlsx'))
writetable(abs_data_table_beta_ols_corrected,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Student-t_beta_ols_corrected_Abs_deltax_', req_date,'.xlsx'))
writetable(abs_data_table_beta_ewpo_corrected,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Student-t_beta_ewpo_corrected_Abs_deltax_', req_date,'.xlsx'))
writetable(abs_data_table_ols_correction,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Student-t_beta_ols_correction_Abs_deltax_', req_date,'.xlsx'))
writetable(abs_data_table_ewpo_correction,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Student-t_beta_ewpo_correction_Abs_deltax_', req_date,'.xlsx'))

%{
writetable(noabs_data_table_beta,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Student-t_beta_deltax_', req_date,'.xlsx'))
writetable(noabs_data_table_beta_corrected,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Student-t_beta_corrected_deltax_', req_date,'.xlsx'))
writetable(noabs_data_table_beta_ols,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Student-t_beta_ols_deltax_', req_date,'.xlsx'))
writetable(noabs_data_table_beta_ewpo,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Student-t_beta_ewpo_deltax_', req_date,'.xlsx'))
writetable(noabs_data_table_beta_ols_corrected,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Student-t_beta_ols_corrected_deltax_', req_date,'.xlsx'))
writetable(noabs_data_table_beta_ewpo_corrected,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Student-t_beta_ewpo_corrected_deltax_', req_date,'.xlsx'))
writetable(noabs_data_table_ols_correction,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Student-t_beta_ols_correction_deltax_', req_date,'.xlsx'))
writetable(noabs_data_table_ewpo_correction,strcat('Results/Excels/ewpo_testing_fullpairwise',sorted,'_Student-t_beta_ewpo_correction_deltax_', req_date,'.xlsx'))
%}