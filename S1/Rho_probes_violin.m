% script to violin plot R score distributions 
% SM Dec 20th 2022 

%% RhoA2G section
clc; clear; 

root1 = 'F:\Seph\data\data_210331 - Trial 2 RhoA, RhoA2G cyto correct\RHOA2G\cyto bias correlation-YFP'; 
root2 = 'F:\Seph\data\data_221006 - Trial 4 RhoA, RhoA2G cyto\RhoA2G\cyto bias correlation-YFP'; 

 load([root1, filesep, 'R values.mat']); 
 R_arr_A2G = R_arr'; 
 
load([root2, filesep, 'R values.mat']); 
R_arr_A2G = [R_arr_A2G, R_arr']; 
%% RhoA 
root1 = 'F:\Seph\data\data_210331 - Trial 2 RhoA, RhoA2G cyto correct\RHOA\cyto bias correlation-YFP'; 
root2 = 'F:\Seph\data\data_221006 - Trial 4 RhoA, RhoA2G cyto\RhoA\cyto bias correlation-YFP'; 

 load([root1, filesep, 'R values.mat']); 
 R_arr_A = R_arr'; 
 
load([root2, filesep, 'R values.mat']); 
R_arr_A = [R_arr_A, R_arr']; 

%% RhoB, 
root1 = 'D:\221207_40x 2x2 bin_RhoB_ezrin_cyto\cyto bias correlation-YFP'; 
root2 = 'D:\221209 - 40x 2x2 bin_RhoB_cyto\cyto bias correlation-YFP-new'; 
root3= 'D:\221219_40x_2x2bin_RhoB_cyto\cyto bias correlation-YFP'; 
 
load([root1, filesep, 'R values.mat']); 
 R_arr_B = R_arr'; 
 
load([root2, filesep, 'R values.mat']); 
R_arr_B = [R_arr_B, R_arr']; 

load([root3, filesep, 'R values.mat']); 
R_arr_B = [R_arr_B, R_arr']; 

%% violin plot that shit 

s = struct; 
s.RhoA2G = R_arr_A2G; 
s.RhoA = R_arr_A; 
s.RhoB = R_arr_B; 

f= figure; 

grid on; 
violinplot(s); 

yticks(-1:0.25:0.5); 


 y = [ R_arr_A2G,R_arr_A, R_arr_B]; 
 group = repelem(1:3, 1, [numel(R_arr_A2G),numel(R_arr_A),numel(R_arr_B)]);

 %%
 [p, anova, stats] = anova1(y, group); 
 
 [compare]= multcompare(stats) 


save(['C:\Users\gmarsh8\OneDrive - McGill University\research paper\poster presentation\Rho probe comparison figures', filesep, 'R values all 3 probes.mat']);  


