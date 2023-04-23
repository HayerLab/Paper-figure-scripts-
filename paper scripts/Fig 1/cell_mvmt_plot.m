clear; clc; 

f = figure; 

hold on; 

y = [6.39, -5.93]; 

error = [0.61, 0.63]; 

x = [1,2]; 
xticks([1 2]); 


xticklabels(['protrusions'; 'retractions'])


errorbar(x,y,error, 'o', 'LineStyle', 'none' )

xlim([0 3]); 
ylim([-8  8]);
yline(0, '--'); 

%% violin plots

load(['C:\Users\gmarsh8\OneDrive - McGill University\research paper\results good_Feb2023\Fig 1\cell mvmt statistics\compiled_stats.mat']); 

S = struct; 
% number_protrusions=number_protrusions; %./60; 
% number_retractions=number_retractions; %./60; 
%retr_speed=abs(retr_speed); 
% prot_perimeter = (prot_perimeter)*100; 
% retr_perimeter = (retr_perimeter)*100; 

S.prot_size = prot_size; 
S.retr_size = retr_size; 

f= figure; 
grid on; 
%ylim([2 15]); 
ylabel(['event size pixels']); 
violinplot(S); 