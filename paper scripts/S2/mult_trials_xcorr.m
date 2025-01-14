%% Figure 2 plotting 
% Seph Marshall, Nov 2020
%can combine multiple  independant trial data, plots 2 graphs: individual
%cell traces w average bolded, and trial average w 95 conf intervals 

%% Combined Xcorrelation, 2 independant trials w SEM intervals 
 % data initialization
 
clear; clc; 

close all; 

root = 'C:\Users\marsh\OneDrive - McGill University\research paper\results good_Feb2023\Fig 2\RhoA2G\compiled Xcorr depth 1.95 um'; 

load([root, filesep,'Edge Vel. RhoA2G T1 depth 6.mat']); 

averages_T1 = cell_arr; 

for i=1:size(averages_T1,2)
    
    plotting_table(i,:) =averages_T1{1,i};
 
end 

cell_arr =[]; 

load([root, filesep,'Edge Vel. RhoA2G T2 depth 6.mat'])

averages_T2 = cell_arr; 

for d = 1:size(averages_T2,2)
    
      plotting_table(d+ size(averages_T1,2),:) = averages_T2{1,d};
   
end 

cell_arr =[]; 

load([root, filesep,'Edge Vel. RhoA2G T4 depth 6.mat'])

averages_T4= cell_arr; 

for d = 1:size(averages_T4,2)
    
      plotting_table(d+ size(averages_T1,2)+size(averages_T2,2),:) = averages_T4{1,d};
   
end 
stats_arr = zeros(3,41); 

 for k = 1:41 % -20:20 lags, respectively 
     
     % fitting normal dist to each lag 
    pd = fitdist(plotting_table(:,k),'Normal'); 
    stats_arr(1,k) = pd.mu; 
    ci = paramci(pd); 
    
    %95 % confidence interval, upper and lower 
    stats_arr(2,k) = ci(1,1); 
    stats_arr(3,k)=ci(2,1); 
    stats_arr(4,k) = pd.sigma; 
 end 
 
%% plotting 


f1 = figure; 

hold on;
ylim([-0.7 0.4]);
xline(0,'--'); 
yline(0,'--'); 
for x= 1:size(plotting_table,1)
    
    plot(-20:20,plotting_table(x,:));
    
end 

plot(-20:20,(stats_arr(1,:)),'Color',[0,0,0], 'LineWidth', 3 ); 

xticks([-19.2 -16.8 -14.4 -12.0 -9.6 -7.2 -4.8 -2.4 0 2.4 4.8 7.2 9.6 12.0 14.4 16.8 19.2])
xticklabels({'-8' '' '-6' '' '-4' '' '-2' '' '0' '' '2' '' '4' '' '6' '' '8'}); 

 hold off; 
 
 f2= figure; 
 
 hold on; 
 ylim([-0.3 0.6]);
xline(0,'--'); 
yline(0,'--'); 

plot(-20:20,(stats_arr(1,:)),'Color',[1,0,0], 'LineWidth', 3 ); 
plot(-20:20,(stats_arr(2,:)),'Color',[1,0,0] );
plot(-20:20,(stats_arr(3,:)),'Color',[1,0,0] );

xticks([-19.2 -16.8 -14.4 -12.0 -9.6 -7.2 -4.8 -2.4 0 2.4 4.8 7.2 9.6 12.0 14.4 16.8 19.2])
xticklabels({'-8' '' '-6' '' '-4' '' '-2' '' '0' '' '2' '' '4' '' '6' '' '8'}); 





 saveas(f1,[root,filesep,'Edge vel RhoA2G_xcorr_n=40_depth6.svg']); 
% 
% saveas(f2,[root,filesep,'Rho_avg+SEM_n=27_depth6.svg']); 

