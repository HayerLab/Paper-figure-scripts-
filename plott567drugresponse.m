
clc; clear; 
root = 'F:\231113_RhoB_ezrin_ WTvsT567_10umNOC\cropped\t567'; 
datadir = 'F:\231113_RhoB_ezrin_ WTvsT567_10umNOC\activation_comparison'; 

t567_control = []; 
t567_NOC = []; 
for i =1:20
   datadir_cntrl = [root, filesep,  filesep,'control', filesep, num2str(i), filesep, 'output\ezrin_data']; 
  
 load([datadir_cntrl, filesep, 'CytoRatioData']);
 t567_control(i,:) = ezrin_mean_norm(1,1:35); 
     
end 


for i =1:17
    
 if i ==11 
     continue; 
 end 
   datadir_NOC = [root, filesep,  filesep,'+10umNOC', filesep, num2str(i), filesep, 'output\ezrin_data']; 
  
 load([datadir_NOC, filesep, 'CytoRatioData']);
 t567_NOC(i,:) = ezrin_mean_norm(1,1:35); 
     
end 
t567_NOC(11, :) = []; 

%get statistics on the distributions of the data 
for k = 1:35 % -20:20 lags, respectively 
     
     % fitting normal dist to each lag 
    pd = fitdist(t567_control(:,k),'Normal'); 
    stats_arr_cntrl_t567(1,k) = pd.mu; 
    ci = paramci(pd); 
    
    %95 % confidence interval, upper and lower 
    stats_arr_cntrl_t567(2,k) = ci(1,1); 
    stats_arr_cntrl_t567(3,k)=ci(2,1); 
    stats_arr_cntrl_t567(4,k) = pd.sigma; 
    
     pd2 = fitdist(t567_NOC(:,k),'Normal'); 
    stats_arr_NOC_t567(1,k) = pd2.mu; 
    ci = paramci(pd2); 
    
    %95 % confidence interval, upper and lower 
    stats_arr_NOC_t567(2,k) = ci(1,1); 
    stats_arr_NOC_t567(3,k)=ci(2,1); 
    stats_arr_NOC_t567(4,k) = pd2.sigma; 
end 
 
f1= figure; 

hold on; 
ylim([0.8 1.3]); 

plot(1:35,(stats_arr_cntrl_t567(1,:)),'Color',[1,0,0], 'LineWidth', 3 ); 
plot(1:35,(stats_arr_cntrl_t567(2,:)),'Color',[1,0,0] );
plot(1:35,(stats_arr_cntrl_t567(3,:)),'Color',[1,0,0] );

plot(1:35,(stats_arr_NOC_t567(1,:)),'Color',[0,0,1], 'LineWidth', 3 ); 
plot(1:35,(stats_arr_NOC_t567(2,:)),'Color',[0,0,1] );
plot(1:35,(stats_arr_NOC_t567(3,:)),'Color',[0,0,1] );

save([datadir, filesep, 't567data.mat'], 't567_control', 't567_NOC', 'stats_arr_cntrl', 'stats_arr_NOC'); 
saveas(f1,[datadir,filesep,'t567 cntrl vs NOC ezrin ratio signal.fig']);

%% 
load('F:\231113_RhoB_ezrin_ WTvsT567_10umNOC\activation_comparison\WTdata.mat'); 

ff= figure; 

hold on; 
ylim([0.8 1.3]); 

plot(1:35,(stats_arr_cntrl_t567(1,:)),'Color',[1,0,0], 'LineWidth', 3 ); 
plot(1:35,(stats_arr_cntrl_t567(2,:)),'Color',[1,0,0] );
plot(1:35,(stats_arr_cntrl_t567(3,:)),'Color',[1,0,0] );

plot(1:35,(stats_arr_NOC_t567(1,:)),'Color',[0,0,1], 'LineWidth', 3 ); 
plot(1:35,(stats_arr_NOC_t567(2,:)),'Color',[0,0,1] );
plot(1:35,(stats_arr_NOC_t567(3,:)),'Color',[0,0,1] );

plot(1:35,(stats_arr_cntrl(1,:)),'Color',[0,1,0], 'LineWidth', 3 ); 
plot(1:35,(stats_arr_cntrl(2,:)),'Color',[0,1,0] );
plot(1:35,(stats_arr_cntrl(3,:)),'Color',[0,1,0] );

plot(1:35,(stats_arr_NOC(1,:)),'Color',[0,0,0], 'LineWidth', 3 ); 
plot(1:35,(stats_arr_NOC(2,:)),'Color',[0,0,0] );
plot(1:35,(stats_arr_NOC(3,:)),'Color',[0,0,0] );

