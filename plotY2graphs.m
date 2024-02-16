
clc; clear; 
%% cntrl loading 

 load("F:\230914 - Y2 motility 1\cropped\cntrl_LIS\stats_25px_thresh.mat"); 
 cntrl_1_25 = stats; 
 
%  load("F:\230914 - Y2 motility 1\cropped\cntrl_LIS\stats_50px_thresh.mat"); 
%   cntrl_1_50 = stats;
  
  
  load("F:\230919 - Y2motility 2\cropped\cntrl_LIS\stats_25px_thresh.mat"); 
 cntrl_2_25 = stats; 
 
%  load("F:\230919 - Y2motility 2\cropped\cntrl_LIS\stats_50px_thresh.mat"); 
%   cntrl_2_50 = stats;
%   
  
  
  
  
%% Y2 loading   
 load("F:\230914 - Y2 motility 1\cropped\20uM_Y2\stats_25px_thresh.mat"); 
 Y2_1_25 = stats; 
 
%   load("F:\230914 - Y2 motility 1\cropped\20uM_Y2\stats_50px_thresh.mat"); 
%  Y2_1_50 = stats; 
%  
  load("F:\230919 - Y2motility 2\cropped\20uM_Y2\stats_25px_thresh.mat"); 
 Y2_2_25 = stats; 
 
%   load("F:\230919 - Y2motility 2\cropped\20uM_Y2\stats_50px_thresh.mat"); 
%  Y2_2_50 = stats; 
 
 
 
 %% plotting 
 
 f = figure; 
 hold on; 
 
 title('Av protrusion speed'); 
 ylabel('um/min'); 
xlim([0,4]); 
ylim([6.6 7.6]); 
xticks([1,3]); 
xticklabels({'Control', '20uM Y27632'}); 
x = [1,3]; 
 trial_1 = [ nanmean(cntrl_1_25.prot_speed), nanmean(Y2_1_25.prot_speed)]; 
 trial_2 = [ nanmean(cntrl_2_25.prot_speed), nanmean(Y2_2_25.prot_speed)]; 
 
 %plot(x, trial_1, 'b-o', x, trial_2, 'r-o'); 
 
 hold off; 

 
 %% 
 
 f2 = figure; 
 ylim([4 11]); 
 struct.Cntrl = [abs(cntrl_1_25.prot_speed), abs(cntrl_2_25.prot_speed)]; 
 struct.Y = [abs(Y2_1_25.prot_speed), abs(Y2_2_25.prot_speed)]; 
 
 violinplot(struct); 
 
 p = ranksum(struct.Cntrl,struct.Y);  
 
 
 %% saving 
 root = 'C:\Users\marsh\OneDrive - McGill University\research paper\results good_Feb2023\Fig 4\control vs Y2 motility graphs'; 
 saveas(f,[root, filesep, 'protrusion_speed_25_px.fig']); 
 
 