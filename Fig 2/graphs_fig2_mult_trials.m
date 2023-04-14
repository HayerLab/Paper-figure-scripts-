% Fig 2 graphs triple replicates 

%%
clc; 
depth =6;

root1 = 'F:\Seph\data\data_221006 - Trial 4 RhoA, RhoA2G cyto\RhoA2G\cropped\graphs depth_6\overall averages Rho\statistics.mat'; 
root2 = 'F:\Seph\data\data_210331 - Trial 2 RhoA, RhoA2G cyto correct\RHOA2G\cropped\graphs depth 6\overall averages Rho\statistics.mat';  
%root3= 'F:\Seph\data\data_200116 - Trial 3 Rho, Myosin\cropped\graphs depth 6\overall averages Rho\statistics.mat'; 

datadir = 'F:\Seph\RhoA2G fig 2 data'; 

load(root1); 
prot_arr = protrusionsfinal; 
retr_arr = retractionsfinal; 

mean_prot = stats2.mu; 
mean_retr = stats.mu; 

load(root2); 
prot_arr = [prot_arr, protrusionsfinal]; 
retr_arr = [retr_arr, retractionsfinal]; 

mean_prot = [mean_prot, stats2.mu]; 
mean_retr = [mean_retr, stats.mu]; 

%  load(root3); 
%  prot_arr = [prot_arr, protrusionsfinal]; 
%  retr_arr = [retr_arr, retractionsfinal]; 
% % 
% mean_prot = [mean_prot, stats2.mu]; 
% mean_retr = [mean_retr, stats.mu]; 

%%

f5=figure; 
 

 hold on; 
 ylim([0.6 1.4]); 
 

 
 [line1, xi]=ksdensity(prot_arr,'kernel','normpdf','npoints',size(prot_arr,2)); 
 [line2, ai]=ksdensity(retr_arr,'kernel','normpdf','npoints',size(retr_arr,2)); 
 
 plot(line1,xi); 
 plot(line2,ai); 
  

%   

  saveas(f5,strcat(datadir, filesep,'RhoA2G_3trials_depth',num2str(depth),'.svg') ); 

 
 
%   histogram(protrusionsfinal,'Normalization','probability', 'DisplayStyle', 'stairs');
%   histogram(protrusionsfinal,'Normalization','probability', 'DisplayStyle', 'stairs');
% 
%% 


f4=figure; 
   
   hold on 
    
    xlim([0 3])
   yline(1,'--'); 
  % ytitle(Normalized FRET Intensity); 
  xticks([2 4])
  xticklabels({'Protrusions', 'Retractions'}); 
  
  C = [prot_arr retr_arr]; 
grp = [zeros(1,size(prot_arr,2)),ones(1,size(retr_arr,2))];
  boxplot(C,grp,'labels',{'Protrusions', 'Retractions'}); 
  title(strcat('RhoA2G- Depth',num2str(depth))); 
  
  ylim([0.6 1.4]);
    

   
  % errorbar(2,mean(protrusionsfinal),std(protrusionsfinal),'-sk','CapSize',15,'LineWidth',3); 
 %  errorbar(4,mean(retractionsfinal), std(retractionsfinal),'-sk', 'CapSize',15,'LineWidth',3); 
   
  [p_val]= ranksum(prot_arr, retr_arr) % 'Vartype', 'unequal'); %look at the ttest2 matlab documentation, chose option with enqual variance her
   
  group = {[1,2]};
  %sigstar(group, p_val)
   
   hold off; 
   
    saveas(f4,strcat(datadir,'\graphs depth_',num2str(depth),'error bar.fig'));
    saveas(f4,strcat(datadir,'\graphs depth_',num2str(depth),'error bar sig.png'));
% 
    save(strcat(datadir,'\graphs depth_',num2str(depth),'statistics.mat'),'mean_prot','mean_retr','p_val','prot_arr','retr_arr'); 
%    
   %% 
    
f6= figure;    
   hold on; 
   
   xlim([0,3]); 
   xticks([0:1:3]); 
   ylim([0.7 1.3]); 
   
   yline(1,'--'); 
   
   for i = 1:2
       
 x = [1 2]; 
 y = [mean_prot(1,i), mean_retr(1,i)]; 
 
 plot(x,y,'-+');
   end 
  
    saveas(f6,strcat(datadir,'\graphs depth_',num2str(depth),'trial_replicates.fig'));
