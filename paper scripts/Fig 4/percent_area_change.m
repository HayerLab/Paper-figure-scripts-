%% percent area change 
% script to calculate changes in cell size pre/post drug treatment 
clc; clear; 
%root1 = 'F:\Seph\data\data_221104 - Trial 5 RhoB Myosin drug treatments\data'; 
 %root1  ='F:\Seph\data\data_210312 - RhoB drug treatments 20x\data-redone'; 
 root1 = 'F:\Seph\data\data_210317 - Trial 2 RhoB drug treatments 20x\data-redone'; 
 
 %% 
 
 
 change_control = []; 
 change_Y2 = []; 
 change_ML7 = []; 
 change_Y2_ML7=[]; 
 
 k=0;
for row=1:2
%     
    for col=1
        for site=1:4


            k=k+1;
            position_cntrl{k}=[num2str(row),'_',num2str(col),'_',num2str(site)];
        end
    end
end

  k=0;
for row=1:2
%     
    for col=1
        for site=5:8


            k=k+1;
            position_Y2{k}=[num2str(row),'_',num2str(col),'_',num2str(site)];
        end
    end
end

 k=0;
for row=1:2
%     
    for col=1
        for site=9:12

            k=k+1;
            position_ML7{k}=[num2str(row),'_',num2str(col),'_',num2str(site)];
        end
    end
end
 
 k=0;
for row=1:2
%     
    for col=1
        for site=13:16

            k=k+1;
            position_Y2_ML7{k}=[num2str(row),'_',num2str(col),'_',num2str(site)];
        end
    end
end
 
%% control
percent_change_cntrl = []; 

for i = 1: length(position_cntrl)
    load([root1, filesep,strcat(position_cntrl{i}, '_RatioData_raw.mat')]); 
    mask_pre = []; 
    mask_post = []; 
   mask_pre = maskFinal{1,5}(maskFinal{1,5} ==1); 
   mask_post = maskFinal{1,35}(maskFinal{1,35} ==1); 
   
   percent_change = ((size(mask_post,1)-size(mask_pre,1))/size(mask_pre,1))*100; 
   
   percent_change_cntrl = [percent_change_cntrl,percent_change]; 

end 
 
%% 
percent_change_Y2 = []; 


for i = 1: length(position_Y2)
    mask_pre =[]; 
    mask_post =[]; 
    load([root1, filesep,strcat(position_Y2{i}, '_RatioData_raw.mat')]); 
    
   mask_pre = maskFinal{1,5}(maskFinal{1,5} ==1); 
   mask_post = maskFinal{1,35}(maskFinal{1,35} ==1); 
   
   percent_change = ((size(mask_post,1)-size(mask_pre,1))/size(mask_pre,1))*100; 
   
   percent_change_Y2 = [percent_change_Y2,percent_change]; 

end 

%% 
percent_change_ML7 = []; 


for i = 1: length(position_ML7)
    mask_pre =[]; 
    mask_post =[]; 
    load([root1, filesep,strcat(position_ML7{i}, '_RatioData_raw.mat')]); 
    
   mask_pre = maskFinal{1,5}(maskFinal{1,5} ==1); 
   mask_post = maskFinal{1,35}(maskFinal{1,35} ==1); 
   
   percent_change = ((size(mask_post,1)-size(mask_pre,1))/size(mask_pre,1))*100; 
   
   percent_change_ML7 = [percent_change_ML7,percent_change]; 

end 
 
%% 
percent_change_Y2_ML7 = []; 


for i = 1: length(position_Y2_ML7)
    mask_pre =[]; 
    mask_post =[]; 
    load([root1, filesep,strcat(position_Y2_ML7{i}, '_RatioData_raw.mat')]); 
    
   mask_pre = maskFinal{1,5}(maskFinal{1,5} ==1); 
   mask_post = maskFinal{1,35}(maskFinal{1,35} ==1); 
   
   percent_change = ((size(mask_post,1)-size(mask_pre,1))/size(mask_pre,1))*100; 
   
   percent_change_Y2_ML7 = [percent_change_Y2_ML7,percent_change]; 

end 

%% 
 save('C:\Users\gmarsh8\OneDrive - McGill University\research paper\results good_Feb2023\Fig 4\area_change_statistics\trial3.mat', 'percent_change_cntrl','percent_change_Y2', 'percent_change_ML7', 'percent_change_Y2_ML7');  
 
 
 %% 
 clear; 
 cntrl =[];
 Y2 =[]; 
 ML7 =[]; 
 Y2_ML7 = []; 
 
 load('C:\Users\gmarsh8\OneDrive - McGill University\research paper\results good_Feb2023\Fig 4\area_change_statistics\trial1.mat'); 
 
 cntrl = [cntrl, percent_change_cntrl]; 
 Y2 =[Y2, percent_change_Y2]; 
 ML7 =[ML7, percent_change_ML7]; 
 Y2_ML7 = [Y2_ML7, percent_change_Y2_ML7];
 
 load('C:\Users\gmarsh8\OneDrive - McGill University\research paper\results good_Feb2023\Fig 4\area_change_statistics\trial2.mat'); 
 
 cntrl = [cntrl, percent_change_cntrl]; 
 Y2 =[Y2, percent_change_Y2]; 
 ML7 =[ML7, percent_change_ML7]; 
 Y2_ML7 = [Y2_ML7, percent_change_Y2_ML7];
 
  
 load('C:\Users\gmarsh8\OneDrive - McGill University\research paper\results good_Feb2023\Fig 4\area_change_statistics\trial3.mat'); 
 
 cntrl = [cntrl, percent_change_cntrl]; 
 Y2 =[Y2, percent_change_Y2]; 
 ML7 =[ML7, percent_change_ML7]; 
 Y2_ML7 = [Y2_ML7, percent_change_Y2_ML7];
 
 %% 
 s = struct; 



s.Y2_ML7 = Y2_ML7; 
s.ML7 = ML7; 
s.DY2 = Y2; 
s.CNTRL = cntrl; 

f= figure; 

grid on; 
violinplot(s); 

yticks(-100:10:40); 
ylim([-100, 40]); 

% 
  y = [ cntrl, Y2, ML7, Y2_ML7]; 
  group = repelem(1:4, 1, [numel(cntrl),numel(Y2),numel(ML7), numel(Y2_ML7)]);
 
  %% 
 [p, anova, stats] = anova1(y, group); 
 
 [compare]= multcompare(stats) 

 save('C:\Users\gmarsh8\OneDrive - McGill University\research paper\results good_Feb2023\Fig 4\area_change_statistics\compiled.mat', 'compare', 'Y2', 'ML7', 'cntrl', 'Y2_ML7'); 
 
 
 
 
 
 