% percent area change nocodazloe nada data 

% for i = 1:80
%     mask_temp=mask(:,:,i); 
%     mask_temp =mask_temp(mask_temp==1); 
%    area(1,i) = size(mask_temp,1); 
%    
% end 
% 
% plot(1:80, area); 

clc; 
root = 'F:\nada nocodazole data\trial'; 

cntrl= [4,5,6,10,11,12,16,17,18]; 
NOC = [1,2,3,7,8,9,13,14,15]; 
trial1_date ='07-04-23-'; 
trial2_date ='07-05-23-'; 
trial3_date='07-06-23-'; 
trial4_date = '07-11-23-'; 

%% control section

percent_change_cntrl = []; 

for trial = 1:4 

for i = 1:9 
        mask_temp_pre = []; 
    mask_temp_post = []; 
    if trial ==1 
    if cntrl(1,i) <10
    load([root,num2str(trial), filesep, strcat(trial1_date, '0', num2str(cntrl(1,i))),'-01-Mask.mat']); 
    
    else 
    load([root,num2str(trial), filesep, strcat(trial1_date, num2str(cntrl(1,i))),'-01-Mask.mat']); 
          
    end 
    end 
     if trial ==2 
    if cntrl(1,i) <10
    load([root,num2str(trial), filesep, strcat(trial2_date, '0', num2str(cntrl(1,i))),'-01-Mask.mat']); 
    
    else 
    load([root,num2str(trial), filesep, strcat(trial2_date, num2str(cntrl(1,i))),'-01-Mask.mat']); 
          
    end 
     end 
     if trial ==3
    if cntrl(1,i) <10
    load([root,num2str(trial), filesep, strcat(trial3_date, '0', num2str(cntrl(1,i))),'-01-Mask.mat']); 
    
    else 
    load([root,num2str(trial), filesep, strcat(trial3_date, num2str(cntrl(1,i))),'-01-Mask.mat']); 
          
    end 
     end 
    
     if trial ==4
    if cntrl(1,i) <10
    load([root,num2str(trial), filesep, strcat(trial4_date, '0', num2str(cntrl(1,i))),'-01-Mask.mat']); 
    
    else 
    load([root,num2str(trial), filesep, strcat(trial4_date, num2str(cntrl(1,i))),'-01-Mask.mat']); 
          
    end 
    end 
     mask_temp_pre=mask(:,:,10); 
    mask_temp_pre =mask_temp_pre(mask_temp_pre==1); 
     mask_temp_post=mask(:,:,40); 
    mask_temp_post =mask_temp_post(mask_temp_post==1); 
    
    
%     mask_pre = mask(:,:,10)(mask(:,:,10) ==1); 
%    mask_post = mask(:,:,40)(mask(:,:,40) ==1); 
   
   percent_change = ((size(mask_temp_post,1)-size(mask_temp_pre,1))/size(mask_temp_pre,1))*100; 
   
   percent_change_cntrl = [percent_change_cntrl,percent_change]; 
    
end 

end 

%% %% control section

percent_change_NOC = []; 

for trial = 1:4 

for i = 1:9 
    mask_temp_pre = []; 
    mask_temp_post = []; 
    if trial ==1 
    if NOC(1,i) <10
    load([root,num2str(trial), filesep, strcat(trial1_date, '0', num2str(NOC(1,i))),'-01-Mask.mat']); 
    
    else 
    load([root,num2str(trial), filesep, strcat(trial1_date, num2str(NOC(1,i))),'-01-Mask.mat']); 
          
    end 
    end 
     if trial ==2 
    if NOC(1,i) <10
    load([root,num2str(trial), filesep, strcat(trial2_date, '0', num2str(NOC(1,i))),'-01-Mask.mat']); 
    
    else 
    load([root,num2str(trial), filesep, strcat(trial2_date, num2str(NOC(1,i))),'-01-Mask.mat']); 
          
    end 
     end 
     if trial ==3
    if NOC(1,i) <10
    load([root,num2str(trial), filesep, strcat(trial3_date, '0', num2str(NOC(1,i))),'-01-Mask.mat']); 
    
    else 
    load([root,num2str(trial), filesep, strcat(trial3_date, num2str(NOC(1,i))),'-01-Mask.mat']); 
          
    end 
     end 
    
     if trial ==4
    if NOC(1,i) <10
    load([root,num2str(trial), filesep, strcat(trial4_date, '0', num2str(NOC(1,i))),'-01-Mask.mat']); 
    
    else 
    load([root,num2str(trial), filesep, strcat(trial4_date, num2str(NOC(1,i))),'-01-Mask.mat']); 
          
    end 
     end 
    
     mask_temp_pre=mask(:,:,10); 
    mask_temp_pre =mask_temp_pre(mask_temp_pre==1); 
     mask_temp_post=mask(:,:,40); 
    mask_temp_post =mask_temp_post(mask_temp_post==1); 
    
    
%     mask_pre = mask(:,:,10)(mask(:,:,10) ==1); 
%    mask_post = mask(:,:,40)(mask(:,:,40) ==1); 
   
   percent_change = ((size(mask_temp_post,1)-size(mask_temp_pre,1))/size(mask_temp_pre,1))*100; 
   
   percent_change_NOC = [percent_change_NOC,percent_change]; 
    
end 

end 

%%

s = struct; 


s.NOC = percent_change_NOC; 
s.CNTRL = percent_change_cntrl; 

p_val=ranksum(percent_change_NOC, percent_change_cntrl); 

f= figure; 

grid on; 
violinplot(s); 

yticks(-100:10:40); 
ylim([-100, 40]); 
