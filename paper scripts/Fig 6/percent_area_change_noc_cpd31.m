clear; clc; 

root = 'H:\230817\Masks'; 
cntrl_key = '-wT-CTRLD-Mask.mat'; 
NOC_key = '-WT-NCDZL-Mask.mat'; 
NOC_key_2 = '-WT-NCRES-Mask.mat'; 
NOC_cpd31_key ='-WT-CPNOC-Mask.mat'; 
cpd31_key = '-WT-CPD31-Mask.mat'; 


% control section 

percent_change_ctrl = []; 
percent_change_NOC = [];
percent_change_NOC_cpd31 = [];
percent_change_cpd31 = []; 


for row = 1:3

if row  ==2
  continue; 
end 

for   well = 1:4 %well=21:24

load([root, filesep, '230817-0', num2str(row), '-','0', num2str(well), cntrl_key]); 
%load([root, filesep, '230817-0', num2str(row), '-', num2str(well), cntrl_key]); 



mask_pre = mask(:,:,14); 
mask_post = mask(:,:,44); 

 mask_pre = mask_pre(mask_pre ==1); 
   mask_post = mask_post(mask_post ==1); 


percent_change_ctrl = [percent_change_ctrl, ((size(mask_post,1)-size(mask_pre,1))/size(mask_pre,1))*100; ]; 
end 

for well= 5:8 %well =25:28 

load([root, filesep, '230817-0', num2str(row), '-','0', num2str(well), NOC_key])
%load([root, filesep, '230817-0', num2str(row), '-', num2str(well), NOC_key])

    mask_pre = mask(:,:,14); 
    mask_post = mask(:,:,44); 

    mask_pre = mask_pre(mask_pre ==1); 
    mask_post = mask_post(mask_post ==1); 


    percent_change_NOC = [percent_change_NOC, ((size(mask_post,1)-size(mask_pre,1))/size(mask_pre,1))*100; ]; 
     
end 
 
for well= 9:12 %well =29:32
    
    
if well <10
load([root, filesep, '230817-0', num2str(row), '-','0', num2str(well), NOC_key_2])
end 
if well >=10
load([root, filesep, '230817-0', num2str(row), '-', num2str(well), NOC_key_2])
end 
mask_pre = mask(:,:,14); 
mask_post = mask(:,:,44); 

 mask_pre = mask_pre(mask_pre ==1); 
   mask_post = mask_post(mask_post ==1); 


percent_change_NOC = [percent_change_NOC, ((size(mask_post,1)-size(mask_pre,1))/size(mask_pre,1))*100; ]; 
    
    
end 

for  well=13:16   % well= 33:36
    
load([root, filesep, '230817-0', num2str(row), '-', num2str(well), NOC_cpd31_key])
 
mask_pre = mask(:,:,14); 
mask_post = mask(:,:,44); 

 mask_pre = mask_pre(mask_pre ==1); 
   mask_post = mask_post(mask_post ==1); 


percent_change_NOC_cpd31 = [percent_change_NOC_cpd31, ((size(mask_post,1)-size(mask_pre,1))/size(mask_pre,1))*100; ]; 
    
    
end 

for  well= 17:20 %well = 37:40
    
load([root, filesep, '230817-0', num2str(row), '-', num2str(well), cpd31_key])

mask_pre = mask(:,:,14); 
mask_post = mask(:,:,44); 

 mask_pre = mask_pre(mask_pre ==1); 
   mask_post = mask_post(mask_post ==1); 


percent_change_cpd31 = [percent_change_cpd31, ((size(mask_post,1)-size(mask_pre,1))/size(mask_pre,1))*100; ]; 
    
    
end 
end 

cntrl = mean(percent_change_ctrl)
NOC = mean(percent_change_NOC)
NOC_cpd31= mean(percent_change_NOC_cpd31)
cpd31 = mean(percent_change_cpd31)

save(['H:\230817\percent_change_area_stats', filesep, '30 mins post treat Rho WT.mat'], 'percent_change_ctrl',... 
'percent_change_NOC','percent_change_NOC_cpd31',   'percent_change_cpd31'); 

 

