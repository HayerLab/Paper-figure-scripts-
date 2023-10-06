% Quantify protrusion/retraction dynamics
clc; clear; 
root = 'F:\230914 - Y2 motility 1\cropped\20uM_Y2'; 
 
%cells =[2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19]; 
%start_arr = [30,30,50,30,60,30,40,40,40,35,40,40,50,30,60,40,35]; 
stats =struct; 
index =1; 
for i= 1:25%size(cells, 2)
   
     
 if i == 1 || i == 9 || i ==17  %|| i == 23 || i == 24 || i ==27
          continue; 
end 
  

   load([root, filesep, num2str(i),filesep,'output', filesep, 'edge_vels', filesep, 'edge vel mapping_3', filesep, 'protrusionlist_over50.mat']); 
   load([root, filesep, num2str(i),filesep,'output', filesep, 'edge_vels', filesep, 'edge vel mapping_3', filesep, 'retractionlist_over50.mat']);
   load([root, filesep, num2str(i),filesep, 'output', filesep, 'edge_vels', filesep, 'edge vel mapping_3',  filesep, 'protrusion and FRET values.mat'], 'protvalsWindowF');
    
  % accounts for some cells spreading first, where events aren't counted.
  % This standardizes events to number occuring in a 1hr time frame 
   stats.cellNum(index) = i; 
  scale_factor = size(protvalsWindowF,2);   %-start_arr(1,i); 
  
  if size(retractions,1) ==0
      stats.cellNum(index) = i; 
      stats.num_retr(index)= 0;
      stats.ret_avg_size(index) =nan; 
      stats.ret_speed(index) = nan; 
      stats.retr_perimeter(index) = nan; 
    stats.retr_timelength(index) =nan; 
    
       stats.num_prots(index) =size(protrusions,1)*(144/scale_factor); 
       stats.prot_avg_size(index) = mean([protrusions.Area]);
        stats.prot_speed(index) = mean([protrusions.mspeed])/0.641; 
        
   for j = 1:size(protrusions,1)
       temp(1,j) = (protrusions(j).coorend - protrusions(j).coorstart)/180; 
       temp2(1,j) = (protrusions(j).frameend - protrusions(j).framestart)*(25/60); % turns frames into minutes 
   end 
   
    stats.prot_perimeter(index) = mean(temp); 
    stats.prot_timelength(index) = mean(temp2); 
 index = index+1; 
  else 
  
    stats.num_prots(index) =size(protrusions,1)*(144/scale_factor); 
   stats.num_retr(index)= size(retractions,1)*(144/scale_factor); 
   
   stats.prot_avg_size(index) = mean([protrusions.Area]);
   stats.ret_avg_size(index) = mean([retractions.Area]); 
   
   stats.prot_speed(index) = mean([protrusions.mspeed])/0.641;  %/1.282; % this converts px/25sframe to um/min @ 0.325 um/px
   stats.ret_speed(index) = mean([retractions.mspeed])/0.641;     %/1.282; % 0.641 is for 0.65 um/px
   
   for j = 1:size(protrusions,1)
       temp(1,j) = (protrusions(j).coorend - protrusions(j).coorstart)/180; 
       temp2(1,j) = (protrusions(j).frameend - protrusions(j).framestart)*(25/60); % turns frames into minutes 
   end 
   
    stats.prot_perimeter(index) = mean(temp); 
    stats.prot_timelength(index) = mean(temp2); 
    
    for k = 1:size(retractions,1)
       temp3(1,k) = (retractions(k).coorend - retractions(k).coorstart)/180; 
       temp4(1,k) = (retractions(k).frameend - retractions(k).framestart)*(25/60); % turns frames into minutes 
   end 
   
    stats.retr_perimeter(index) = mean(temp3); 
    stats.retr_timelength(index) = mean(temp4); 
   
   index =index+1; 
end 
end 
 save([root, filesep, 'stats_50px_thresh'],  'stats'); 
