% Quantify protrusion/retraction dynamics

root = 'F:\Seph\data\data_200116 - Trial 3 Rho, Myosin\cropped'; 
datadir = 
cells =[2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19]; 
start_arr = [30,30,50,30,60,30,40,40,40,35,40,40,50,30,60,40,35]; 
stats =struct; 

for i= 1:size(cells, 2)
    
   load([root, filesep, 'cell_', num2str(cells(1,i)), filesep, 'edge vel mapping_6', filesep, 'protrusionlist.mat']); 
   load([root, filesep, 'cell_', num2str(cells(1,i)), filesep, 'edge vel mapping_6', filesep, 'retractionlist.mat']);
    load([root, filesep, 'cell_', num2str(cells(1,i)), filesep, 'edge vel mapping_6', filesep, 'protrusion and FRET values.mat'], 'protvalsWindowF');
    
  % accounts for some cells spreading first, where events aren't counted.
  % This standardizes events to number occuring in a 1hr time frame 
   stats.cellNum(i) = cells(1,i); 
  
  scale_factor = size(protvalsWindowF,2)-start_arr(1,i); 
  
    stats.num_prots(i) =size(protrusions,1)*(144/scale_factor); 
   stats.num_retr(i)= size(retractions,1)*(144/scale_factor); 
   
   stats.prot_avg_size(i) = mean([protrusions.Area]);
   stats.ret_avg_size(i) = mean([retractions.Area]); 
   
   stats.prot_speed(i) = mean([protrusions.mspeed])/1.282; % this converts px/25sframe to um/min
   stats.ret_speed(i) = mean([protrusions.mspeed])/1.282; 
   
   for j = 1:size(protrusions,1)
       temp(1,j) = (protrusions(j).coorend - protrusions(j).coorstart)/180; 
       temp2(1,j) = (protrusions(j).frameend - protrusions(j).framestart)*(25/60); % turns frames into minutes 
   end 
   
    stats.prot_perimeter(i) = mean(temp); 
    stats.prot_timelength(i) = mean(temp2); 
    
    for k = 1:size(retractions,1)
       temp3(1,k) = (retractions(k).coorend - retractions(k).coorstart)/180; 
       temp4(1,k) = (retractions(k).frameend - retractions(k).framestart)*(25/60); % turns frames into minutes 
   end 
   
    stats.retr_perimeter(i) = mean(temp3); 
    stats.retr_timelength(i) = mean(temp4); 
   
   
end 