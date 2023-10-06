  cells = [1,3,4,5,6,7,9,10,11,12,13]; 

 depths = [3]; 
 
 root = 'F:\230919 - Y2motility 2\cropped\20uM_Y2';
 for i = 2:26
    
%      if i == 7 || i == 14 || i ==20
%          continue; 
%      end 
     for j= 1:size(depths,2)
         
         
      celldir = ([root, filesep, strcat(num2str(i)),filesep, strcat('output', filesep,'edge_vels', filesep, 'edge vel mapping_', num2str(3))]);       
         
         
     getEdgeVelStats_edits(celldir,1,2.5,-2.5,25);
     
     
     end 
     
     
     
     
 end 