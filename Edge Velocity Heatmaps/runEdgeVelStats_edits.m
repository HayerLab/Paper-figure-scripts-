  cells = [1,3,4,5,6,7,9,10,11,12,13]; 

 depths = [6]; 
 
 root = 'E:\seph backup\LOK SLK ERM KD\cropped\siERM';
 for i = 1:size(cells,2)
    
     for j= 1:size(depths,2)
         
         
      celldir = ([root, filesep, strcat(num2str(cells(1,i))),filesep, strcat('edge_vels', filesep, 'edge vel mapping_', num2str(6))]);       
         
         
     getEdgeVelStats_edits(celldir,4,-4,-2,25);
     
     
     end 
     
     
     
     
 end 