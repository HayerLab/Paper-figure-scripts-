 
root = 'C:\Users\gmarsh8\OneDrive - McGill University\research paper\poster presentation\ezrin rho timelapse'; 
for i =62:102
    
     window = windowCoors{1,i};
   % image=imFRETOutline{1,i}; 
   colorRange = [0.7 1.3]; 
    image=ratio2RGB(imRatio{1,i},colorRange); 
    
   h = figure;  %('visible','off');
      hold on; 
      axis ij; 
     imagesc(image);
     
      for num = 60:15:75
           plot(window(num,2),window(num,1),'.', 'Color','r');
      end 
      
       frame = getframe(h);
       im=frame2im(frame);
%      [imind, cm] = rgb2ind(im,256);

      imwrite(im,[root,filesep,'rho_cell.tif'],'WriteMode','append','Compression','none');
end 