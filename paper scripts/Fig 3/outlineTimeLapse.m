%root ='F:\Seph\research paper\Fig 1'; 

colour_array =jet(26);
 
 f1 = figure('visible','on');
      hold on; 
       axis ij;
 
 %change these depending on actual size of image - so it lines up perfectly       
 ylim([0 342]) 
xlim([0 456]);


   for t = [51,58,61,63,65,68,75]%51:75
  
    temp = edgeCoors{1,t};
   plot(temp(:,2), temp(:,1),'Color',colour_array(t-50,:))
 
   end
   
   set(gca,'XColor', 'none','YColor','none')
hold off; 
  %  title('T = 28 min'); 
    
    
% CAN SEE ACTUAL CELL OUTLINE IF NEEDED HERE     
  f2=figure;  
  
 
  
         
 ylim([0 342])
xlim([0 456]);
  axis ij; 
image = imFRETOutline{1,51};
imagesc(image);
hold off; 

 
    


%     ax2 = subplot(2,2,2);
%      ylim([0 250])
% xlim([0 250]);
% hold on; 
%     title('T = 28 min'); 
%     for t = 1:66
%   
%     temp = windowCoors{1,t};
%    plot(temp(:,2), temp(:,1),'Color',colour_array(t,:))
%  
%     end
%     
%     ax3= subplot(2,2,3);
%     title('T = 42 min');
%      ylim([0 250])
% xlim([0 250]);
% hold on; 
% for t = 1:99
%   
%     temp = windowCoors{1,t};
%    plot(temp(:,2), temp(:,1),'Color',colour_array(t,:))
%  
%    
% end
% 
% ax4 = subplot(2,2,4); 
% title('T = 55 min'); 
%  ylim([0 250])
% xlim([0 250]);
% 
% hold on; 
% for t = 1:130
%   
%     temp = windowCoors{1,t};
%    plot(temp(:,2), temp(:,1),'Color',colour_array(t,:))
%  
%     end 
% 
% hold off; 

%saveas(f1,strcat(root,'\','B_i-14min.svg'))