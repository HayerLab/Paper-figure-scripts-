%Cross correlation Averages 

%% use this section to get the averages out of an entire set 
clc; clear; 
root='F:\old data\data_210225 -  Trial 1 RhoA, RhoA2G\RhoA2G\cropped';
%removed 2 and 35 here to test whats going on 
%cells=[6,7,9,10,12,13,14,15,16,17,18]; %trial 1
%cells = [2,3,5,6,7,8,9,11,12,13]; % trial 2
%cells =[2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19]; %trial 3
%cells = [12,18,19,20,24,26,27,28,29,30,31,32,33]; %CDC42 t1
%cells= [2,3,4,5,7,8,9,11]; % rac trial 1 
%cells =[1,2,4,6,7,8,10,12,13,14,15,16,17,18,19]; %cdc42 t2
%cells = [1, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 17,18, 21, 22]; 
 %cells =[1,2,3,4,6,7,12,13,14,16,17,18,19,20];
 %cells = [1,2,4,5,8,9,10,11]; 
% cells =[1,2,3,6,7,13,14,16,17,18,19] ; 
 %cells =[2,3,4,5,6,7,8,9,10,11,12,13,14,16,28,20,21,23,24,26,27,28,29];
% cells= [1,2,3,4,5,6,7,8,9,10,11,13,14,15,16]; 
% no longer necessary - see google sheets for summary 
cells = [1,2,3,4,5,6,7,8,9,10,11,13,14,15,16];
%cells = [10,15,19]; 
%startarr = [40,30,35]; 
%startarr = [40,40,40,50,40,50,30,35,40,40,30]; %trial 1
%startarr = [30,45,40,40,40,60,50,45,50,40]; %trial 2
%startarr = [30,30,50,30,60,30,40,40,40,35,40,40,50,30,60,40,35]; %trial 3
%startarr = [2,2,2,2,2,2,2,2,2,2,25,2,2]; %CDC42 t1 
%startarr= [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]; 
%startarr =[2,2,2,30,2,2,2,2];% rac trial 1 
%startarr=[30,45,40,40,40,60,50,45,50,40]; 
%startarr= [2,2,2,2,2,2,2,2,2,2,2]; 
cell_arr=cell(1,size(cells,2),1);

%%

depth = 6; 
for loop=1:size(cells,2)
    
     
    
    fileKey=strcat(num2str(cells(1,loop)));
    
    start =  2; %startarr(1,loop);
    
    load([root,filesep,fileKey, filesep, 'edge vel mapping_',num2str(depth),filesep,'Protrusion and FRET values.mat'],'fretvalsF','protvalsWindowF'); %'myosinF','myosin', 'cytoF')
 
% This maps velocity vector from 1-2 with frame 2 of protein expression, etc    
    %edgeVel_arr =fretvalsF(:,start:end);% use this one when FRET myosin is being compared 
    edgeVel_arr =protvalsWindowF(:,start-1:end);   %use this one for when edge vel is the first variable  % can also do a -1 here 
  protExp_arr=fretvalsF(:,start:end);%here change either FRET or myosin
    
 edgeVel_arr(isnan(edgeVel_arr))=0;
  protExp_arr(isnan(protExp_arr))=0; 
    
 %velocity vector created by averaging vectors before and after protein
 %expression frame 
    edgeVel_arr =protvalsWindowF;  
    temporary = NaN(180,size(protvalsWindowF,2)-1);
    for row = 1: size(protvalsWindowF,1)
        for col = 1:size(protvalsWindowF,2)-1
          temporary(row,col)=(protvalsWindowF(row,col)+protvalsWindowF(row,col+1))/2;
        end 
    end 
    edgeVel_arr_adjusted = temporary(:,start-1:end); 
    protExp_arr=fretvalsF(:,start:end-1); %because edge Velocity map now 2 frames smaller than protein expression (first and last frame cut off);
    
   protExp_arr(isnan(protExp_arr))=0; 
    
    
    %using CLT to shift each distribution to a normal distribution
    % Z_edge=( (edgeVel_arr-nanmean(edgeVel_arr(:))) / std(edgeVel_arr(:)));
    Z_edge=( (edgeVel_arr_adjusted-nanmean(edgeVel_arr_adjusted(:))) / std(edgeVel_arr_adjusted(:)));
   
     Z_prot=( (protExp_arr-nanmean(protExp_arr(:))) / std(protExp_arr(:)));
     
    %generating cross correlation for each specific coordinate window 
    
   
    avg_table = NaN(180,size(fretvalsF,2)); 
     for i = 1:size(Z_edge,1) % 
  
     [c,lags] = xcorr((Z_edge(i,:)),Z_prot(i,:), 'coeff');
   
     avg_table(i,1:size(lags,2))=c; 
    
     end 
     
     %compiling average for a given cell  
      for c = 1:size(avg_table,2)
      meanval(1,c) = nanmean(avg_table(:,c)); %seeing here if can be nan mean or just mean 
      end 
      
      %discarding xcorr data outside the range of [-20,20] 
      for i= 1:size(lags,2) 
     
          if (abs(lags(1,i)) > 20)
              meanval(1,i)= NaN;  
   
          end 
          
      end 
      
      %now all mean vals should be centred around zero with range [-20,20] 
      meanval=meanval(~isnan(meanval));
     
      %inputting this into cell array that stores 
      cell_arr{1,loop}= meanval;
  
end 

save(['C:\Users\marsh\OneDrive - McGill University\research paper\results good_Feb2023\Fig 2\RhoA2G\compiled Xcorr depth 1.95 um\Edge Vel. RhoA2G T1 depth 6.mat'],'cell_arr');
   
%%

overall_avg = zeros(1,41); 
for c =1:41
  
    counter =0;
    for cell =1:size(cell_arr,2)
        
 temp = cell_arr{1,cell};
        
  counter = counter + temp(1,c);   
    end 
    
    overall_avg(1,c) = counter/size(cell_arr,2); 
end 

f1=figure; 


hold on; 
%grid on; 
    title('Edge Vel Vs. RhoA2G XCorr');
    xlabel('Lag (min)','FontWeight','bold');
    ylabel('Correlation Coefficient', 'FontWeight','bold'); 
    xline(0, '--'); 
    yline(0,'--');
   
for x= 1:size(cell_arr,2)
    
    plot([-20:20],cell_arr{1,x});
    
end 

plot([-20:20],overall_avg,'Color',[0,0,0], 'LineWidth', 3 ); 

yline(0,'--');
xline(0,'--'); 
xticks([-19.2 -16.8 -14.4 -12.0 -9.6 -7.2 -4.8 -2.4 0 2.4 4.8 7.2 9.6 12.0 14.4 16.8 19.2])
xticklabels({'-8' '' '-6' '' '-4' '' '-2' '' '0' '' '2' '' '4' '' '6' '' '8'}); 

hold off; 

%saveas(f1,strcat(root,'\','graphs depth_',num2str(depth),'\CrossCorrelation','\','RhoA vs. cyto adjusted.png'))
%save([root, '\graphs depth_', num2str(depth),'\CrossCorrelation\RhoA cyto adjusted.mat'],'overall_avg');
  % saveas(f1,strcat('F:\Seph\research paper\Fig 4','\','EdgeVel Actin xcorr+Y2 adjusted.svg'))
  % save(['F:\Seph\research paper\Fig 4\EdgeVel Actin xcorr+Y2 adjusted.mat'],'overall_avg');