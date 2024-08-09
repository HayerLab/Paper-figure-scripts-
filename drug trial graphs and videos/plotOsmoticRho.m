% Script that plots FRET activity for the different osmotic conditions 
% Seph Marshall 200228

%% Control section 
clc; clear; 
control_root = 'E:\omaima\TEST'; 
num_frames = 20; 
% control_root1= 'F:\Seph\data\230313_20x_2x2_Rho_thrombin_NSC';
%control_root2= 'F:\Seph\data\data_220104 - Trial 5 RhoB Myosin drug treatments';

datadir =('E:\omaima\TEST\graphs');  
%cellFiles=getFilenames([control_root],'_RatioData.mat');

control_arr=zeros(4,num_frames); %here put number of sites and frames
 
k = 0; 
for row = 1
    

    for col =1 
        for site = 1:4

            
           
            k=k+1; 
            position{k}=[num2str(row),'_',num2str(col),'_',num2str(site)];
        end 
    end 
end 


for i = 1:size(position,2)    %(size(cellFiles,1))
   imRatio_raw={}; 
    load([control_root,filesep, 'data', filesep, position{i}, '_RatioData_raw.mat'],'imRatio_raw'); 
    
    for j = 1:num_frames
        
        temp_arr=imRatio_raw{1,j};
        
        control_arr(i,j)= nanmean(temp_arr,'all'); 
         
    end 
    
end

% position = {}; 
% k= 0; 
% for row = 1
% 
% 
%     for col =1 
%         for site = 1:9
% 
% 
% 
%             k=k+1; 
%             position{k}=[num2str(row),'_',num2str(col),'_',num2str(site)];
%         end 
%     end 
% end 
% 
% for i = 1:size(position,2)    %(size(cellFiles,1))
%    imRatio_raw={}; 
%     load([control_root1,filesep, 'data', filesep, position{i}, '_RatioData_raw.mat'],'imRatio_raw'); 
% 
%     for j = 1:60 
%         index = 2*j-1; 
%         temp_arr=imRatio_raw{1,index};
% 
%         control_arr(i+5,j)= nanmean(temp_arr,'all'); 
% 
%     end 
% 
% end 
% k= 0; 
% 
% for row = 1:2
%     
% 
%     for col =1 
%         for site = 1:5
%             
%                       
%            
%             k=k+1; 
%             position{k}=[num2str(row),'_',num2str(col),'_',num2str(site)];
%         end 
%     end 
% end 
% 
% for i = 1:size(position,2)    %(size(cellFiles,1))
%    imRatio_raw={}; 
%     load([control_root2,filesep, 'data', filesep, position{i}, '_RatioData_raw.mat'],'imRatio_raw'); 
%     
%     for j = 1:65 
%         
%         temp_arr=imRatio_raw{1,j};
%         
%         control_arr(i+15,j)= nanmean(temp_arr,'all'); 
%          
%     end 
%     
% end 
% 


%% drug addition setting /
hypo_root = control_root;
% hypo_root1= control_root1;
%hypo_root2= 'F:\Seph\data\data_220104 - Trial 5 RhoB Myosin drug treatments'; 
position = []; 
k = 0; 

hypo_arr=zeros(4,num_frames); 

for row = 1
    
    
    for col =1 
        for site = 5:8

          
          
            k=k+1; 
            position{k}=[num2str(row),'_',num2str(col),'_',num2str(site)];
        end 
    end 
end 




for i = 1:(size(position,2))
   
    load([hypo_root,filesep, 'data', filesep, position{i}, '_RatioData_raw.mat'],'imRatio_raw'); 
    
    
    for j = 1:num_frames
        
        temp_arr2=imRatio_raw{1,j};
        
        hypo_arr(i,j)= nanmean(temp_arr2,'all'); 
         
    end 
    
end 


% position = {}; 
% k= 0; 
% for row = 1
% 
% 
%     for col =1 
%         for site = 10:24
% 
%             k=k+1; 
%             position{k}=[num2str(row),'_',num2str(col),'_',num2str(site)];
%         end 
%     end 
% end 
% 
% 
% 
% 
% for i = 1:(size(position,2))
% 
%     load([hypo_root1,filesep, 'data', filesep, position{i}, '_RatioData_raw.mat'],'imRatio_raw'); 
% 
% 
%     for j = 1:60
% 
%         temp_arr2=imRatio_raw{1,2*j-1};
% 
%         hypo_arr(i+28,j)= nanmean(temp_arr2,'all'); 
% 
%     end 
% 
% end 
% 

% k= 0; 
% for row = 1:2
%     
%     
%     for col =1 
%         for site = 6:10
%           
%             k=k+1; 
%             position{k}=[num2str(row),'_',num2str(col),'_',num2str(site)];
%         end 
%     end 
% end 
% 
% 
% 
% 
% for i = 1:(size(position,2))
%    
%     load([hypo_root2,filesep, 'data', filesep, position{i}, '_RatioData_raw.mat'],'imRatio_raw'); 
%     
%     
%     for j = 1:65
%         
%         temp_arr2=imRatio_raw{1,j};
%         
%         hypo_arr(i+16,j)= nanmean(temp_arr2,'all'); 
%          
%     end 
%     
% end 


%% optional for another drug at same time section 
% hyper_root = 'F:\Seph\data\data_210121 - Trial 1 Rho Y27632_hypotonic\Run2';
% 
% %cellFiles=getFilenames([hyper_root],'_RatioData.mat');
% k= 0; 
% for row = 1
%     for col =1 
%         for site = 8:11
%             k= k+1; 
%             position{k}=[num2str(row),'_',num2str(col),'_',num2str(site)];
%         end 
%     end 
% end 
% hyper_arr=zeros(size(position,2),150); 
% 
% for i = 1:(size(position,2))
%    
%     load([hyper_root,filesep, 'data', filesep, position{i}, '_RatioData_raw.mat'],'imRatio_raw'); 
%     
%     for j = 1:150
%         
%         temp_arr3=imRatio_raw{1,j};
%         
%         hyper_arr(i,j)= nanmean(temp_arr3,'all'); 
%          
%     end 
%     
% end

%% normalize the responses, get averages
% 
frames_preaddition=10; 
 for i=1:size(control_arr,1)
     meanval = mean(control_arr(i,1:frames_preaddition));
     control_arr(i,:)=control_arr(i,:)/meanval; 
 end 
 
 mean_control = mean(control_arr,1); 
 
 for i=1:size(control_arr,1)
    control_arr(i,:) = control_arr(i,:)./mean_control; 
 end 
 
stats_arrCNTRL = zeros(3,num_frames); 

 for k = 1:num_frames
     
     % fitting normal dist to each lag 
    pd = fitdist(control_arr(:,k),'Normal'); 
    stats_arrCNTRL(1,k) = pd.mu; 
    ci_CNTRL = paramci(pd); 
    
    %95 % confidence interval, upper and lower 
  stats_arrCNTRL(2,k) = ci_CNTRL(1,1); 
    stats_arrCNTRL(3,k)=ci_CNTRL(2,1); 
   
 end 



for i=1:size(hypo_arr,1)
    meanval = mean(hypo_arr(i,1:frames_preaddition));
    hypo_arr(i,:)=hypo_arr(i,:)/meanval; 
end 


 
 for i=1:size(hypo_arr,1)
    hypo_arr(i,:) = hypo_arr(i,:)./mean_control; 
 end 
 
 
stats_arrTREAT = zeros(3,num_frames); 

 for k = 1:num_frames 
     
     % fitting normal dist to each lag 
    pd = fitdist(hypo_arr(:,k),'Normal'); 
    stats_arrTREAT(1,k) = pd.mu; 
    ci_TREAT = paramci(pd); 
    
    %95 % confidence interval, upper and lower 
   
    stats_arrTREAT(2,k) = ci_TREAT(1,1); 
    stats_arrTREAT(3,k) = ci_TREAT(2,1); 
 end 

%% plot the data

f1=figure; 


%ylim([0.9 1.9])
xlim([1 num_frames]); 


 xline(frames_preaddition,'--'); 
 ylim([0.8 1.8]);
 hold on; 
 title('RhoA, 20nM u466 at 10 mins'); 
 xlabel('TimePoint (1 min)'); 
 ylabel('Norm. RhoA'); 

for a=1:size(control_arr,1)
    plot([1:num_frames],control_arr(a,:), 'Color','k'); 
    
end 


for c=1:size(hypo_arr,1)
    plot([1:num_frames],hypo_arr(c,:), 'Color','r'); 
end 
 
hold off; 

f2= figure; 
 title('RhoA 20nM U466'); 
 xlabel('TimePoint (1 min)'); 
 ylabel('Norm. RhoA'); 
 hold on; 
 ylim([0.8 1.8]);
 xlim([1 num_frames]);
xline(frames_preaddition,'--'); 




plot(1:num_frames,(stats_arrCNTRL(1,:)),'Color',[0,0,0], 'LineWidth', 3 ); 

plot(1:num_frames,stats_arrCNTRL(2,:),'Color',[0,0,0] );
plot(1:num_frames,stats_arrCNTRL(3,:),'Color',[0,0,0] );

plot(1:num_frames,(stats_arrTREAT(1,:)),'Color',[1,0,0], 'LineWidth', 3 ); 
plot(1:num_frames,(stats_arrTREAT(2,:)),'Color',[1,0,0] );
plot(1:num_frames,(stats_arrTREAT(3,:)),'Color',[1,0,0] );
% 
% hold off; 
% for a=1:size(control_arr,1)
%     plot([1:100],hyper_arr(a,:), 'Color','r','DisplayName','Hyper-osmotic'); 
% end 
% Thrombin = hypo_arr; 
%   save([datadir,filesep,'CNTRL vs. 1 U thrombin @5[0.8 2.2 y axis].mat'],'control_arr', 'Thrombin', 'stats_arrCNTRL', 'stats_arrTREAT'); 
%    saveas(f1, ([datadir, filesep, 'CNTRL vs. 1 U thrombin @5[0.8 2.2 y axis].fig'])); 
%    saveas(f1, ([datadir, filesep, 'CNTRL vs. 1 U thrombin @5[0.8 2.2 y axis].svg'])); 
%   saveas(f2, ([datadir, filesep, 'CNTRL vs. 1 U thrombin @5 avgSD[0.8 2.2 y axis].fig'])); 
%   saveas(f2, ([datadir, filesep, 'CNTRL vs. 1 U thrombin @5 avgSD_[0.8 2.2 y axis].svg'])); 