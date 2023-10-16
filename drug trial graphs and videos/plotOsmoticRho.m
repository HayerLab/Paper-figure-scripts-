% Script that plots FRET activity for the different osmotic conditions 
% Seph Marshall 200228

%% Control section 
clc; clear; 
control_root = 'F:\2023105-KomoesinMA-KOmoesin-U46-ROSA-U46-ROSA-MA'; 
%control_root1= 'F:\Seph\data\230313_20x_2x2_Rho_thrombin_NSC';
%control_root2= 'F:\Seph\data\data_220104 - Trial 5 RhoB Myosin drug treatments';

datadir =('F:\2023105-KomoesinMA-KOmoesin-U46-ROSA-U46-ROSA-MA\DRUG RESPOSNE GRAPHS');  
%cellFiles=getFilenames([control_root],'_RatioData.mat');

control_arr=zeros(4,120); 
 
k = 0; 
for row = 1
    

    for col =1 
        for site = 9:12
%             
%             if row ==1 && site==2
%                 continue; 
%             end 
            
           
            k=k+1; 
            position{k}=[num2str(row),'_',num2str(col),'_',num2str(site)];
        end 
    end 
end 


for i = 1:size(position,2)    %(size(cellFiles,1))
   imRatio_raw={}; 
    load([control_root,filesep, 'data', filesep, position{i}, '_RatioData_raw.mat'],'imRatio_raw'); 
    
    for j = 1:120
        
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


%% hypo section /
hypo_root = control_root;
%hypo_root1= control_root1;
%hypo_root2= 'F:\Seph\data\data_220104 - Trial 5 RhoB Myosin drug treatments'; 
position = []; 
k = 0; 

hypo_arr=zeros(8,120); 

for row = 1
    
    
    for col =1 
        for site = 13:20
%             
%              
%             if row ==1 && site <=10
%                 continue; 
%             end 
%             
% %             
            if row ==2 && site >=19
                continue; 
            end 
%           
          
            k=k+1; 
            position{k}=[num2str(row),'_',num2str(col),'_',num2str(site)];
        end 
    end 
end 




for i = 1:(size(position,2))
   
    load([hypo_root,filesep, 'data', filesep, position{i}, '_RatioData_raw.mat'],'imRatio_raw'); 
    
    
    for j = 1:120
        
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


%% hyper section 
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

%% correct the hyper and hypo, get averages
% 
 for i=1:size(control_arr,1)
     meanval = mean(control_arr(i,1:10));
     control_arr(i,:)=control_arr(i,:)/meanval; 
 end 
 
 mean_control = mean(control_arr,1); 
 
 for i=1:size(control_arr,1)
    control_arr(i,:) = control_arr(i,:)./mean_control; 
 end 
 
stats_arrCNTRL = zeros(3,120); 

 for k = 1:120
     
     % fitting normal dist to each lag 
    pd = fitdist(control_arr(:,k),'Normal'); 
    stats_arrCNTRL(1,k) = pd.mu; 
    ci_CNTRL = paramci(pd); 
    
    %95 % confidence interval, upper and lower 
  stats_arrCNTRL(2,k) = ci_CNTRL(1,1); 
    stats_arrCNTRL(3,k)=ci_CNTRL(2,1); 
   
 end 



for i=1:size(hypo_arr,1)
    meanval = mean(hypo_arr(i,1:10));
    hypo_arr(i,:)=hypo_arr(i,:)/meanval; 
end 


 
 for i=1:size(hypo_arr,1)
    hypo_arr(i,:) = hypo_arr(i,:)./mean_control; 
 end 
 
 
stats_arrTREAT = zeros(3,120); 

 for k = 1:120 
     
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
xlim([0 120]); 


 xline(10,'--'); 
 ylim([0.8  1.7]);
 hold on; 
 title('Hs578t ROSA DORA RhoB Moesin KO, cntrl MA vs 20 nM U46619 @ 10 minutes'); 
 xlabel('TimePoint (30s)'); 
 ylabel('Norm. Rho'); 

for a=1:size(control_arr,1)
    plot([1:120],control_arr(a,:), 'Color','k','DisplayName','CNTRL'); 
    
end 
%     
% for b=1:size(hyper_arr,1)
%     plot([1:150],hyper_arr(b,:), 'Color','r','DisplayName','hypertonic buffer @15'); 
% end

for c=1:size(hypo_arr,1)
    plot([1:120],hypo_arr(c,:), 'Color','r','DisplayName','+ML7 @5'); 
end 
 
hold off; 

f2= figure; 
 title('RhoA,1 U Thrombin'); 
 xlabel('TimePoint (1 min)'); 
 ylabel('Norm. RhoA'); 
 hold on; 
 ylim([0.8 2.2]);
 xlim([0 120]);
xline(5,'--'); 

% for i=1:65 
%     
%     sd_up_CNTRL(1,i) =stats_arrCNTRL(1,i)+stats_arrCNTRL(2,i); 
%     sd_down_CNTRL(1,i) =stats_arrCNTRL(1,i)-stats_arrCNTRL(2,i);
%     
%     sd_up_TREAT(1,i) =stats_arrTREAT(1,i)+stats_arrTREAT(2,i); 
%     sd_down_TREAT(1,i) =stats_arrTREAT(1,i)-stats_arrTREAT(2,i);
%     
% end 
 

plot(1:120,(stats_arrCNTRL(1,:)),'Color',[0,0,0], 'LineWidth', 3 ); 

plot(1:120,stats_arrCNTRL(2,:),'Color',[0,0,0] );
plot(1:120,stats_arrCNTRL(3,:),'Color',[0,0,0] );

plot(1:120,(stats_arrTREAT(1,:)),'Color',[1,0,0], 'LineWidth', 3 ); 
plot(1:120,(stats_arrTREAT(2,:)),'Color',[1,0,0] );
plot(1:120,(stats_arrTREAT(3,:)),'Color',[1,0,0] );
% 
% hold off; 
% for a=1:size(control_arr,1)
%     plot([1:100],hyper_arr(a,:), 'Color','r','DisplayName','Hyper-osmotic'); 
% end 
 U46619 = hypo_arr; 
   save([datadir,filesep,'CNTRL vs. 20uM U46619 @10_MoesinKO_[0.8 1.6 y axis].mat'],'control_arr', 'U46619', 'stats_arrCNTRL', 'stats_arrTREAT'); 
   saveas(f1, ([datadir, filesep, 'CNTRL vs. 20uM U46619 @10_MoesinKO_[0.8 1.6 y axis].fig'])); 
    saveas(f1, ([datadir, filesep, 'CNTRL vs. 20uM U46619 @10_MoesinKO_[0.8 1.6 y axis].svg'])); 
%   saveas(f2, ([datadir, filesep, 'CNTRL vs. 1 U thrombin @5 avgSD[0.8 2.2 y axis].fig'])); 
%   saveas(f2, ([datadir, filesep, 'CNTRL vs. 1 U thrombin @5 avgSD_[0.8 2.2 y axis].svg'])); 