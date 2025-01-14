%% calculates granularity of actin or myosin of an ROI before/ after drug treatment 
%SM December 2022 
clear; clc; 

root = 'F:\231113_RhoB_ezrin_ WTvsT567_10umNOC\cropped\WT\control\2\output'; 
datadir= 'F:\231113_RhoB_ezrin_ WTvsT567_10umNOC\cropped\WT\control\2\output'; 
%datadir1 = ([datadir, filesep, 'data']); 
%% control scection
k = 0; 
for row =1
    for col=1
        for site=1:16

            k=k+1;
            position{k}=[num2str(row),'_',num2str(col),'_',num2str(site)];
        end
    end
end
 
for i= 1 %1:size(position,2)
    
     average =[]; 
     average_m = []; 
    % load([root, filesep, position{i}, '_RatioData_raw.mat'])
     load([root, filesep, 'RatioData_raw.mat'])
     
 
 % number of control frames you have before drug addition
     for k = 1:5
         

       % FRET_mean = nanmean(imRatio_raw{1,k}, 'all'); 
      ezrin_mean = nanmean(imEzrin_raw{1,k}, 'all'); 
        average = [average; ezrin_mean];
       % 
       % myosin_mean = nanmean(im_mRuby_raw{1,k}, 'all'); 
       %  average_m = [average_m;myosin_mean];
     end
     
     mean_FRET = mean(average);
     %mean_m = mean(average_m); 
        
        colorRange1 = [0.5 1.5]; 
        colorRange2 = [0 3]; 
       
        %correct the images by normalizing them to pre-drug addition
        %averages 
    for j=1:size(imEzrin_raw,2)      %size(imRatio_raw,2)
    temp = imEzrin_raw{1,j};
            ezrin_norm{1,j}=temp./mean_FRET;
            tempest= ezrin_norm{1,j}; 
        
        
%             temp = imRatio_raw{1,j};
%             FRET_norm{1,j}=temp./mean_FRET;
%             tempest= FRET_norm{1,j}; 
            

            %temp2 = im_mRuby_raw{1,j}; 
            % mRuby_norm{1,j} = temp2./mean_m; 
            % tempest1= mRuby_norm{1,j}; 
            % tempest1(tempest1 <=0) =0.001; 
           
          tempRATIOforstack=ratio2RGB(tempest,colorRange1);%Cdc42
         %   tempmRubyforstack=ratio2RGB(tempest1,colorRange2);%Cdc42            
           % tempmRubyforstack(tempmRubyforstack < 0) = 0.01; 
    
         
          imwrite(tempRATIOforstack,[datadir,filesep, position{i},  '_ezrin_',num2str(colorRange1(1)),'_',num2str(colorRange1(2)),'.tif'],'WriteMode','append','Compression','none');     %
     %      imwrite(tempmRubyforstack,[datadir, filesep, position{1,1}, '_mRuby__redo',num2str(colorRange2(1)),'_',num2str(colorRange2(2)),'.tif'],'WriteMode','append','Compression','none');    %
    end

       % save([datadir,filesep,'data',filesep, position{1},'_RatioData.mat'],'FRET_norm', 'mRuby_norm', 'colorRange1', 'colorRange2');
end 
        
%% specify ROI of interest of larger image 

cropdir = ([datadir, filesep, 'granularity_ROIs']); 
 
% as specified in ImageJ ROI manager
    x_range = [725,955]; 
    y_range = [228,551]; 
  
  crop_Rho = cell(1,65); 
  crop_myosin = cell(1, 65); 
   crop_myosin_filter = cell(1,65); 
 
     colorRange1 = [0 2]; 
     colorRange2 = [0 3]; 
       
  for j=1:size(imRatio_raw,2)
         temp = FRET_norm{1,j};
%          temp=temp./mean_FRET;  
         crop_Rho{1,j} = temp(y_range(1,1):y_range(1,2), x_range(1,1):x_range(1,2));
         temp2 = mRuby_norm{1,j}; 
       %  temp2 = temp2./mean_m; 
         crop_myosin{1,j} = temp2(y_range(1,1):y_range(1,2), x_range(1,1):x_range(1,2));
    
         % gaussian filtering to reduce salt and peppper noise 
         crop_myosin_filter{1,j} = imgaussfilt(crop_myosin{1,j}, 1); 
         
          tempRATIOforstack1=ratio2RGB(crop_Rho{1,j},colorRange1);%Cdc42
         tempmRubyforstack1=ratio2RGB( crop_myosin{1,j},colorRange2);
         tempmRubyforstack1(tempmRubyforstack1 < 0) = 0; 
         
         tempmRubyforstack2=ratio2RGB( crop_myosin_filter{1,j},colorRange2);
         tempmRubyforstack2(tempmRubyforstack1 < 0) = 0; 
         
          imwrite(tempRATIOforstack1,[cropdir,filesep,position{1,1}, '_ROI_1_RhoFRET ','[',num2str(colorRange1(1)),'_',num2str(colorRange1(2)),']','.tif'],'WriteMode','append','Compression','none')
          imwrite(tempmRubyforstack1,[cropdir,filesep,position{1,1}, '_ROI_1_mRuby ','[',num2str(colorRange2(1)),'_',num2str(colorRange2(2)),']','.tif'],'WriteMode','append','Compression','none')
          imwrite(tempmRubyforstack2,[cropdir,filesep,position{1,1}, '_ROI_1_mRuby_gaussian ','[',num2str(colorRange2(1)),'_',num2str(colorRange2(2)),']','.tif'],'WriteMode','append','Compression','none')
  end
  
  %frames = [5, 15, 25, 35, 45];
  frames = [5, 6, 10, 15, 25];
  
  granularity_data = cell(5,4); 
  granularity_data_gaussian = cell(5,4); 
  
  
  for i = 1:size(frames,2)
       % no gaussian blur 
      granularity_data{i,1}= strcat('Frame_', num2str(frames(1,i))); 
      %removing nan values from the array 
     temp =isnan (crop_myosin{1,frames(1,i)}); 
     final = crop_myosin{1,frames(1,i)}(~temp); 
  %standard deviation of pixel distribution
  granularity_data{i,2} = std(final); 
  
  %entropy of the image 
  granularity_data{i,3} = entropy(final); 
 
  %gray co-occurence matrix 
  crop_myosin{1,frames(1,i)}(isnan(crop_myosin{1,frames(1,i)})) =0; 
  gcm = graycomatrix(crop_myosin{1,frames(1,i)}); 
  granularity_data{i,4} = graycoprops(gcm,'all');  
  
   %%%%%%gaussian blur applied
      granularity_data_gaussian{i,1}= strcat('Frame_', num2str(frames(1,i))); 
      %removing nan values from the array 
     temp =isnan (crop_myosin_filter{1,frames(1,i)}); 
     final = crop_myosin_filter{1,frames(1,i)}(~temp); 
  %standard deviation of pixel distribution
  granularity_data_gaussian{i,2} = std(final); 
  
  %entropy of the image 
  granularity_data_gaussian{i,3} = entropy(final); 
 
  %gray co-occurence matrix 
  crop_myosin_filter{1,frames(1,i)}(isnan(crop_myosin_filter{1,frames(1,i)})) =0; 
  gcm = graycomatrix(crop_myosin_filter{1,frames(1,i)}); 
  granularity_data_gaussian{i,4} = graycoprops(gcm,'all'); 
  
  FRET_average(1, i) = nanmean(crop_Rho{1,frames(1,i)}, 'all'); 
  end 
  
 save([cropdir,filesep,position{1},'_ROI_1_data.mat'],'x_range','y_range','granularity_data','granularity_data_gaussian','FRET_average','crop_Rho', 'crop_myosin');
  