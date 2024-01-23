%% master script
% intialization
%ND2 to TIFF
% tiff stacks 
% background alignment 
% create RatioData_raw and bleaching correction curves 

% to run properly, store ND2 files in folder called 'ND2 files' with
% background.nd2, run1.nd2, run2.nd2 etc. This puts them in alphabetical
% order so the code sorts them properly 
% Seph Marshall Feb 2023

%% setup and initialization. 
clc; 
clear; 

root= 'Y:\Baishali\GCAMP single cell edge tracking 17_01_2024\seph version'; 

% if you want to create Tiff Stacks 
makeTiffStacks = 1; 
 channels = {'gCAMP'};
 channel =  ["gCAMP"]; 
channel1= {'gCAMP'};

frames = 150; 

bgdir = [root, filesep, 'background']; 
if  ~exist(bgdir)
    mkdir(bgdir)
end 

rawdir = [root, filesep, 'raw']; 
if  ~exist(rawdir)
    mkdir(rawdir)
end 

datadir = [root, filesep, 'data']; 
if  ~exist(datadir)
    mkdir(datadir)
end 

ND2dir = [root, filesep, 'ND2 files']; 
if  ~exist(ND2dir)
    mkdir(ND2dir)
end 


if makeTiffStacks ==1 
    TIFFdir = [root, filesep, 'tiff_stacks']; 
if  ~exist(TIFFdir)
    mkdir(TIFFdir)
end 
end 
%% background shit 
 ND2files = getFilenames([ND2dir],'.nd2');
% 
bgfilepath = [ND2dir, filesep, ND2files{1}]; 
counter =0; 


finfo = nd2finfo(bgfilepath);
num_sites = finfo.img_seq_count; 

for row = 1
    for col = 1
        for site =1:num_sites
            
            
          for timept = 1:1
                
            
            shot = [num2str(row),'_',num2str(col),'_',num2str(site)];
            
        if size (channels, 1) == 1
             [gCAMP]= nd2read_hayer_1chan(bgfilepath,finfo, counter,1); 
             imwrite(gCAMP,[bgdir,filesep,shot,'_gCAMP_',num2str(timept),'.tif'],'TIFF','Compression','None');
        end 
  
      
      counter = counter+1;
      
            end   
         end 
        end 
    end 
    

%% files
num_sites = 0; 
 for i =1 % :(size(ND2files,1)-1)
filepath = [ND2dir, filesep, ND2files{i+1}];
 
counter =0; 


finfo = nd2finfo(filepath);
num_sites = finfo.img_seq_count/frames;  

for timept = 1:frames
for row =  i %: %i
    for col = 1
        for site =1:num_sites
       
            
            shot = [num2str(row),'_',num2str(col),'_',num2str(site)];
            
       if size (channels, 1) == 1
             [gCAMP]= nd2read_hayer_1chan(filepath,finfo, counter,1); 
             imwrite(gCAMP,[rawdir,filesep,shot,'_gCAMP_',num2str(timept),'.tif'],'TIFF','Compression','None');
       end
      
      counter = counter+1;
      
            end   
        end 
    end 
end 
 end 
 


%% tiff stacks 
if makeTiffStacks == 1

    
 for row=1 % 1:(size(ND2files,1)-1)
    
    for col = 1
      
writeTiffStacks(root, channel, row, col, num_sites, frames); 
       
    end 
end   
end 

 %% generate background images 


for chan=1:numel(channel1')
    calculate_bg_img_rm_blobs(bgdir,channel1{chan});
    disp(channel1(chan));
    disp(num2str(chan));
end

%% crop file setup


 tiffDir = ([root, filesep,'tiff_stacks']); 

cropdir=[root,filesep,'cropped', filesep, 'cntrl'];
if ~exist(cropdir)
         mkdir(cropdir); 
end  
bgpath=[root,filesep,'background'];
   
    %% crop cell files and accompanying background 
% currently configured for just FRET and CFP, add in channels as needed 

% run this section as many times as needed, just change the cell number 
% which specifies which folder to save under and filekey which specifies 
% which orginal tiff stack to draw from 
   clc; 
cell=1;
filekey = '1_1_1'; 


 cellDir = ([cropdir,filesep, num2str(cell)]); 
     if ~exist(cellDir)
         mkdir(cellDir); 
     end   
     
     
bggCAMP = ([root, filesep, 'background\AVG_bggCAMP.tif']); 
gCAMP= ([tiffDir,filesep, strcat(filekey,'_gCAMP_stacked.tif')]); 

bg_gCAMP_image = double(readTIFFstack(bggCAMP)); 

            cropSite = 0;

                     


 gCAMP_stack=double(readTIFFstack(gCAMP));
         
        gCAMP_stack_1 = gCAMP_stack(:,:,1); 
               
         fg = figure;
              
         axis image;  
             
            
      [stackgCAMP, cropArea] = serimcropold(gCAMP_stack,mean(gCAMP_stack,3));
      
                
      Stack2TIFF(stackgCAMP, [cellDir, filesep,'gCAMP.tif']);
            
    
        gCAMP_bg_crop = imcrop(uint16(bg_gCAMP_image), [cropArea(1), cropArea(2), cropArea(3), cropArea(4)]); 
                 
       imwrite(gCAMP_bg_crop,[cellDir, filesep, 'gCAMP_bg.tif'] , "WriteMode", "overwrite", "Compression", "none");
            
            
           
                    %close(fg);
            clear stack;
  close all; 
  

  %% cell_mask generation
numCells = 1; % specify the number of cells in your current CropDir

Threshold = 2; 
  for i = 1:numCells

 CellDir = ([cropdir, filesep, num2str(i)]);
 dataDir = ([CellDir, filesep, 'output']); 
 if  ~exist(dataDir)
    mkdir(dataDir)
end 

 getMaskData_stacked(i, CellDir,dataDir, Threshold)

  end 



  
 