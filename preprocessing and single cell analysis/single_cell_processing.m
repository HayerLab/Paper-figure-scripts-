% single cell analysis before edge tracking
%% initiatialization
clear; clc; 
 root = 'F:\230919 - Y2motility 2';
 addpath(root); 
 tiffDir = ([root, filesep,'tiff_stacks']); 
 
%  cellDir = ([root, filesep, 'cropped',filesep,'siERM', num2str(cell)]); 
%      if ~exist(cellDir)
%          mkdir(cellDir); 
%      end      
     
cropdir=[root,filesep,'cropped', filesep, 'cntrl_LIS'];
if ~exist(cropdir)
         mkdir(cropdir); 
     end  
bgpath=[root,filesep,'background'];
   
    %% crop cell files and accompanying background 
% currently configured for just FRET and CFP, add in channels as needed 

% run this section as many times as needed, just change the cell number 
% which specifies which folder to save under and filekey which specifies 
% which orginal tiff stack to draw from 
   
cell= 25;
filekey = '1_1_1'; 


 cellDir = ([cropdir,filesep, num2str(cell)]); 
     if ~exist(cellDir)
         mkdir(cellDir); 
     end   
     
     
bgCFP = ([root, filesep, 'background\AVG_bgCFP.tif']); 
bgFRET = ([root, filesep, 'background\AVG_bgFRET.tif']);
          
FRET = ([tiffDir,filesep, strcat(filekey,'_FRET_stacked.tif')]);  
CFP= ([tiffDir,filesep, strcat(filekey,'_CFP_stacked.tif')]); 

bg_FRET_image = double(readTIFFstack(bgFRET)); 
bg_CFP_image = double(readTIFFstack(bgCFP)); 
     
            cropSite = 0;

                     
         FRET_stack=double(readTIFFstack(FRET));
         
         FRET_stack_1 = FRET_stack(:,:,1); 
               
         fg = figure;
              
         axis image;  
             
            
      [stackFRET, cropArea] = serimcropold(FRET_stack,mean(FRET_stack,3));
      
                
      Stack2TIFF(stackFRET, [cellDir, filesep,'FRET.tif']);
                    
                    stackCFP = readFileToStack(CFP); 
                    stackCFP = imcrop3(stackCFP, [cropArea(1), cropArea(2), 1,...
                    cropArea(3), cropArea(4), size(stackFRET,3)-1]);                    
                    Stack2TIFF(stackCFP, [cellDir, filesep, 'CFP.tif']);
            
    
                    FRET_bg_crop = imcrop(uint16(bg_FRET_image), [cropArea(1), cropArea(2), cropArea(3), cropArea(4)]); 
                   % saveas(FRET_bg_crop, [cellDir, filesep, 'FRET_bg.tif']); 
                    imwrite(FRET_bg_crop,[cellDir, filesep, 'FRET_bg.tif'] , "WriteMode", "overwrite", "Compression", "none");
                    
                     CFP_bg_crop = imcrop(bg_CFP_image, [cropArea(1), cropArea(2), cropArea(3), cropArea(4)]); 
                 %   saveas(CFP_bg_crop, [cellDir, filesep, 'CFP_bg.tif']); 
            imwrite(uint16(CFP_bg_crop),[cellDir, filesep, 'CFP_bg.tif'] , "WriteMode", "overwrite", "Compression", "none");
            
                    
                    close(fg);
            clear stack;
  
  
    
    
%% background alignment 
    
for i=1:3

channels={'CFP' 'FRET'};
 
    cellPath=strcat(num2str(i));
    cellFiles=getFilenames([cropdir,filesep,cellPath],'.tif');
    
    %note: need to change numbers here depending on if you have an extra
    %channel like myosin- changes order of files in cropped folder
    CFP_stack=double(readTIFFstack([cropdir,filesep,cellPath,filesep,cellFiles{1}]));
    FRET_stack=double(readTIFFstack([cropdir,filesep,cellPath,filesep,cellFiles{3}]));
  
    
alignStack(:,:,2)=(FRET_stack(:,:,75)); % choose arbitrary framenumber, here 75, to generate the alignment parameters 
alignStack(:,:,1)=(CFP_stack(:,:,75));
[pX,pY,dxMat1,dyMat1]=dualviewComputeAlignmentFromGridImages(alignStack);
figure;
subplot(1,2,1); imagesc(dxMat1); colorbar
subplot(1,2,2); imagesc(dyMat1); colorbar

save([cropdir,filesep,cellPath,filesep,'alignment parameters pX pY.mat'],'pX','pY','dxMat1','dyMat1');

alignStack = []; 
CFP_stack= []; 
FRET_stack = []; 
dxMat1=[]; 
dyMat1 = []; 
end 


%% FRET data 

% clc; clear; 
% root = 'E:\seph backup\LOK SLK ERM KD';
% cellNum=1;% for now manually select cell folder 
% 
% bleachdir=([root,filesep,'data']);
% bgdir=[root,filesep,'background'];
% load([bleachdir,filesep,'bleachingcurve.mat']);
%load([bleachdir,filesep,'bleachingcurve_mRuby.mat']);
% load([bleachdir,filesep,'bleachingcurve_cyto.mat']);


%% Parallel loop
% number of cells you have in a for loop 
for k=2
    rawdir=[root,filesep,'cropped', filesep, strcat(num2str(k))]; 
    load([rawdir,filesep,'alignment parameters pX pY.mat']);
    
   datadir=[rawdir,filesep,'output'];

 if ~exist(datadir)
    mkdir(datadir);
 end
 
 %choose which one you need here 
 
 % this one has FRET/CFP configured, commented out are options for a 3rd
 % and 4th channel if you want 
getFRETDataHCS_stacked(k,rawdir,datadir,1.5); % last number the threshold for segmentation
% getFRETDataHCS_stacked_3chan(k,rawdir,datadir); % FRET, CFP, mRuby
% getFRETDataHCS_stacked_4chan(k,rawdir,datadir); % FRET, CFP, mRuby


% choose which one you want - this is if no drug treatments 
%correctBleachingExp2_stacked_YFP_cyto(fitpara,datadir); %fitpara_mRuby
%  correctBleachingExp2_stacked( fitpara, datadir); %fitparamRuby  % does FRET, and mRuby
  correctBleachingExp_FRET_stacked(fitpara, datadir); %only does FRET
% correctBleachingExp2_cyto_ratio_stacked(datadir, fitpara_mRuby, fitpara_cyto); % for ezrin ratio calculations 
 
    
end
disp('done!');

