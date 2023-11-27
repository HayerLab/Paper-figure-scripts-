% single cell analysis before edge tracking
%% initiatialization
clear; clc; 
 root = 'D:\221209 - 40x 2x2 bin_RhoB_cyto';
 
 tiffDir = ([root, filesep,'tiff_stacks']); 
 
%  cellDir = ([root, filesep, 'cropped',filesep,'siERM', num2str(cell)]); 
%      if ~exist(cellDir)
%          mkdir(cellDir); 
%      end      
     
cropdir=[root,filesep,'cropped'];
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
cell=22;
filekey = '2_1_8'; 


 cellDir = ([cropdir,filesep, num2str(cell)]); 
     if ~exist(cellDir)
         mkdir(cellDir); 
     end   
     
     
bgCFP = ([root, filesep, 'background\AVG_bgCFP.tif']); 
%bgCFP = ([bgpath, filesep, 'AVG-BG-CFP.tif']); 
bgFRET = ([root, filesep, 'background\AVG_bgFRET.tif']);
%bgFRET = ([bgpath, filesep, 'AVG-BG-YFP-FRET.tif']);
bgmRuby=([root, filesep, 'background\AVG_bgmRuby.tif']);
          
 FRET = ([tiffDir,filesep, strcat(filekey,'_FRET_stacked.tif')]);  
 CFP= ([tiffDir,filesep, strcat(filekey,'_CFP_stacked.tif')]); 
  mRuby= ([tiffDir,filesep, strcat(filekey,'_mRuby_stacked.tif')]); 

%FRET = ([tiffDir,filesep, strcat('230816-03-20-WT-CPD31-YFP-FRET.tif')]);  
%CFP= ([tiffDir,filesep, strcat('230816-03-20-WT-CPD31-CFP.tif')]); 

bg_FRET_image = double(readTIFFstack(bgFRET)); 
bg_CFP_image = double(readTIFFstack(bgCFP)); 
bg_mRuby_image = double(readTIFFstack(bgmRuby)); 
     
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
       
       stackmRuby = readFileToStack(mRuby); 
       stackmRuby = imcrop3(stackmRuby, [cropArea(1), cropArea(2), 1,...
       cropArea(3), cropArea(4), size(stackFRET,3)-1]);                    
       Stack2TIFF(stackmRuby, [cellDir, filesep, 'mRuby.tif']);
            
    
                    FRET_bg_crop = imcrop(uint16(bg_FRET_image), [cropArea(1), cropArea(2), cropArea(3), cropArea(4)]); 
                   % saveas(FRET_bg_crop, [cellDir, filesep, 'FRET_bg.tif']); 
                    imwrite(FRET_bg_crop,[cellDir, filesep, 'FRET_bg.tif'] , "WriteMode", "overwrite", "Compression", "none");
                    
                     CFP_bg_crop = imcrop(uint16(bg_CFP_image), [cropArea(1), cropArea(2), cropArea(3), cropArea(4)]); 
                 %   saveas(CFP_bg_crop, [cellDir, filesep, 'CFP_bg.tif']); 
            imwrite(CFP_bg_crop,[cellDir, filesep, 'CFP_bg.tif'] , "WriteMode", "overwrite", "Compression", "none");
            
              mRuby_bg_crop = imcrop(uint16(bg_mRuby_image), [cropArea(1), cropArea(2), cropArea(3), cropArea(4)]); 
                 %   saveas(CFP_bg_crop, [cellDir, filesep, 'CFP_bg.tif']); 
            imwrite(mRuby_bg_crop,[cellDir, filesep, 'mRuby_bg.tif'] , "WriteMode", "overwrite", "Compression", "none");
                    
                    %close(fg);
            clear stack;
  close all; 
  
    
    
%% background alignment 
    
 for i=4  %par

channels={'CFP' 'FRET'};
 
    cellPath=strcat(num2str(i));
    cellFiles=getFilenames([cropdir,filesep,cellPath, filesep, 'drift_correct'],'.tif');
    
    %note: need to change numbers here depending on if you have an extra
    %channel like myosin- changes order of files in cropped folder
    CFP_stack=double(readTIFFstack([cropdir,filesep,cellPath,filesep,cellFiles{3}]));
    FRET_stack=double(readTIFFstack([cropdir,filesep,cellPath,filesep,cellFiles{4}]));
  
    
alignStack(:,:,2)=(FRET_stack(:,:,60)); % choose arbitrary framenumber, here 75, to generate the alignment parameters 
alignStack(:,:,1)=(CFP_stack(:,:,60));
[pX,pY,dxMat1,dyMat1]=dualviewComputeAlignmentFromGridImages(alignStack);
figure;
subplot(1,2,1); imagesc(dxMat1); colorbar
subplot(1,2,2); imagesc(dyMat1); colorbar

%parsave([cropdir,filesep,cellPath,filesep,'alignment parameters pX pY.mat'],pX,pY,dxMat1,dyMat1);
save([cropdir,filesep,cellPath,filesep,'drift_correct', filesep, 'alignment parameters pX pY.mat'],'pX','pY','dxMat1','dyMat1');


alignStack = []; 
CFP_stack= []; 
FRET_stack = []; 
dxMat1=[]; 
dyMat1 = []; 
end 



%% FRET data 

clc;
root = 'D:\221209 - 40x 2x2 bin_RhoB_cyto';
cellNum=1;% for now manually select cell folder 

bleachdir=([root,filesep,'data']);
bgdir=[root,filesep,'background'];
load([bleachdir,filesep,'bleachingcurve.mat']);
load([bleachdir,filesep,'bleachingcurve_mRuby.mat']);
load([bgpath, filesep, 'alignment parameters pX pY.mat']); 
% load([bleachdir,filesep,'bleachingcurve_cyto.mat']);


%% Parallel loop
% number of cells you have in a for loop 
for k=4
    rawdir=[root,filesep,'cropped',  filesep, strcat( num2str(k)), filesep, 'drift_correct']; 
    %load([rawdir,filesep,'alignment parameters pX pY.mat']);
    
   datadir=[rawdir,filesep,'output'];

 if ~exist(datadir)
    mkdir(datadir);
 end
 
 %choose which one you need here 
 
 % this one has FRET/CFP configured, commented out are options for a 3rd
 % and 4th channel if you want 
getFRETDataHCS_stacked(k,rawdir,datadir,4); 
 %getFRETDataHCS_stacked_3chan(k,rawdir,datadir,3.5, pX, pY); % FRET, CFP, mRuby
% getFRETDataHCS_stacked_4_chan(k,rawdir,datadir); % FRET, CFP, mRuby


% choose which one you want 
%correctBleachingExp2_stacked_YFP_cyto(fitpara,datadir); %fitpara_mRuby
correctBleachingExp2_stacked(fitpara, datadir, fitpara_mRuby); %  % does FRET, and mRuby
  %correctBleachingExp2_FRET_stacked(fitpara, datadir); %only does FRET
% correctBleachingExp2_cyto_ratio_stacked(datadir, fitpara_mRuby, fitpara_cyto); % for ezrin ratio calculations 
 
    
end
disp('done!');

%% stupid little function to let you save in a parrallel loop 
function parsave(fname, w,x,y,z)
  save(fname, 'w', 'x', 'y', 'z')
end

