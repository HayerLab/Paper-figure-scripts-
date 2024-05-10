% single cell analysis before edge tracking
%% initiatialization
clear; clc; 
 root = 'F:\Arnold 240508 single frame immunos T2\t567';
 
 tiffDir = ([root, filesep,'raw']); 
 
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
cell=29; 
filekey = '29_1_1'; 


 cellDir = ([cropdir,filesep, num2str(cell)]); 
     if ~exist(cellDir)
         mkdir(cellDir); 
     end   
     
     
%bgcyto = ([bgpath,filesep, 'AVG_bgcyto.tif']); 
%bgCAAX= ([bgpath,filesep, 'AVG_bgCAAX.tif']); 
%bgCFP = ([bgpath, filesep, 'AVG-BG-CFP.tif']); 
%bgFRET = ([root, filesep, 'background\AVG_bgFRET.tif']);
%bgFRET = ([bgpath,filesep,  'AVG_bgFRET.tif']);
%bgmRuby=([bgpath, filesep,  'AVG_bgmRuby.tif']);
          
 %CAAX = ([tiffDir,filesep, strcat(filekey,'_CAAX_stacked.tif')]);  
 % CFP= ([tiffDir,filesep, strcat(filekey,'_CFP_stacked.tif')]); 
  mRuby= ([tiffDir,filesep, strcat(filekey,'_mRuby_1.tif')]); 
    cyto= ([tiffDir,filesep, strcat(filekey,'_cyto_1.tif')]);  
    pERM = ([tiffDir,filesep, strcat(filekey,'_pERM_1.tif')]);  

%FRET = ([tiffDir,filesep, strcat('230816-03-20-WT-CPD31-YFP-FRET.tif')]);  
%CFP= ([tiffDir,filesep, strcat('230816-03-20-WT-CPD31-CFP.tif')]); 

 %bg_CAAX_image = double(readTIFFstack(bgCAAX)); 
 %bg_CFP_image = double(readTIFFstack(bgCAAX)); 
%bg_mRuby_image = double(readTIFFstack(bgmRuby)); 
%bg_cyto_image = double(readTIFFstack(bgcyto)); 
     
            cropSite = 0;

%                      
%          CAAX_stack=double(readTIFFstack(CAAX));
%          
%          FRET_stack_1 = CAAX_stack(:,:,1); 
%                
%          fg = figure;
%               
%          axis image;  
%              
%             
%       [stackCAAX, cropArea] = serimcropold(CAAX_stack,mean(CAAX_stack,3));
%       
%                 
%       Stack2TIFF(stackCAAX, [cellDir, filesep,'CAAX.tif']);

 mRuby_stack=double(readTIFFstack(mRuby));
         
        mRuby_stack_1 = mRuby_stack(:,:,1); 
               
         fg = figure;
              
         axis image;  
             
            
      [stackmRuby, cropArea] = serimcropold(mRuby_stack,mean(mRuby_stack,3));
      
                
      Stack2TIFF(stackmRuby, [cellDir, filesep,'mRuby.tif']);
                    
%        stackCFP = readFileToStack(CFP); 
%        stackCFP = imcrop3(stackCFP, [cropArea(1), cropArea(2), 1,...
%        cropArea(3), cropArea(4), size(stackCAAX,3)-1]);                    
%        Stack2TIFF(stackCFP, [cellDir, filesep, 'CFP.tif']);
       
       stackcyto = readFileToStack(cyto); 
       stackcyto = imcrop3(stackcyto, [cropArea(1), cropArea(2), 1,...
       cropArea(3), cropArea(4), size(stackcyto,3)-1]);                    
       Stack2TIFF(stackcyto, [cellDir, filesep, 'cyto.tif']);
       
        stackpERM = readFileToStack(pERM); 
       stackpERM = imcrop3(stackpERM, [cropArea(1), cropArea(2), 1,...
       cropArea(3), cropArea(4), size(stackpERM,3)-1]);                    
       Stack2TIFF(stackpERM, [cellDir, filesep, 'pERM.tif']);
            
       
%         stackcyto = readFileToStack(cyto); 
%        stackcyto = imcrop3(stackcyto, [cropArea(1), cropArea(2), 1,...
%        cropArea(3), cropArea(4), size(stackcyto,3)-1]);                    
%        Stack2TIFF(stackcyto, [cellDir, filesep, 'cyto.tif']);
    
                   %  CAAX_bg_crop = imcrop(uint16(bg_CAAX_image), [cropArea(1), cropArea(2), cropArea(3), cropArea(4)]); 
%                    % saveas(FRET_bg_crop, [cellDir, filesep, 'FRET_bg.tif']); 
                  %   imwrite(CAAX_bg_crop,[cellDir, filesep, 'CAAX_bg.tif'] , "WriteMode", "overwrite", "Compression", "none");
%                     
 %                     CFP_bg_crop = imcrop(uint16(bg_CFP_image), [cropArea(1), cropArea(2), cropArea(3), cropArea(4)]); 
%                  %   saveas(CFP_bg_crop, [cellDir, filesep, 'CFP_bg.tif']); 
%             imwrite(CFP_bg_crop,[cellDir, filesep, 'CFP_bg.tif'] , "WriteMode", "overwrite", "Compression", "none");
%             
         %     mRuby_bg_crop = imcrop(uint16(bg_mRuby_image), [cropArea(1), cropArea(2), cropArea(3), cropArea(4)]); 
                 %   saveas(CFP_bg_crop, [cellDir, filesep, 'CFP_bg.tif']); 
      %      imwrite(mRuby_bg_crop,[cellDir, filesep, 'mRuby_bg.tif'] , "WriteMode", "overwrite", "Compression", "none");
            
            
       %       cyto_bg_crop = imcrop(uint16(bg_cyto_image), [cropArea(1), cropArea(2), cropArea(3), cropArea(4)]); 
                 %   saveas(CFP_bg_crop, [cellDir, filesep, 'CFP_bg.tif']); 
       %     imwrite(cyto_bg_crop,[cellDir, filesep, 'cyto_bg.tif'] , "WriteMode", "overwrite", "Compression", "none");
                    
                    %close(fg);
            clear stack;
  close all; 
  
  disp('lets go'); 
    
    
%% background alignment 
    
 for i=4  %par

channels={'CFP' 'FRET'};
 
    cellPath=strcat(num2str(i));
    cellFiles=getFilenames([cropdir,filesep,cellPath, filesep, 'drift_correct'],'.tif');
    
    %note: need to change numbers here depending on if you have an extra
    %channel like myosin- changes order of files in cropped folder
    CFP_stack=double(readTIFFstack([cropdir,filesep,cellPath,filesep,cellFiles{3}]));
    cyto_stack=double(readTIFFstack([cropdir,filesep,cellPath,filesep,cellFiles{4}]));
  
    
alignStack(:,:,2)=(cyto_stack(:,:,60)); % choose arbitrary framenumber, here 75, to generate the alignment parameters 
alignStack(:,:,1)=(CFP_stack(:,:,60));
[pX,pY,dxMat1,dyMat1]=dualviewComputeAlignmentFromGridImages(alignStack);
figure;
subplot(1,2,1); imagesc(dxMat1); colorbar
subplot(1,2,2); imagesc(dyMat1); colorbar

%parsave([cropdir,filesep,cellPath,filesep,'alignment parameters pX pY.mat'],pX,pY,dxMat1,dyMat1);
save([cropdir,filesep,cellPath,filesep,'drift_correct', filesep, 'alignment parameters pX pY.mat'],'pX','pY','dxMat1','dyMat1');


alignStack = []; 
CFP_stack= []; 
cyto_stack = []; 
dxMat1=[]; 
dyMat1 = []; 
end 



%% FRET data 

clc;
%root = 'F:\240301_t567_motility_noco';
cellNum=1;% for now manually select cell folder 

bleachdir=([root,filesep,'data']);
bgdir=[root,filesep,'background'];
  load([bleachdir,filesep,'bleachingcurve.mat']);
  load([bleachdir,filesep,'bleachingcurve_mRuby.mat']);
 load([bgpath, filesep, 'alignment parameters pX pY.mat']); 
 load([bleachdir,filesep,'bleachingcurve_cyto.mat']);


%% Parallel loop
% number of cells you have in a for loop 
for k=1:29
   % rawdir=[root,filesep,'cropped\t567_Cpd31', filesep, strcat( num2str(k))]; 
   rawdir=[cropdir,filesep, strcat( num2str(k))]; 
    %load([rawdir,filesep,'alignment parameters pX pY.mat']);
    
   datadir=[rawdir,filesep,'output'];

 if ~exist(datadir)
    mkdir(datadir);
 end
 
 %choose which one you need here 
 
 % this one has FRET/CFP configured, commented out are options for a 3rd
 % and 4th channel if you want 
%getFRETDataHCS_stacked(k,rawdir,datadir,4); 
 %getFRETDataHCS_stacked_3chan(k,rawdir,datadir,3.5, pX, pY); % FRET, CFP, mRuby
 %getFRETDataHCS_stacked_4_chan(k,rawdir,datadir,1.7, pX, pY); % FRET, CFP, mRuby
getFRETDataHCS_stacked_ezrin_cyto_caax(k,rawdir,datadir,3)


% choose which one you want 
%correctBleachingExp2_stacked_YFP_cyto(fitpara,datadir); %fitpara_mRuby
%correctBleachingExp2_stacked(fitpara, datadir, fitpara_mRuby); %  % does FRET, and mRuby
 % correctBleachingExp2_FRET_stacked(fitpara, datadir); %only does FRET
%correctBleachingExp2_cyto_ratio_stacked(datadir, fitpara_mRuby, fitpara_cyto); % for ezrin ratio calculations 
 
    
end
disp('done!');

%% stupid little function to let you save in a parrallel loop 
function parsave(fname, w,x,y,z)
  save(fname, 'w', 'x', 'y', 'z')
end

