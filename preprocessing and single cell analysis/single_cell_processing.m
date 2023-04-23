% single cell analysis before edge tracking
%% initiatialization
clear; clc; 
 root = 'E:\seph backup\LOK SLK ERM KD';
 
 tiffDir = ([root, filesep,'tiff_stacks']); 
 
 cellDir = ([root, filesep, 'cropped',filesep,'siERM', num2str(cell)]); 
     if ~exist(cellDir)
         mkdir(cellDir); 
     end      
     
cropdir=[root,filesep,'cropped\siNT'];
bgpath=[root,filesep,'background'];
   
    %% crop cell files and accompanying background 
% currently configured for just FRET and CFP, add in channels as needed 

% run this section as many times as needed, just change the cell number 
% which specifies which folder to save under and filekey which specifies 
% which orginal tiff stack to draw from 
   
cell= 1;
filekey = '1_1_14';  

bgCFP = ([root, filesep, 'background\AVG_bgCFP.tif']); 
bgFRET = ([root, filesep, 'background\AVG_bgFRET.tif']);
          
FRET = ([tiffDir,filesep, strcat(filekey,'_FRET_stacked.tif')]);  
CFP= ([tiffDir,filesep, strcat(filekey,'_CFP_stacked.tif')]); 
     
            cropSite = 0;

            
     
         FRET_stack=double(readTIFFstack(FRET));
         
         FRET_stack_1 = FRET_stack(:,:,1); 
               
         fg = figure;
              
         axis image;  
             
            
      [stackFRET, cropArea] = serimcropold(FRET_stack,mean(FRET_stack,3));
      
                
      Stack2TIFF(stackFRET, ['FRET.tif']);
                    
                    stackCFP = readFileToStack(CFP); 
                    stackCFP = imcrop3(stackCFP, [cropArea(1), cropArea(2), 1,...
                    cropArea(3), cropArea(4), size(stack,3)-1]);                    
                    Stack2TIFF(stackCFP, ['CFP.tif']);
            
                    
            close(fg);
            clear stack;
  
  
    close(uifig);
    
%% background alignment 
    
for i= 1:20

channels={'CFP' 'FRET'};
 
    cellPath=strcat(num2str(i));
    cellFiles=getFilenames([cropdir,filesep,cellPath],'.tif');
    
    %note: need to change numbers here depending on if you have an extra
    %channel like myosin- changes order of files in cropped folder
    CFP_stack=double(readTIFFstack([cropdir,filesep,cellPath,filesep,cellFiles{1}]));
    FRET_stack=double(readTIFFstack([cropdir,filesep,cellPath,filesep,cellFiles{2}]));
  
    
alignStack(:,:,2)=(FRET_stack(:,:,75)); % choose arbitrary framenumber, here 75, to generate the alignment parameters 
alignStack(:,:,1)=(CFP_stack(:,:,75));
[pX,pY,dxMat1,dyMat1]=dualviewComputeAlignmentFromGridImages(alignStack);
figure;
subplot(1,2,1); imagesc(dxMat1); colorbar
subplot(1,2,2); imagesc(dyMat1); colorbar

save([cropdir,filesep,cellPath,filesep,'alignment parameters pX pY.mat'],'pX','pY','dxMat1','dyMat1');


end 


%% FRET data 

clc; clear; 
root = 'E:\seph backup\LOK SLK ERM KD';
cellNum=1;% for now manually select cell folder 

bleachdir=([root,filesep,'data']);
bgdir=[root,filesep,'background'];
load([bleachdir,filesep,'bleachingcurve.mat']);
%load([bleachdir,filesep,'bleachingcurve_mRuby.mat']);
% load([bleachdir,filesep,'bleachingcurve_cyto.mat']);


%% Parallel loop
% number of cells you have in a for loop 
for k=20
    rawdir=[root,filesep,'cropped', filesep,'siNT', filesep, strcat( num2str(k))]; 
    load([rawdir,filesep,'alignment parameters pX pY.mat']);
    
   datadir=[rawdir,filesep,'output'];

 if ~exist(datadir)
    mkdir(datadir);
 end
 
 %choose which one you need here 
 
 % this one has FRET/CFP configured, commented out are options for a 3rd
 % and 4th channel if you want 
getFRETDataHCS_stacked(k,rawdir,datadir); 
getFRETDataHCS_stacked_3chan(k,rawdir,datadir); % FRET, CFP, mRuby

% choose which one you want 
%correctBleachingExp2_stacked_YFP_cyto(fitpara,datadir); %fitpara_mRuby
 correctBleachingExp2_stacked( fitpara, datadir); %fitparamRuby  % does FRET, and mRuby
 correctBleachingExp_FRET_stacked(fitpara, datadir); %only does FRET
correctBleachingExp2_cyto_ratio_stacked(datadir, fitpara_mRuby, fitpara_cyto); % for ezrin ratio calculations 
 
    disp(cellNum);
end
disp('done!');

