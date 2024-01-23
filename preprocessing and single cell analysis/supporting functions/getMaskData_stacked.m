function getMaskData_stacked(cellNum,CellDir,datadir, threshold)


% Used non built-in subfunctions: 
%
%   dualviewAlignFromFittedSurface
%   getBGMask
%   subBG
%   getCellMask
%   DrawMaskOutline
%   ratio2RGB
%   vect
%
% Arnold, 24 May 2015
if ~exist([datadir,filesep,'MaskData_raw.mat'])
%%%%%% Set up


%%%%%% Call background images

gCAMPbg_raw=double(imread([CellDir,filesep,'gCAMP_bg.tif']));
        



%cellFiles=double(imread([CellDir,filesep,'gCAMP.tif']));

%check here order of cellFiles!!! 
 gCAMP_stack=double(readTIFFstack([CellDir,filesep,'gCAMP.tif']));
 
 
  
%%%%%% Loop through frames
maskFinal={};cellCoors={};  
for frameNum=1:size(gCAMP_stack,3)
    disp(num2str(frameNum));
    imgCAMP_raw=gCAMP_stack(:,:,frameNum);
 
    %%%%%% Background-subtract 
    bgmask=getBGMask(imgCAMP_raw);
   

 imgCAMPbg=subBG(imgCAMP_raw,bgmask,gCAMPbg_raw);

 
    %%%%%% Get mask from raw FRET image

    s= cellNum; 
   [mask cellCoorsTemp]=getCellMaskCyto_2_stacked((2*imgCAMPbg),1000, frameNum,s, threshold); % see here if its better taking away the imCFPbg
   
    maskFinal{frameNum}=mask;
   cellCoors{frameNum}=cellCoorsTemp;
    %%%%%% Detrmine ratio
    

    gCAMPOutline{frameNum}=DrawMaskOutline(imgCAMP_raw,mask);
   imwrite(gCAMPOutline{frameNum},[datadir,filesep,'Outline_prelim.tif'],'WriteMode','append','Compression','none');
   
end
%%%%%% Bleaching correction: Detrmine linear fit parameters for FRET/CFP decay

for frameNum= 1:length(imgCAMP_raw)
   bleach_raw(frameNum)=nanmean(vect(imgCAMP_raw(frameNum)));
   
end
save([datadir,filesep,'MaskData_raw.mat'],'maskFinal','cellCoors','gCAMPOutline','-v7.3'); 
save([datadir,filesep,'Bleach_raw.mat'],'bleach_raw'); 
end

