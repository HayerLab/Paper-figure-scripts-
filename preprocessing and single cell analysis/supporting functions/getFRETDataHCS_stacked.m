function getFRETDataHCS_stacked(cellNum,rawdir,datadir)


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
if ~exist([datadir,filesep,'RatioData_raw.mat'])
%%%%%% Set up

load([rawdir,filesep,'alignment parameters pX pY.mat'],'pX','pY');

%%%%%% Call background images
binning=1; % relevant if alingment images and data images were acquired using distinct binning settings
CFPbg_raw=double(imread([rawdir,filesep,'bgCFP.tif']));
FRETbg_raw=double(imread([rawdir,filesep,'bgFRET.tif']));

bg1(:,:,1)=CFPbg_raw; bg1(:,:,2)=FRETbg_raw;
bg2=dualviewAlignFromFittedSurface(bg1,pX,pY,binning);
CFPbg=bg2(:,:,1);
FRETbg=bg2(:,:,2);


cellPath=strcat('cell_',num2str(cellNum));
cellFiles=getFilenames([rawdir],'.tif');

%check here order of cellFiles!!! 
 CFP_stack=double(readTIFFstack([rawdir,filesep,cellFiles{1}]));
 FRET_stack=double(readTIFFstack([rawdir,filesep,cellFiles{2}]));
 
 
  
%%%%%% Loop through frames
imRatio_raw={};maskFinal={};cellCoors={};  
for frameNum=1:size(CFP_stack,3)
    disp(num2str(frameNum));
    imCFP_raw=CFP_stack(:,:,frameNum);
    imFRET_raw=FRET_stack(:,:,frameNum);

    %%%%%% Align CFP/FRET images
    imstack(:,:,1)=imCFP_raw; imstack(:,:,2)=imFRET_raw;
    imaligned=dualviewAlignFromFittedSurface(imstack,pX,pY,1);
    imCFP_raw=imaligned(:,:,1);
    imFRET_raw=imaligned(:,:,2);
    
    %%%%%% Background-subtract CFP/FRET images
    bgmask=getBGMask(imCFP_raw+imFRET_raw);
   

 imCFPbg=subBG(imCFP_raw,bgmask,CFPbg);
 imFRETbg=subBG(imFRET_raw,bgmask,FRETbg);
 
    %%%%%% Get mask from raw FRET image
    
    s= cellNum; 
   [mask cellCoorsTemp]=getCellMaskCyto_2_stacked((2*imFRETbg),1000, frameNum,s); % see here if its better taking away the imCFPbg
   
    maskFinal{frameNum}=mask;
   cellCoors{frameNum}=cellCoorsTemp;
    %%%%%% Detrmine ratio
    imFRETbg(~mask)=nan;
    imFRET=ndnanfilter(imFRETbg,fspecial('disk',3),'replicate');
    imCFPbg(~mask)=nan;
    imCFP=ndnanfilter(imCFPbg,fspecial('disk',3),'replicate');
    imRatioTemp=imFRET./imCFP;
    imRatioTemp(~mask)=nan;
    imRatio_raw{frameNum}=imRatioTemp;
    

    %%%%%% Determine scaling for representation
    if frameNum==1
       colorRange = [round(prctile(imRatioTemp(:),1),1),round(prctile(imRatioTemp(:),98),1)];
    end
    %%%%%% Generate and write files for raw ratio and outlined objects - optional
    %%%%%% 
   % tempRATIO=ratio2RGB(imRatioTemp,colorRange);
    imFRETOutline{frameNum}=DrawMaskOutline(imFRET_raw,mask);
   imwrite(imFRETOutline{frameNum},[datadir,filesep,'Outline_prelim.tif'],'WriteMode','append','Compression','none');
   
end
%%%%%% Bleaching correction: Detrmine linear fit parameters for FRET/CFP decay

for frameNum= 1:length(imRatio_raw)
   bleach_raw(frameNum)=nanmean(vect(imRatio_raw{frameNum}));
end
save([datadir,filesep,'RatioData_raw.mat'],'maskFinal','cellCoors','imRatio_raw','imFRETOutline','-v7.3'); 
save([datadir,filesep,'Bleach_raw.mat'],'bleach_raw'); 
end

