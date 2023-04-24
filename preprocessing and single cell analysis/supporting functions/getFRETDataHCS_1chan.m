function getFRETDataHCS_1chan(position,bgdir,rawdir,datadir,threshold)
%function getFRETDataHCS( row,col,site )
% image processing for FRET data analysis 
% row=3;col=3;site=2;
% Input: CFP/FRET images, background images,  and alignment parameters for
% CFP/FRET images. 
%
% Generates first a mask for background subtraction, subtracts background,
% and calculates a refined mask before taking a ratio between background-
% subtracted CFP and FRET images. Bleaching correction assumes constant 
% mean FRET values per frame across the entire time-lapse series.
%
% Output: 
%
%   maskFinal   Cell array with binary mask for each frame
%
%   cellCoors   Cell array with centroid x,y coordinates of detected
%               objects and object size in pixels
%
%   imRatioraw  Cell array with masked raw ratio values, before bleaching
%               correction
%
%   imRatio     Cell array with bleaching-corrected masked ratio values
%
%   imFRETOUtline_row_col.tif 
%               RGB tif stack of FRET channel images with outlined detected
%               objects
%
%   imRatio_row_col_colorRangeLow_colorRangeHigh.tif
%               RGB tif stack with false-colored (jet-black) scaled ratio
%               images, bleaching corrected.
%
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
key=[datadir,filesep,position,'_RatioData_raw.mat'];
if ~exist(key)
%%%%%% Set up


%%%%%% Call background images
 binning=1; % relevant if alingment images and data images were acquired using distinct binning settings
%mCitbg=double(imread([bgdir,filesep,'AVG_bgmCit.tif']));


CFP_files=dir([rawdir,filesep,position,'_mCit_*']);

%%%%%% Loop through frames
im_mCit={};maskFinal={};cellCoors={}; imOutline = {};   
for frameNum=1:length(CFP_files)
   
    disp([position,'__',num2str(frameNum)]);
    im_mCit_raw=double(imread([rawdir,filesep,position,'_mCit_',num2str(frameNum),'.tif']));
   
  

    %%%%%% Background-subtract CFP/FRET images
%     bgmask=getBGMask(im_mCit_raw);
%     im_mCitbg = subBG(im_mCit_raw,bgmask,mCitbg); 
  
  
    im_mCitbg = im_mCit_raw; 
    %%%%%% Get mask from raw FRET image
    [mask cellCoorsTemp]=getCellMaskCyto_2(im_mCitbg,1500, frameNum, position); 
   
    % changing this from imFRET_raw and imCFP_raw, which are required if you are using getCellMaskCyto_edits
    maskFinal{frameNum}=mask;
    cellCoors{frameNum}=cellCoorsTemp;
   
     im_mCitbg(~mask)=nan;
     % if you want to filter the channel 
  %  im_mCit_raw_=ndnanfilter(im_mCitbg,fspecial('disk',3),'replicate');
   im_mCit{frameNum}=im_mCitbg;
    
    
    

    %%%%%% Generate and write files for raw ratio and outlined objects
   
    imOutline{frameNum}=DrawMaskOutline(im_mCit_raw,mask);
    
    imwrite(imOutline{frameNum},[datadir,filesep,position,'_Outline_prelim.tif'],'WriteMode','append','Compression','none');
end


%%%%%% Bleaching correction: Detrmine linear fit parameters for FRET/CFP decay
for frameNum=1:length(imRatio_raw)
   bleach_raw(frameNum)=nanmean(vect(im_mCit{frameNum}));

save([datadir,filesep,position,'_RatioData_raw.mat'],'maskFinal','cellCoors','im_mCit','imOutline','-v7.3');  

save([datadir,filesep,position,'_Bleach_raw.mat'],'bleach_raw' );  
end  

end 
end 

