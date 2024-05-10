function correctBleachingExp2_FRET_stacked(fitpara,datadir)
% What it does
% Normalizes imRatio_raw for each frame by the median intensity per frame of the entire
% time series. 
% bleaching correction by division by linear bleaching function normalized
% by the median of its values. 

%remember to add in/remove 'position' tag for all the file names depending
%on if calling the stacked versions or not 

load([datadir,filesep, 'RatioData_raw.mat'],'imRatio_raw');
load([datadir,filesep, 'Bleach_raw.mat']);
timepts=1: length(imRatio_raw);

 corr=feval(fitpara,timepts); %use this one
 corr_norm=corr./median(corr); %use this one
 normfact=nanmedian(bleach_raw); %use this one
 



 imRatio= cell(1,size(imRatio_raw,2));
 
for frameNum=1:length(imRatio_raw)
   
    imRatio{frameNum}=imRatio_raw{frameNum}./(normfact*corr_norm(frameNum));
 colorRange = [0.7 1.3]; 
 
 colorRange3 = [0.8 1.2]; 
  %  imRatio{frameNum}=tempRATIO_corr;
   % added this in to make some max intensity projections 
   %tempRATIO_corr(tempRATIO_corr <1.1) = 0; 
    tempRATIOforstack=ratio2RGB( imRatio{frameNum},colorRange);%Cdc42
    
   
  imwrite(tempRATIOforstack,[datadir,filesep,  'Rho-FRET','.tif'],'WriteMode','append','Compression','none');
 %  imwrite(tempRATIOforstack,[datadir,filesep, 'Rho-FRET','max intensity proj2','.tif'],'WriteMode','append','Compression','none');
   

end



load([datadir,filesep, 'RatioData_raw.mat'],'maskFinal','cellCoors','imFRETOutline');
save([datadir,filesep, 'RatioData.mat'],'-v7.3','maskFinal','cellCoors','imRatio','imFRETOutline');
end    