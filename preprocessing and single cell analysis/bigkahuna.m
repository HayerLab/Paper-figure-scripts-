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

root= 'F:\Arnold 240508 single frame immunos T2\t567'; 

% if you want to create Tiff Stacks 
makeTiffStacks = 1; 

%channels = {'cyto'; 'CAAX';'mRuby'};
channels = {'hoechst';'pERM'; 'mRuby'; 'cyto'};
%channels = {'CFP'}; %'mRuby'; 'blank'}; 

%channel = ["cyto"; "CAAX"; "mRuby"];
channel = [ "hoechst"; "pERM"; "mRuby"; "cyto"];

% channel = ["CFP"]; % "mRuby";"blank"]; 
%channel1= {'cyto' 'CAAX' 'mRuby'};
channel1= {'hoechst'  'pERM' 'mRuby' 'cyto'};

%channel1= {'CFP'}; % 'mRuby'}; 


%  channels = {'mCit'}; 
%  channel = ["mCit"]; 
%  channel1= {'mCit'}; 

threshold = 4; % for threshold based segmentation
frames = 1; 

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

for row = 3
    for col = 1
        for site =1:num_sites
            
            
          for timept = 1:1
                
            
            shot = [num2str(row),'_',num2str(col),'_',num2str(site)];
            
        if size (channels, 1) == 1
             [mCit]= nd2read_hayer_1chan(bgfilepath,finfo, counter,1); 
             imwrite(mCit,[bgdir,filesep,shot,'_mCit_',num2str(timept),'.tif'],'TIFF','Compression','None');
        end 
            
       if size(channels, 1) ==2
         [FRET, CFP]= nd2read_hayer_2chan(bgfilepath,finfo, counter,1); 
         imwrite(FRET,[bgdir,filesep,shot,'_FRET_',num2str(timept),'.tif'],'TIFF','Compression','None');
         imwrite(CFP,[bgdir,filesep,shot,'_CFP_',num2str(timept),'.tif'],'TIFF','Compression','None');
       end 
       
       if size (channels,1) == 3
              [mRuby, cyto,CAAX]= nd2read_hayer(bgfilepath,finfo, counter,1); 
              imwrite(cyto,[bgdir,filesep,shot,'_cyto_',num2str(timept),'.tif'],'TIFF','Compression','None');
        imwrite(CAAX,[bgdir,filesep,shot,'_CAAX_',num2str(timept),'.tif'],'TIFF','Compression','None');
        imwrite(mRuby,[bgdir,filesep,shot,'_mRuby_',num2str(timept),'.tif'],'TIFF','Compression','None');
  
       end 
            
       if size (channels,1) == 4
              [FRET,CFP,mRuby,cyto]= nd2read_hayer(bgfilepath,finfo, counter,1); 
        imwrite(CFP,[bgdir,filesep,shot,'_CFP_',num2str(timept),'.tif'],'TIFF','Compression','None');
        imwrite(FRET,[bgdir,filesep,shot,'_FRET_',num2str(timept),'.tif'],'TIFF','Compression','None');
        imwrite(mRuby,[bgdir,filesep,shot,'_mRuby_',num2str(timept),'.tif'],'TIFF','Compression','None');
         imwrite(cyto,[bgdir,filesep,shot,'_cyto_',num2str(timept),'.tif'],'TIFF','Compression','None');
            
      if size (channels,1) == 6
             [mRuby, blank, cyto, blank2, FRET, CFP ] = nd2read_hayer_6chan(filepath,finfo,counter,1);
    imwrite(CFP,[bgdir,filesep,shot,'_CFP_',num2str(timept),'.tif'],'TIFF','Compression','None');
        imwrite(FRET,[bgdir,filesep,shot,'_FRET_',num2str(timept),'.tif'],'TIFF','Compression','None');
        imwrite(mRuby,[bgdir,filesep,shot,'_mRuby_',num2str(timept),'.tif'],'TIFF','Compression','None');
        imwrite(cyto,[bgdir,filesep,shot,'_cyto_',num2str(timept),'.tif'],'TIFF','Compression','None');
      end 
      %
      
      
    
      
      counter = counter+1;
      
            end   
         end 
        end 
    end 
end 

%% files
% num_sites = 0; 
 for i = 1:(size(ND2files,1))
filepath = [ND2dir, filesep, ND2files{i}];
%filepath = [ND2dir, filesep, ND2files{i}]; 
counter =0; 


finfo = nd2finfo(filepath);
num_sites = finfo.img_seq_count/frames;  

for timept = 1:frames
for row = i
    for col = 1
        for site =1:num_sites
       
            
            shot = [num2str(row),'_',num2str(col),'_',num2str(site)];
            
       if size (channels, 1) == 1
             [CFP]= nd2read_hayer_1chan(filepath,finfo, counter,1); 
             imwrite(CFP,[rawdir,filesep,shot,'_CFP_',num2str(timept),'.tif'],'TIFF','Compression','None');
        end 
            
       if size(channels, 1) ==2
         [FRET, CFP]= nd2read_hayer_2chan(filepath,finfo, counter,1); 
         imwrite(FRET,[rawdir,filesep,shot,'_FRET_',num2str(timept),'.tif'],'TIFF','Compression','None');
         imwrite(CFP,[rawdir,filesep,shot,'_CFP_',num2str(timept),'.tif'],'TIFF','Compression','None');
       end 
       
       if size (channels,1) == 3
              [mRuby, cyto,CAAX]= nd2read_hayer(filepath,finfo, counter,1); 
        imwrite(cyto,[rawdir,filesep,shot,'_cyto_',num2str(timept),'.tif'],'TIFF','Compression','None');
        imwrite(CAAX,[rawdir,filesep,shot,'_CAAX_',num2str(timept),'.tif'],'TIFF','Compression','None');
        imwrite(mRuby,[rawdir,filesep,shot,'_mRuby_',num2str(timept),'.tif'],'TIFF','Compression','None');
       end 
            
       if size (channels,1) == 4
              [hoechst, pERM, mRuby, cyto]= nd2read_hayer(filepath,finfo, counter,1); 
           imwrite(hoechst,[rawdir,filesep,shot,'_hoechst_',num2str(timept),'.tif'],'TIFF','Compression','None');
        imwrite(pERM,[rawdir,filesep,shot,'_pERM_',num2str(timept),'.tif'],'TIFF','Compression','None');
        imwrite(mRuby,[rawdir,filesep,shot,'_mRuby_',num2str(timept),'.tif'],'TIFF','Compression','None');
        imwrite(cyto,[rawdir,filesep,shot,'_cyto_',num2str(timept),'.tif'],'TIFF','Compression','None');
       end 
            
      if size (channels,1) == 6
             [mRuby, blank, cyto, blank2, FRET, CFP ] = nd2read_hayer_6chan(filepath,finfo,counter,1);
     imwrite(mRuby,[rawdir,filesep,shot,'_mRuby_',num2str(timept),'.tif'],'TIFF','Compression','None');
     imwrite(cyto,[rawdir,filesep,shot,'_cyto_',num2str(timept),'.tif'],'TIFF','Compression','None');
     imwrite(FRET,[rawdir,filesep,shot,'_FRET_',num2str(timept),'.tif'],'TIFF','Compression','None');
       imwrite(CFP,[rawdir,filesep,shot,'_CFP_',num2str(timept),'.tif'],'TIFF','Compression','None');
      end 
      %
      
      
    
      
      counter = counter+1;
      
            end   
        end 
    end 
end 
 end 
 


%% tiff stacks 
if makeTiffStacks == 1
    
    
 for row=1:(size(ND2files,1))
    
    for col = 1
        for site =1
            
          
           
writeTiffStacks(root, channel, row, col, site, frames); 
        end 
    end 
end   
end 

 %% generate background images 


for chan=1:numel(channel1')
    calculate_bg_img_rm_blobs(bgdir,channel1{chan});
    disp(channel1(chan));
    disp(num2str(chan));
end
%% CFP/FRET channel alignment 

if size (channels, 1)> 1 
CFPStack=[];
FRETStack=[];

k=0;
for row=3:4 % (size(ND2files,1)-1)
    
   
    for col=1
        for site=1:32
            k=k+1;
            shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
            CFP_temp=imread([rawdir,filesep,shot,'_CFP_4.tif']);
            FRET_temp=imread([rawdir,filesep,shot,'_FRET_4.tif']);
        
            CFPStack(:,:,k)=CFP_temp;
            FRETStack(:,:,k)=FRET_temp;
        
            
        end
    end
end


CFP_AV=uint16(mean(CFPStack,3));
FRET_AV=uint16(mean(FRETStack,3));

imwrite(CFP_AV,[bgdir,filesep,'AVG_rawdata_mRuby.tif'],'TIFF','Compression','None');
imwrite(FRET_AV,[bgdir,filesep,'AVG_rawdata_cyto.tif'],'TIFF','Compression','None');

alignStack(:,:,2)=imread([bgdir,filesep,'AVG_rawdata_cyto.tif']);
alignStack(:,:,1)=imread([bgdir,filesep,'AVG_rawdata_mRuby.tif']);
[pX,pY,dxMat1,dyMat1]=dualviewComputeAlignmentFromGridImages(alignStack);
figure;
subplot(1,2,1); imagesc(dxMat1); colorbar
subplot(1,2,2); imagesc(dyMat1); colorbar

save([bgdir,filesep,'alignment parameters pX pY_day2.mat'],'pX','pY','dxMat1','dyMat1');
end 
%% generate raw FRET Data baby 

%load([bgdir,filesep,'alignment parameters pX pY.mat']);
threshold = 4 ; 
k=0;
for row=3  % (size(ND2files,1)-1)
%     
    for col=1
        for site=17:32%num_sites

%             if row ==1 && site ==8 
%     continue; 
% end 
            
            k=k+1;
            
            position{k}=[num2str(row),'_',num2str(col),'_',num2str(site)];
        end
    end
end 

%pick which one you need 
for k=1:length(position)
  %getFRETDataHCS_1chan(position{k},bgdir,rawdir,datadir,threshold); 
%getFRETDataHCS(position{k},bgdir,rawdir,datadir,1.3); 
getFRETDataHCS_2ChanFRET(position{k},bgdir,rawdir,datadir,1.17); 
  %getFRETDataHCS_3chan(position{k},bgdir,rawdir,datadir,5)
 %getFRETDataHCS_4chan(position{k},bgdir,rawdir,datadir,threshold)
end
disp('done!');
%clc; clear; 

%% bleaching correction

bleach_raw_all=[];
for k=1:length(position)
    load([datadir,filesep,position{k},'_Bleach_raw.mat']);
    bleach_raw_all(:,k)=bleach_raw/nanmedian(bleach_raw);
end
bleach_mean=nanmedian(bleach_raw_all(1:75, :),2);
plot(1:length(bleach_mean),bleach_raw_all(1:75, :));hold on
plot(1:length(bleach_mean),bleach_mean,'linewidth',3);

title('Bleaching Curve'); 
xlabel('TimePoint');
ylabel('Mean Intensity');

timepts=1:length(bleach_mean);
fitpara=fit(timepts',bleach_mean,'exp2');
corr=feval(fitpara,1:length(bleach_mean));
plot(1:length(bleach_mean),corr,'g');
axis([0 75 0.7 1.5]);

save([datadir,filesep,'bleachingcurve.mat'],'fitpara');
pause; 
bleach_raw_all=[]; for k=1:length(position)
    load([datadir,filesep,position{k},'_Bleach_raw.mat']);
    bleach_raw_all(:,k)=bleach_raw_mRuby/nanmedian(bleach_raw_mRuby);
end
bleach_mean=nanmedian(bleach_raw_all (1:75, :),2);
plot(1:length(bleach_mean),bleach_raw_all(1:75, :) );hold on
plot(1:length(bleach_mean),bleach_mean,'linewidth',3);

title('Bleaching Curve'); xlabel('TimePoint'); ylabel('Mean Intensity');

timepts=1:length(bleach_mean);
fitpara_mRuby=fit(timepts',bleach_mean,'exp2');
corr=feval(fitpara_mRuby,1:length(bleach_mean));
plot(1:length(bleach_mean),corr,'g'); axis([0 75 0 3]);

 save([datadir,filesep,'bleachingcurve_mRuby.mat'],'fitpara_mRuby');
 
 pause; 
 bleach_raw_all=[]; for k=1:length(position)
    load([datadir,filesep,position{k},'_Bleach_raw.mat']);
    bleach_raw_all(:,k)=bleach_raw_cyto/nanmedian(bleach_raw_cyto);
end
bleach_mean=nanmedian(bleach_raw_all (1:75, :),2);
plot(1:length(bleach_mean),bleach_raw_all(1:75, :) );hold on
plot(1:length(bleach_mean),bleach_mean,'linewidth',3);

title('Bleaching Curve'); xlabel('TimePoint'); ylabel('Mean Intensity');

timepts=1:length(bleach_mean);
fitpara_cyto=fit(timepts',bleach_mean,'exp2');
corr=feval(fitpara_cyto,1:length(bleach_mean));
plot(1:length(bleach_mean),corr,'g'); axis([0 75 0 3]);

 save([datadir,filesep,'bleachingcurve_cyto.mat'],'fitpara_cyto');

%% corrected FRET data 
for k=1:length(position)
   
   correctBleachingExp2_FRET(position{k},fitpara, datadir);
   % correctBleachingExp2(position{k},fitpara, datadir, fitpara_mRuby);
     %correctBleachingExp2_cyto_ratio(position{k}, datadir, fitpara_mRuby, fitpara_cyto);
end 

