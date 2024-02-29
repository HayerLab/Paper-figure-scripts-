%% Test protrusion/FRET measurements
% Required functions: 
%       functions located folder "trackingcode"  
%       parametrizeCellEdge.m
%       computeProtrusionValues.m
%       getWindowLabelMap.m       
%       mapSmoothedEdgeCoorsToCellBoundary.m
%       matrixCircularlyPermuteRows.m
%       computeNormalVectorsFromParametrizedCellEdge.m
%       CMAP_blue_grey_yellow.mat
%       ndnanfilter.m
%       vect.m

%saves following figures: labelled cell mask with coordinates, area change,
%and edge velocity heat maps 
clc; clear; 
cells = [1 2 3];   

for place=1:size(cells,2)
cells = [1 2 3];    


%for zwster = 1:6

%clear; close all; clc;
 %for 20x no binning, or 40x 2x2 bin etc. 0.325 um/px
depths = [3]; % 6,10,15,20,25]; 

% for 60x, 2x2 binning
%depths = [5,9,15,23,30,38]; 
 
root='I:\Nada\sephsversion';
rawdir=([root,filesep,'cropped',filesep, strcat(num2str(cells(1,place))),filesep,'output']); %
datadir=([rawdir,filesep,'edge_vels', filesep,  strcat('edge vel mapping_',num2str(3))]); %'cropped', filesep,  strcat(num2str(cells(1,place))), filesep,
if ~exist(datadir)
    mkdir(datadir)
end 


% Parameters for cell edge parametrization
nFretWindows=180;                                   % Number of windows to use for FRET vs edge correlation measurements
edgeOversamplingParam=5;                            % How many times more points should the edge of the cell be tracked at ("subwindow")
nPointsParam=nFretWindows*edgeOversamplingParam;    % Number of points to track on the cell edge
pdSmoothing=10;%was 10                                      % Used with imclose to make the selection of points for tracking less dependent on noise or wrinkles in the cell edge

edgeDepthDist= 3;  %depths(1, zwster);   % Number of pixels deep for the windows for computing FRET values.

startFrame=1;
endFrame=50;

binning=1;            %only change if binning is changed while using same objective magnification!!

%colormap('parula');


% Load previously determined sequence of masks


load([rawdir,filesep,'RatioData_raw.mat']);
%load([rawdir,filesep,'CytoRatioData.mat']);

% perform tracking of mask centroid
% min(cellfun(@(x) size(x,1),cellCoors))    % changed it from min to max ? 
%% %% delete flash glitch frames that mess FRET data - if needed 
% % % % 
% flashFrame=1:29; 
% maskFinal(flashFrame)=[];
% imFRETOutline(flashFrame)=[];
% imRatio_raw(flashFrame)=[];
% cellCoors(flashFrame)=[]; 
% % im_mRuby(flashFrame) = []; 
% % % 
% flashFrame=114; 
% maskFinal(flashFrame)=[];
% imFRETOutline(flashFrame)=[];
% imRatio(flashFrame)=[];
% cellCoors(flashFrame)=[]; 
% % % im_mRuby(flashFrame) = []; 
% % % % 
% flashFrame=53; 
% maskFinal(flashFrame)=[];
% imFRETOutline(flashFrame)=[];
% imRatio(flashFrame)=[];
% cellCoors(flashFrame)=[]; 
% % % im_mRuby(flashFrame)=[]; 
% 
% flashFrame=1:19; 
% maskFinal(flashFrame)=[];
% imFRETOutline(flashFrame)=[];
% imRatio(flashFrame)=[];
% cellCoors(flashFrame)=[]; 

%%

%first checks if any trajectories exist, and then skips through any empty
%cells at the beginning 

%added iterations that deal with empty frames before cells are dropped,
%begins edge tracking at the first non-empty cell in the array cellCoors


empty_count=0;
if  ~isempty(cellCoors)

    for i=1:size(cellCoors,2)
       temp = cellCoors(1,i);
       contents=temp{1,1};
           if isempty(contents)
             empty_count=empty_count+1;
             continue; 

           else  %here need to add in 'mem' arg with number of frames you want it to skip, otherwise you need to take it out 
            traj=ultTrackAnnSearch(cellCoors(1,empty_count+1:end),'pairrule','fwdbckmtch','maxdisp',100,'minlength',5,'verbose',false); %wanting to add in a min length? - see if it changes anything 
         end %min(cellfun(@(x) size(x,1),cellCoors))>0
    end 
else
    traj={};
end
fprintf('%i trajectories.\n',length(traj));






%% Manually select cell for analysis 
selectedCell=1; %input which trajectory
isConnect = false; %true if you had to connect broken trajectories
%deleteFrame =46; % manually input frame that needs to be deleted
thisTraj=traj{selectedCell};
centroid_coordinates=zeros(2,size(imFRETOutline,2));
cell_area=0; 
if ~(isConnect) %modifying this Traj code if you needed to correct broken trajectores
    start =thisTraj(1,5);
else
    
    start=thisTraj(1,4);

    for z = start:start+size(thisTraj,1) -1
        thisTraj(z,4)=z;
    end 
%     maskFinal(deleteFrame)=[]; 
%     imFRETOutline(deleteFrame)=[];
%     imRatio(deleteFrame)=[];
end

index = 1;
for imnum=start:start+size(thisTraj,1) -1
   
        % Make a one cell mask
        objects=regionprops(maskFinal{empty_count+imnum},'PixelIdxList','PixelList','Centroid','BoundingBox', 'Area');

        cellCent=round(thisTraj(find(thisTraj(:,end)==imnum,1),1:2));
        centroid_coordinates(:,imnum) = cellCent'; 
       
       if size(objects,1) == 1
        cell_area =cell_area +objects.Area; 
      end 

      
       if size(objects,1) > 1
        
           cell_areas = [objects.Area]; 
           cell_area = cell_area+max(cell_areas); 
           
       end 
           

            cellNum=find(round(arrayfun(@(x) x.Centroid(1),objects))==cellCent(1) & round(arrayfun(@(x) x.Centroid(2),objects))==cellCent(2));
            thisMask_raw=false(size(maskFinal{imnum}));
            thisMask_raw(objects(cellNum).PixelIdxList)=true;        
            thisMask_raw=false(size(maskFinal{imnum})); % sets a frame of the size of the mask to zero
            thisMask_raw(objects(cellNum).PixelIdxList)=true; % sets values of object 2 to 1

            thisMask = thisMask_raw;
            % Smooth mask - move this into the getCellCyto function, so you
            % can see what the cell actually looks like 
%             thisMask=imerode(thisMask_raw,strel('disk',4));
%             thisMask=imdilate(thisMask,strel('disk',4));
%             thisMask=bwareaopen(thisMask,300);% Clean up smoothing artifacts (remove objects smaller than 500 pixels). 

                % Parametrize cell edge and compute protrusion values
                if index ==1
                    [edgeCoors{index}, edgeCoorsSmoothed{index}]=parametrizeCellEdge(thisMask,nPointsParam,round(pdSmoothing/binning));

                 %   plot(edgeCoors{1,index}(:,1),edgeCoors{1,index}(:,2)); hold on; 
                  %  plot(edgeCoorsSmoothed{1,index}(:,1),edgeCoorsSmoothed{1,index}(:,2))
                else
                    [edgeCoors{index},edgeCoorsSmoothed{index}]=parametrizeCellEdge(thisMask,nPointsParam,round(pdSmoothing/binning),edgeCoors{index-1});
                    %windowCoors{frameCount}=edgeCoors{frameCount}((edgeOversamplingParam+1)/2:edgeOversamplingParam:end,:);
                    protvals(:,index-1)=vect(computeProtrusionValues(edgeCoors{index-1},edgeCoorsSmoothed{index-1},edgeCoors{index}));
                end
                windowCoors{index}=edgeCoors{index}((edgeOversamplingParam+1)/2:edgeOversamplingParam:end,:);

            % Compute edge regions and local fret ratio values
            labelMask{index}=getWindowLabelMap(thisMask,windowCoors{index},(edgeDepthDist/binning));
         
            for k=1:size(windowCoors{index},1)
                
                 fretvals(k,index)=mean(im_YFP_raw{index+empty_count}(labelMask{index}==k));
               %myosin(k,index)=mean(im_mRuby{index+empty_count}(labelMask{index}==k));
               %  cyto(k,index)=mean(ezrin_ratio{index+empty_count}(labelMask{index}==k));
          
  
            end
           MaskThisCell{index}=thisMask;  
   
           disp(num2str(imnum));

           index=index+1;             
end

protvalsWindow=zeros(nFretWindows,size(protvals,2));

 for k=1:edgeOversamplingParam
     
     protvalsWindow=protvalsWindow+protvals(k:edgeOversamplingParam:end,:);
     
 end
 %% calculate total distance travelled by cell centroid and average cell area 
 distance = 0; 
 for i = 2:size(centroid_coordinates,2)
     dx = centroid_coordinates(1,i)-centroid_coordinates(1,i-1); 
     dy = centroid_coordinates(2,i)-centroid_coordinates(2,i-1); 
     
  distance = distance+ sqrt((dx^2)+(dy^2));       
 end 
 
 % make it per minute so can compare different length movies 
 distance = distance/(size(imFRETOutline,2)*(2/3))
 avg_cell_area = cell_area/size(im_YFP_raw,2)
 
 %% filtered protusionvalues
 
   % filtered protrusion values
protvalsWindowF=ndnanfilter(protvalsWindow,fspecial('disk',2),'replicate');
fretvalsF=ndnanfilter(fretvals,fspecial('disk',2),'replicate');
% myosinF=ndnanfilter(myosin,fspecial('disk',2),'replicate');
% cytoF=ndnanfilter(cyto,fspecial('disk',2),'replicate');



protvalsrangeF=[round(prctile(protvalsWindowF(:),1),1),round(prctile(protvalsWindowF(:),99),1)];
%fretvalsrangeF=[round(prctile(fretvalsF(:),1),1),round(prctile(fretvalsF(:),99),1)];
%myosinrangeF=[round(prctile(myosinF(:),1),1),round(prctile(myosinF(:),99),1)];


%% display location of coordinates overlayed on cell mask images - optional for visualization
%  % creates gif file for reference
    for mapper =1:size(thisTraj,1)
     
        window = windowCoors{1,mapper}; 
      if ~(isConnect)
      image = (imFRETOutline{1,thisTraj(mapper,5)+empty_count}); %slightly convoluted but lines up YFP raw image to the correct windowcoors
      else
           image = imread(imFRETOutline{1,thisTraj(mapper,4)+empty_count});
      end 
      h = figure('visible','off');
      hold on; 
      axis ij; 
     imagesc(image); 
   
     %for num = 1:15:size(window,1)
      for num = 100:2:120
         image(window(num,1), window(num,2),1) = 255; 
         image(window(num,1), window(num,2),2) = 0; 
         image(window(num,1), window(num,2),3) = 0; 
         text(window(num,2),window(num,1),num2str(num), 'Color','r', 'FontSize', 8);
      end 
     
     frame = getframe(h);
    im=frame2im(frame);
    [imind, cm] = rgb2ind(im,256);
     
         if mapper ==1
             imwrite(imind,cm,strcat(datadir,'\','Coordinate_Windows'),'gif','Loopcount', inf);
         else
                   imwrite(imind,cm,strcat(datadir,'\','Coordinate_Windows'),'gif','WriteMode','append');
         end 
      imwrite(image,[datadir,filesep,'Outline_label.tif'],'WriteMode','append','Compression','none');

  end 
% % %   

 
%% Plot maps - with thresholds
close all;

protthresh=5;
retthresh=-5;
f1=figure; 


hold on;
% load custom parula colour map for protrusion/retraction visualization
%load('F:\Seph\code\supporting_functions\trackingcode\CMAP_blue_grey_yellow.mat'); 
load('CMAP_blue_grey_yellow.mat');
% do this if imaging done at different time interval than 25s
cmap_15s =cmap; 
% cmap_15s(35:39,:)=[]; %adjusting "grey range" depending on thresholds
protvalrange=[round(prctile(protvalsWindow(:),1),1),round(prctile(protvalsWindow(:),99),1)];
%fretvalsrange=[round(prctile(fretvals(:),1),1),round(prctile(fretvals(:),99),1)];
 %myosinrange=[round(prctile(myosin(:),1),1),round(prctile(myosin(:),99),1)];

ax1=subplot(2,2,3);imagesc(protvalsWindow,[-13,13]);title('Edge Velocity');
colormap(ax1,cmap);
%15s intervals 
%  xticks([40 80 120 160 200]); 
% xticklabels({'0','10','20','30','40','50'}); 
%25 s intervals 
 xticks([0 25 50])
 xticklabels({'0','25','50'});

ax2=subplot(2,2,1);imagesc(fretvalsF,[400 430] );title('YFP');
xlim([0 50]); 
 %ax2 =subplot(2,2,2);imagesc(protvalsWindowF,[-13,13]);title('Edge Velocity');
% colormap(ax2,cmap);
%  rectangle('Position',[84,120,41,15],'LineWidth',2); 
%  rectangle('Position',[62,60,41,15],'LineWidth',2); 
 %15s intervals 
%  xticks([40 80 120 160 200]); 
% xticklabels({'0','10','20','30','40','50'}); 
%25 s intervals 
 xticks([0 25 50])
  xticklabels({'0','25','50'});



protvalsWindowHigh=protvalsWindow>protthresh;

ax3=subplot(2,2,4);imagesc(protvalsWindowF, [-13 13]); title ('Filtered');
colormap(ax3,cmap);
xlim([0 50]); 
%15s intervals 
%  xticks([40 80 120 160 200]); 
% xticklabels({'0','10','20','30','40','50'}); 
%25 s intervals 
 xticks([0 25 50])
  xticklabels({'0','25','50'});
%protvalsWindowFHigh=protvalsWindowF>protthresh;

%  ax4=subplot(2,2,2);imagesc(fretvalsF,[0.3 1.7]); title ('ezxrin');
%  %15s intervals 
% %  xticks([40 80 120 160 200]); 
% % xticklabels({'0','10','20','30','40','50'}); 
%  %25 s intervals 
%   xticks([0 24 48 71 95 120 144])
%   xticklabels({'0','10','20','30','40','50' '60'});
%
 hold off; 
saveas(f1,strcat(datadir,'\','edge_velocity_map.png'))
saveas(f1,strcat(datadir,'\','edge_velocity_mapM.fig'))

%save('protrusionpeaks.mat','protvalsWindowF','protvalsWindowFHigh')

% save all new data into mat file 

save(strcat(datadir,'\','Protrusion and FRET Values.mat'),'protvalsWindow','protvalsWindowF','distance','avg_cell_area','fretvals','fretvalsF') ; %'myosin','myosinF'); % 'cyto', 'cytoF');


%close all; clc;

end 
clear; clc; 

 
%% Polar plot of sum of all protrusions - all code after this is optional
close all;
protsum=sum(protvalsWindowF,2);
polarplot((1:2:360)*pi/180,protsum);rlim([-200 200]);

%% Plot maps. 
% close all;
% % edge velocity/RhoA maps
% load('E:\supporting_functions\trackingcode\CMAP_blue_grey_yellow.mat');
% %protvalrange=[round(prctile(protvalsWindow(:),1),1),round(prctile(protvalsWindow(:),99),1)];
% %fretvalsrange=[round(prctile(fretvals(:),1),1),round(prctile(fretvals(:),99),1)];
% %redvalsrange=[round(prctile(redvals(:),1),1),round(prctile(redvals(:),99),1)];
% ax1=subplot(2,2,1);imagesc(protvalsWindow, [-12 12]);title('protrusion');
% colormap(ax1,cmap);
% %ax2=subplot(2,2,2);imagesc(fretvals,fretvalsrange);title('FRET');
% %subplot(2,3,3);imagesc(redvals,redvalsrange);
% 
% % smoothing
% protvalsWindowF=ndnanfilter(protvalsWindow,fspecial('disk',2),'replicate');
% %fretvalsF=ndnanfilter(fretvals,fspecial('disk',2),'replicate');
% %redvalsF=ndnanfilter(redvals,fspecial('disk',2),'replicate');
% protvalsrangeF=[round(prctile(protvalsWindowF(:),1),1),round(prctile(protvalsWindowF(:),99),1)];
% %fretvalsrangeF=[round(prctile(fretvalsF(:),1),1),round(prctile(fretvalsF(:),99),1)];
% %redvalsrangeF=[round(prctile(redvalsF(:),1),1),round(prctile(redvalsF(:),99),1)];
% 
% ax3=subplot(2,2,3);imagesc(protvalsWindowF, protvalsrangeF);
% colormap(ax3,cmap); title('Protrusiuons-Smoothed');
% % subplot(2,3,5);imagesc(fretvalsF,fretvalsrangeF);
% %ax4=subplot(2,2,4);imagesc(fretvalsF,fretvalsrangeF);
% %subplot(2,3,6);imagesc(redvalsF,redvalsrangeF);
% 

% 
% 
% 
% 
% % Redvalsdiff:
% %redvalsFdiff=diff(redvalsF,[],2);
% %redvalsFdiffrange=[round(prctile(redvalsFdiff(:),1),1),round(prctile(redvalsFdiff(:),99),1)];
% 
% %subplot(2,3,6);imagesc(redvalsFdiff,redvalsFdiffrange);
%% plot area change of cell over time 
%x= linspace(1,1,size(thisTraj,1));
x = thisTraj(:,5);
y=thisTraj(:,3);
f2 =figure;
hold on; 
yyaxis left
title('Area Change vs. Mean Protrusion Length');
xlabel('frame number');
ylabel('cell area');


plot(x,y)
var_coeff = zeros(size(thisTraj, 1),1);
prot_mean=zeros(size(thisTraj, 1),1);
standard_dev=zeros(size(thisTraj, 1),1);
variance = zeros(size(thisTraj, 1),1);

for loop = 1:size(protvalsWindow, 2)
    standard_dev(loop+1,1) = std(protvalsWindowF(:,loop));
    prot_mean(loop+1,1) = mean(protvalsWindowF(:,loop));
    variance(loop+1,1)=var(protvalsWindowF(:,loop)); 
    var_coeff(loop+1,1) = (standard_dev(loop+1,1)/abs(prot_mean(loop+1,1))*100);
end 

yyaxis right
plot (x,prot_mean)
hold off; 
ylabel('mean protrusion length');

saveas(f2, strcat(datadir,'\','Area vs mean.png')); 

%getProbDensity(protvalsWindowF, 25); %here manually input which frame you want

    
%% Plot top 10% of velocity of protrusion, spreading area
for frame=1:(endFrame-start)
    SpeedTemp=protvalsWindowF(:,frame);
    SpeedTemp=sort(SpeedTemp,'descend');
    spreadingArea(frame,:)=sum(vect(MaskThisCell{frame}));
    protvalsFast(frame,:)=median(SpeedTemp(1:ceil(0.1*length(SpeedTemp))));
    %protvalsFast(frame,:)=median(protvalsWindow(:,frame)>prctile(protvalsWindow(:,frame),90),1);
end
yyaxis left;
plot(protvalsFast);

yyaxis right;
plot(spreadingArea);axis([0 90 15000 20000]);


%% quiver plot of edge velocities

for frameNum=1:(endFrame-start-1)
     close all; 
   % imagesc(MaskThisCell{frameNum}); hold on;
    
    image = imFRETOutline{1,frameNum};
           
      h = figure('visible','on');
      hold on; 
      axis ij; 
     imagesc(image);
    
     
    normV=computeNormalVectorsFromParametrizedCellEdge(edgeCoorsSmoothed{frameNum});
    protvals=(computeProtrusionValues(edgeCoors{frameNum},edgeCoorsSmoothed{frameNum},edgeCoors{frameNum+1}));
    %plot(edgeCoors{frameNum}(:,2),edgeCoors{frameNum}(:,1)); quiver(edgeCoors{frameNum}(:,2),edgeCoors{frameNum}(:,1),0.2*normV(:,2),0.2*normV(:,1));
    protvect=normV.*protvals;
    
    % filter protvect with moving average of windowSize 5
    windowSize=3;
    b=(1/windowSize)*ones(1,windowSize);
    a=1;
    
    protvectF(:,1)=filter(b,a,protvect(:,1));
    protvectF(:,2)=filter(b,a,protvect(:,2));
    
    plot(edgeCoors{frameNum}(:,2),edgeCoors{frameNum}(:,1)); quiver(edgeCoors{frameNum}(:,2),edgeCoors{frameNum}(:,1),0.2*protvect(:,2),0.2*protvect(:,1),'Color',[1 0 0]);
    plot(edgeCoors{frameNum}(:,2),edgeCoors{frameNum}(:,1)); quiver(edgeCoors{frameNum}(:,2),edgeCoors{frameNum}(:,1),0.2*protvectF(:,2),0.2*protvectF(:,1),'Color',[0 1 0]);
 
   
    
    
    plot(windowCoors{frameNum}(:,2),windowCoors{1}(:,1));hold on;
    plot(edgeCoors{frameNum}(:,2),edgeCoors{frameNum}(:,1));
    plot(edgeCoorsSmoothed{frameNum}(:,2),edgeCoorsSmoothed{frameNum}(:,1));
    qp.edgePos=windowCoors{frameNum};
    qp.edgeVel=6*(windowCoors{frameNum+1}-windowCoors{frameNum});
    quiver(qp.edgePos(:,2)',qp.edgePos(:,1)',qp.edgeVel(:,2)',qp.edgeVel(:,1)');hold off;%,'Color',[1 1 1]
    pause;  
end


%% overlay maps as greyscales into image
figure;
% fretred(:,:,1)=mat2gray(fretvalsF(:,1:end-1),fretvalsrangeF);
% fretred(:,:,2)=mat2gray(redvalsFdiff,redvalsFdiffrange);
% fretred(:,:,3)=zeros(size(fretred(:,:,1),1),size(fretred(:,:,1),2));
% subplot(1,3,1);imshow(fretred);axis square; title('red=fret, green=F-actin');

edgefret(:,:,1)=mat2gray(protvalsWindowF,protvalsrangeF);
%edgefret(:,:,2)=mat2gray(fretvalsF(:,1:end-1),fretvalsrangeF);
edgefret(:,:,2)=zeros(size(edgefret(:,:,1),1),size(edgefret(:,:,1),2));
imshow(edgefret);axis square; title('red=edge velocity, green=fret');

% edgered(:,:,1)=mat2gray(protvalsWindowF,protvalsrangeF);
% edgered(:,:,2)=mat2gray(redvalsFdiff,redvalsFdiffrange);
% edgered(:,:,3)=zeros(size(edgefret(:,:,1),1),size(edgefret(:,:,1),2));
% subplot(1,3,3);imshow(edgered);axis square; title('red=edge velocity, green=F-actin');

%% Scatter the top 20% of RhoB pixels over edge velocity
close all;

FREThigh=fretvalsF(:,1:end-1)>round(prctile(fretvalsF(:),80),1);
FRETlow=fretvalsF(:,1:end-1)<round(prctile(fretvalsF(:),20),1);
FRETlowhigh=single(or(FREThigh,FRETlow));

FREThigh=single(FREThigh);FREThigh(FREThigh==0)=NaN;
FRETlow=single(FRETlow);FRETlow(FRETlow==0)=NaN;

FRETlowhigh(FRETlowhigh==0)=NaN;
subplot(2,2,1);imagesc(fretvalsF(:,1:end-1),fretvalsrangeF);title('Cdc42 FRET/CFP');ylabel('peripheral windows');
subplot(2,2,2);imagesc((fretvalsF(:,1:end-1).*FRETlowhigh),fretvalsrangeF);title('Cdc42 FRET/CFP, top/bottom 20%');


subplot(2,2,3);imagesc(protvalsWindowF,protvalsrangeF);title('edge velocity');ylabel('peripheral windows');xlabel('frames (1/20s)');
subplot(2,2,4);imagesc((protvalsWindowF.*FRETlowhigh),protvalsrangeF);title('edge velocity, in top/bottom 20% Cdc42 FRET/CFP');xlabel('frames (1/20s)');

%subplot(2,2,4);dscatter(fretvalsF(FREThigh),protvalsWindowF(FREThigh));

figure; 

subplot(1,2,1);dscatter(vect(fretvalsF(:,1:end-1)),vect(protvalsWindowF));xlabel('Rac FRET/CFP');ylabel('edge velocity');
subplot(1,2,2);dscatter(vect(fretvalsF(:,1:end-1).*FRETlowhigh),vect(protvalsWindowF.*FRETlowhigh));xlabel('Rac FRET/CFP');ylabel('edge velocity');
nanmean(vect(protvalsWindowF.*FRETlow))
nanmean(vect(protvalsWindowF.*FREThigh))
lo=~isnan(vect(FRETlowhigh));
prot=vect(protvalsWindowF.*FRETlowhigh);prot=prot(lo);
fret=vect(fretvalsF(:,1:end-1).*FRETlowhigh);fret=fret(lo);
[R,P]=corrcoef(fret,prot)




%% Scatter the top 20% of edge velocity pixels over FRET
close all;

EDGEhigh=protvalsWindowF>round(prctile(protvalsWindowF(:),80),1);
EDGElow=protvalsWindowF<round(prctile(protvalsWindowF(:),20),1);
EDGElowhigh=single(or(EDGEhigh,EDGElow));

EDGEhigh=single(EDGEhigh);EDGEhigh(EDGEhigh==0)=NaN;
EDGElow=single(EDGElow);EDGElow(EDGElow==0)=NaN;

EDGElowhigh(EDGElowhigh==0)=NaN;

subplot(2,2,1);imagesc(protvalsWindowF,protvalsrangeF);title('edge velocity');ylabel('peripheral windows');
subplot(2,2,2);imagesc((protvalsWindowF.*EDGElowhigh),protvalsrangeF);title('edge velocity, top/bottom 20%');

subplot(2,2,3);imagesc(fretvalsF(:,1:end-1),fretvalsrangeF);title('Rac FRET/CFP');ylabel('peripheral windows');xlabel('frames (1/20s)');
subplot(2,2,4);imagesc((fretvalsF(:,1:end-1).*EDGElowhigh),fretvalsrangeF);title('Rac FRET/CFP, in top/bottom 20% edge velocity');xlabel('frames (1/20s)');

%subplot(2,2,4);dscatter(fretvalsF(FREThigh),protvalsWindowF(FREThigh));

figure; 

subplot(1,2,1);dscatter(vect(fretvalsF(:,1:end-1)),vect(protvalsWindowF));
subplot(1,2,2);dscatter(vect(fretvalsF(:,1:end-1).*EDGElowhigh),vect(protvalsWindowF.*EDGElowhigh));
nanmean(vect(fretvalsF(:,1:end-1).*EDGElow))
nanmean(vect(fretvalsF(:,1:end-1).*EDGEhigh))

lo=~isnan(vect(EDGElowhigh));
prot=vect(protvalsWindowF.*EDGElowhigh);prot=prot(lo);
fret=vect(fretvalsF(:,1:end-1).*EDGElowhigh);fret=fret(lo);
[R,P]=corrcoef(fret,prot)



