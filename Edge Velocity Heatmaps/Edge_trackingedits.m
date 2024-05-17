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


%need to add "Statistics and Machine Learning Toolbox" in MATLABS add-on
%section. 

%saves following figures: labelled cell mask with coordinates, area change,
%and edge velocity heat maps 
clc; clear; 
cells = [1:1:20];   

for place=1:size(cells,2)
cells = [1:1:60];   
     


for depth = 1:5 % number of edge depths you wish to analyze 


 %for 20x no binning, or 40x 2x2 bin etc. 0.325 um/px
depths = [3,6,10,15,25]; 



 
root='F:\example_dataset_240522';
rawdir=([root,filesep,'cropped',filesep, strcat(num2str(cells(1,place))),filesep,'output']); %
datadir=([rawdir,filesep,'edge_vels', filesep,  strcat('edge vel mapping_',num2str(depths(1, depth)))]); 
if ~exist(datadir)
    mkdir(datadir)
end 


% Parameters for cell edge parametrization
nFretWindows=180;                                   % Number of windows to use for FRET vs edge correlation measurements
edgeOversamplingParam=5;                            % How many times more points should the edge of the cell be tracked at ("subwindow")
nPointsParam=nFretWindows*edgeOversamplingParam;    % Number of points to track on the cell edge
pdSmoothing=10;                                   % Used with imclose to make the selection of points for tracking less dependent on noise or wrinkles in the cell edge

edgeDepthDist=depths(1, depth);   % Number of pixels deep for the windows for computing FRET values.


binning=1;            %only change if binning is changed while using same objective magnification!!




% Load previously determined sequence of masks


load([rawdir,filesep,'RatioData_raw.mat']);
load([rawdir,filesep,'ezrin_data', filesep,'CytoRatioData.mat']);


%% %% delete unwanted frames (optional!!) 
% ie if cell touches another cell at the end of the movie and you wish to
% exclude it from edge tracking analysis  

flashFrame=141:150; 
maskFinal(flashFrame)=[];
imCAAXOutline(flashFrame) = []; 
ezrin_ratio(flashFrame)=[];
membrane_cyto_ratio(flashFrame)=[];
cellCoors(flashFrame)=[]; 
% % im_mRuby(flashFrame) = []; 


%% determines cell trajectories. 
% Will create multiple if there are multiple masked cells per crop 
% in your loadeddata 

%first checks if any trajectories exist, and then skips through any empty
%cells at the beginning 

%begins edge tracking at the first non-empty cell in the array cellCoors


empty_count=0;
if  ~isempty(cellCoors)

    for i=1:size(cellCoors,2)
       temp = cellCoors(1,i);
       contents=temp{1,1};
           if isempty(contents)
             empty_count=empty_count+1;
             continue; 

           else   
            traj=ultTrackAnnSearch(cellCoors(1,empty_count+1:end),'pairrule','fwdbckmtch','maxdisp',100,'minlength',5,'verbose',false); 
         end 
    end 
else
    traj={};
end
fprintf('%i trajectories.\n',length(traj));






%% Manually select cell for analysis 
selectedCell=1; %input which trajectory
isConnect = false; %true if you had to connect broken trajectories

thisTraj=traj{selectedCell};
centroid_coordinates=zeros(2,size(imCAAXOutline,2));
cell_area=0; 
if ~(isConnect) 
    start =thisTraj(1,5);
else
    
    start=thisTraj(1,4);

    for z = start:start+size(thisTraj,1) -1
        thisTraj(z,4)=z;
    end 

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
            

                % Parametrize cell edge and compute protrusion values
                if index ==1
                    [edgeCoors{index}, edgeCoorsSmoothed{index}]=parametrizeCellEdge(thisMask,nPointsParam,round(pdSmoothing/binning));
                else
                    [edgeCoors{index},edgeCoorsSmoothed{index}]=parametrizeCellEdge(thisMask,nPointsParam,round(pdSmoothing/binning),edgeCoors{index-1});
                   
                    protvals(:,index-1)=vect(computeProtrusionValues(edgeCoors{index-1},edgeCoorsSmoothed{index-1},edgeCoors{index}));
                end
                windowCoors{index}=edgeCoors{index}((edgeOversamplingParam+1)/2:edgeOversamplingParam:end,:);

            % Compute edge regions and local fret ratio values
            labelMask{index}=getWindowLabelMap(thisMask,windowCoors{index},(edgeDepthDist/binning));
         
            for k=1:size(windowCoors{index},1)
                
                %comment out as needed for specific data 
                 fretvals(k,index)=mean(imRatio_raw{index+empty_count}(labelMask{index}==k));
                 myosin(k,index)=mean(im_mRuby{index+empty_count}(labelMask{index}==k));
                 ezrin(k,index)=mean(ezrin_ratio{index+empty_count}(labelMask{index}==k));
                 membrane_cyto(k,index)=mean(membrane_cyto_ratio{index+empty_count}(labelMask{index}==k));
          
  
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
 
 %% filtered protusionvalues
 
   % filtered protrusion values
protvalsWindowF=ndnanfilter(protvalsWindow,fspecial('disk',2),'replicate');
fretvalsF=ndnanfilter(fretvals,fspecial('disk',2),'replicate'); 
myosinF=ndnanfilter(myosin,fspecial('disk',2),'replicate');
ezrinF=ndnanfilter(ezrin,fspecial('disk',1),'replicate');
 membranecytoF=ndnanfilter(membrane_cyto,fspecial('disk',1),'replicate');

 
%% Plot maps - with thresholds
load('desired colormap.mat'); %any existing colormap can be used for visualization
load('CMAP_blue_grey_yellow.mat'); %default one used in paper 
close all;


f1=figure; 


hold on;


%following code creates 2 by 2 subplot where edge velocity, FRET activity, 
% and other protein localization can be observed. Must specify which data
% you want to plot where 

ax1=subplot(2,2,1);imagesc(protvalsWindow,[-13,13]);title('Edge Velocity');
colormap(ax1,BCWOR); 
%25 s intervals 
   xticks([0 24 48 71 95 120 144])
   xticklabels({'0','10','20','30','40','50' '60'});

ax2=subplot(2,2,2);imagesc(ezrinF,[-1 3] );title('Ezrin');
xlim([0 150]); 
 xticks([24 48 71 95 120 144])
  xticklabels({'0','10','20','30','40','50'});


ax3=subplot(2,2,3);imagesc(protvalsWindowF, [-13 13]); title ('Filtered');
colormap(ax3,BCWOR);
xlim([0 150]); 
xticks([ 24 48 71 95 120 144])
xticklabels({'0','10','20','30','40','50'});

 ax4=subplot(2,2,4);imagesc(myosinF,[0 3]); title ('Myosin localization');
 xlim([0 150]); 
 xticks([0 24 48 71 95 120 144])
 xticklabels({'0','10','20','30','40','50' '60'});

 hold off; 
saveas(f1,strcat(datadir,'\','edge_velocity_map.png'))
saveas(f1,strcat(datadir,'\','edge_velocity_mapM.fig'))


% save all new data into mat file 

save(strcat(datadir,'\','Protrusion and FRET Values.mat'),'protvalsWindow','protvalsWindowF','distance', 'ezrin', 'ezrinF', 'membrane_cyto','membranecytoF') ; % 'avg_cell_area'%'myosin','myosinF'); % ;


%close all; clc;

end 
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



