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

%NOTE: need to download Statistical Learning Toolbox  2020 by Dahua Lin

%saves following figures: labelled cell mask with coordinates, area change,
%and edge velocity heat maps 
clc; clear; 
cells = [1];   

for place=1:size(cells,2)
cells = [1];    



%clear; close all; clc;
 %for 20x no binning, or 40x 2x2 bin etc. 0.325 um/px
depths = [3]; % 6,10,15,20,25]; 

% for 60x, 2x2 binning
%depths = [5,9,15,23,30,38]; 
 
root='Y:\Baishali\GCAMP single cell edge tracking 17_01_2024\seph version';
rawdir=([root,filesep,'cropped',filesep,'cntrl',filesep, strcat(num2str(cells(1,place))),filesep,'output']); %
datadir=([rawdir,filesep,'edge_vels', filesep, 'edge vel mapping']); 
if ~exist(datadir)
    mkdir(datadir)
end 


% Parameters for cell edge parametrization
nFretWindows=180;                                   % Number of windows to use for FRET vs edge correlation measurements
edgeOversamplingParam=5;                            % How many times more points should the edge of the cell be tracked at ("subwindow")
nPointsParam=nFretWindows*edgeOversamplingParam;    % Number of points to track on the cell edge
pdSmoothing=10;%was 10                                      % Used with imclose to make the selection of points for tracking less dependent on noise or wrinkles in the cell edge

edgeDepthDist= 3;  %depths(1, zwster);   % Number of pixels deep for the windows for computing FRET values.


binning=1;            %only change if binning is changed while using same objective magnification!!

%colormap('parula');


% Load previously determined sequence of masks


load([rawdir,filesep,'MaskData_raw.mat']);


% perform tracking of mask centroid
% min(cellfun(@(x) size(x,1),cellCoors))    % changed it from min to max ? 
%% %% delete flash glitch frames that mess FRET data - if needed 
% % % % 
% flashFrame=1:29; 
% maskFinal(flashFrame)=[];
% imFRETOutline(flashFrame)=[];
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
centroid_coordinates=zeros(2,size(gCAMPOutline,2));
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
 distance = distance/(size(gCAMPOutline,2)*(2/3))
 avg_cell_area = cell_area/size(gCAMPOutline,2)
 
 %% filtered protusionvalues
 
   % filtered protrusion values
protvalsWindowF=ndnanfilter(protvalsWindow,fspecial('disk',2),'replicate');
protvalsrangeF=[round(prctile(protvalsWindowF(:),1),1),round(prctile(protvalsWindowF(:),99),1)];


%% display location of coordinates overlayed on cell mask images - optional for visualization
%  % creates gif file for reference
%     for mapper =1:size(thisTraj,1)
%      
%         window = windowCoors{1,mapper}; 
%       if ~(isConnect)
%       image = (imFRETOutline{1,thisTraj(mapper,5)+empty_count}); %slightly convoluted but lines up YFP raw image to the correct windowcoors
%       else
%            image = imread(imFRETOutline{1,thisTraj(mapper,4)+empty_count});
%       end 
%       h = figure('visible','off');
%       hold on; 
%       axis ij; 
%      imagesc(image); 
   
      %for num = 1:5:size(window,1)
%       for num = 100:15:115
%          image(window(num,1), window(num,2),1) = 255; 
%          image(window(num,1), window(num,2),2) = 0; 
%          image(window(num,1), window(num,2),3) = 0; 
%          text(window(num,2),window(num,1),num2str(num), 'Color','r', 'FontSize', 8);
%       end 
     
   %   frame = getframe(h);
   %  im=frame2im(frame);
    % [imind, cm] = rgb2ind(im,256);
     
%          if mapper ==1
%              imwrite(imind,cm,strcat(datadir,'\','Coordinate_Windows'),'gif','Loopcount', inf);
%          else
%                    imwrite(imind,cm,strcat(datadir,'\','Coordinate_Windows'),'gif','WriteMode','append');
%          end 
   %    imwrite(image,[datadir,filesep,'Outline_label.tif'],'WriteMode','append','Compression','none');

  %end 
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


ax1=subplot(1,2,1);imagesc(protvalsWindow,[-13,13]);title('Edge Velocity');
axis square; 
colormap(ax1,cmap);
xlabel("time (min)"); 
xline(76, '--', 'LineWidth',2); 
%15s intervals 
%  xticks([40 80 120 160 200]); 
% xticklabels({'0','10','20','30','40','50'}); 
%25 s intervals 
   xticks([0 24 48 71 95 120 144])
   xticklabels({'0','10','20','30','40','50' '60'});

ax2=subplot(1,2,2);imagesc(protvalsWindowF,[-13 13] );title('Edge Vel Filtered');
xlim([0 150]); 
colormap(ax2,cmap);
axis square; 

xlabel("time (min)"); 
xline(76, '--', 'LineWidth',2); 

 %15s intervals 
%  xticks([40 80 120 160 200]); 
% xticklabels({'0','10','20','30','40','50'}); 
%25 s intervals 
 xticks([24 48 71 95 120 144])
  xticklabels({'0','10','20','30','40','50'});


protvalsWindowHigh=protvalsWindow>protthresh;



 hold off; 
saveas(f1,strcat(datadir,'\','edge_velocity_map.png'))
saveas(f1,strcat(datadir,'\','edge_velocity_mapM.fig'))

%save('protrusionpeaks.mat','protvalsWindowF','protvalsWindowFHigh')

% save all new data into mat file 

save(strcat(datadir,'\','Protrusion Values.mat'),'protvalsWindow','protvalsWindowF','distance','avg_cell_area'); 


%close all; clc;

end 
clear; clc; 