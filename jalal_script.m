%% jalal segmentation and cropping script 
clc; clear; 
numStacks = 2; 
for i = 1:numStacks
rawdir = (['G:\imaging 2023\Jalal imaging_Results\Protrusion assay optimization trial 4\2024-04-17\analysis\', strcat('Stack_', num2str(i))]); 


minCellSize = 100; 


    
    cellFiles=getFilenames([rawdir],'.tif');
    
    for cell = 1:size(cellFiles, 1)
        datadir = ([rawdir,filesep, 'Mask Data ', num2str(cell)]); 
        if ~exist(datadir)
         mkdir(datadir)
        end 

%check here order of cellFiles!!! 
 image_stack=double(readTIFFstack([rawdir,filesep,cellFiles{cell}]));
    
 maskFinal = {}; 
    for frameNum = 1:121
        
        
        image_raw=image_stack(:,:,frameNum);
        
        mask = image_raw > 0; 
        mask = bwareafilt(mask,1); % this keeps only the largest segmented mask if multiple in feild of view 
        maskFinal{frameNum}=mask;
        
        objects=regionprops(mask,'BoundingBox','Area','Centroid','PixelIdxList','PixelList');
        objects=objects(arrayfun(@(x) x.Area>minCellSize,objects)); %

cellCoorstemp(:,1)=vect(arrayfun(@(x) x.Centroid(1),objects));
cellCoorstemp(:,2)=vect(arrayfun(@(x) x.Centroid(2),objects));
cellCoorstemp(:,3)=vect(arrayfun(@(x) x.Area,objects));

 cellCoors{frameNum}=cellCoorstemp;


save([datadir,filesep,'MaskData_raw.mat'],'maskFinal','cellCoors','-v7.3'); 
    end 
    end 
end 
    
    %% edge velcoity measurements 
    
   
 numStacks = 2; 
for i = 1:numStacks
rawdir = (['G:\imaging 2023\Jalal imaging_Results\Protrusion assay optimization trial 4\2024-04-17\analysis\', strcat('Stack_', num2str(i))]); 

cellFiles=getFilenames([rawdir],'.tif');
    
    for cell = 1:size(cellFiles, 1)
        datadir = ([rawdir,filesep, 'Mask Data ', num2str(cell)]); 
        if ~exist(datadir)
         mkdir(datadir)
        end 


% Parameters for cell edge parametrization
nFretWindows=180;                                   % Number of windows to use for FRET vs edge correlation measurements
edgeOversamplingParam=5;                            % How many times more points should the edge of the cell be tracked at ("subwindow")
nPointsParam=nFretWindows*edgeOversamplingParam;    % Number of points to track on the cell edge
pdSmoothing=10;%was 10                                      % Used with imclose to make the selection of points for tracking less dependent on noise or wrinkles in the cell edge



binning=1;            %only change if binning is changed while using same objective magnification!!

%colormap('parula');


% Load previously determined sequence of masks


load([datadir,filesep,'MaskData_raw.mat']);


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






%% select cell for analysis 
selectedCell=1; %input which trajectory
isConnect = false; %true if you had to connect broken trajectories
thisTraj=traj{selectedCell};
centroid_coordinates=zeros(2,size(cellCoors,2));
cell_area=0; 
if ~(isConnect) %modifying this Traj code if you needed to correct broken trajectores
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
                    
                    protvals(:,index-1)=vect(computeProtrusionValues(edgeCoors{index-1},edgeCoorsSmoothed{index-1},edgeCoors{index}));
                end
                windowCoors{index}=edgeCoors{index}((edgeOversamplingParam+1)/2:edgeOversamplingParam:end,:);

            % Compute edge regions and local fret ratio values
         %   labelMask{index}=getWindowLabelMap(thisMask,windowCoors{index},(edgeDepthDist/binning));
         
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
 distance = distance/(size(cellCoors,2))*2;
 avg_cell_area = cell_area/size(cellCoors,2); 
 
 %% filtered protusionvalues
 
   % filtered protrusion values
protvalsWindowF=ndnanfilter(protvalsWindow,fspecial('disk',2),'replicate');

%% %% Plot maps - with thresholds
load('C:\Users\marsh\OneDrive - McGill University\Documents\GitHub\Rodrigo_Codes\Colormaps\BCWOR-256.mat'); 

close all;

protthresh=5;
retthresh=-5;
f1=figure; 


hold on;



protvalrange=[round(prctile(protvalsWindow(:),1),1),round(prctile(protvalsWindow(:),99),1)];

ax1=subplot(1,2,1);imagesc(protvalsWindow,[-10, 10]);title('Edge Velocity');
axis square; 
colormap(ax1,BCWOR);


%25 s intervals 
   xticks([0 20 40 60 80 100 120]) % frames 
   xticklabels({'0','10','20','30','40','50' '60'}); % correspdoning minutes 



protvalsWindowHigh=protvalsWindow>protthresh;

ax2=subplot(1,2,2);imagesc(protvalsWindowF, [-10, 10]); title ('Filtered');
colormap(ax2,BCWOR);
axis square;
xlim([0 121]); 
  xticks([ 0 20 40 60 80 100 120])
  xticklabels({'0','10','20','30','40','50'});
 hold off; 

 
saveas(f1,strcat(datadir,'\','edge_velocity_map.png'))
saveas(f1,strcat(datadir,'\','edge_velocity_mapM.fig'))

% save all new data into mat file 

save(strcat(datadir,'\','Protrusion and FRET Values.mat'),'protvalsWindow','protvalsWindowF','distance', 'avg_cell_area') ; %

%% 

%% quantify number of protrusions retractions 

% inputs to this function are:
%file path of cell being analyszed
% which frame you want to start analysis (default 1)
%protrusion threshold in pixels/frame 
%retraction threshold in px/frame 
% minimum event size (default 25)
 getEdgeVelStats_edits_JR(datadir,1,5,-5,25);




    end 
end 

