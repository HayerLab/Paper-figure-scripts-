%% analyze cell motility data from siRNA trials 
clc; clear; 

bare_root = 'E:\seph backup\LOK SLK ERM KD\cropped'; 
datadir = ([bare_root, filesep, 'siRNA motility compare']); 

if ~exist(datadir)
    mkdir(datadir)
end 

%% % siCNTRL 
 root = 'E:\seph backup\LOK SLK ERM KD\cropped\siNT';
 cells = [1,2,3,4,5,6,7,9,10,11,12,13,14,15,16,17,18,19,20];  
 
 distance_arr = []; 
 cell_size = []; 
 num_protrusions = []; 
 num_retractions = []; 
 
 prot_speed = []; 
 retr_speed = []; 
 prot_area = []; 
 retr_area = []; 
 
 for i = 1:length(cells)
     
     load([root, filesep, num2str(cells(1,i)), filesep, 'edge_vels\edge vel mapping_6', filesep, 'Protrusion and FRET Values.mat']); 
          load([root, filesep, num2str(cells(1,i)), filesep, 'edge_vels\edge vel mapping_6', filesep, 'protrusionlist.mat']);
          load([root, filesep, num2str(cells(1,i)), filesep, 'edge_vels\edge vel mapping_6', filesep, 'retractionlist.mat']);
     
   distance_arr = [distance_arr, distance]; 
   cell_size = [cell_size, avg_cell_area]; 
   
   num_protrusions = [num_protrusions, size(protrusions, 1)/((size(protvalsWindow,2))*(2/3))]; 
   num_retractions = [num_retractions,  size(retractions, 1)/((size(protvalsWindow,2))*(2/3))]; 
   
    
   prot_speed = [prot_speed, mean([protrusions.mspeed], 'all')]; 
   retr_speed = [retr_speed, mean([retractions.mspeed], 'all')]; 
   prot_area = [prot_area, mean([protrusions.Area], 'all')]; 
   retr_area = [retr_area, mean([retractions.Area], 'all')]; 
 end    

 %%  
 root = 'F:\Seph\data\LOK SLK Msn KD 1\cropped\siSLKLOK';
 %cells = [1,2,3,4,5,7,8,9,11,12,14,15,16,17,18,19,20]; 
 cells = [1,2,3,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19]; 
 distance_arr_slk = []; 
 cell_size_slk = []; 
 num_protrusions_slk = []; 
 num_retractions_slk = []; 
 
 prot_speed_slk = []; 
 retr_speed_slk = []; 
 
 for i = 1:length(cells)
     
     load([root, filesep, num2str(cells(1,i)), filesep, 'edge_vels\edge vel mapping_6', filesep, 'Protrusion and FRET Values.mat']); 
          load([root, filesep, num2str(cells(1,i)), filesep, 'edge_vels\edge vel mapping_6', filesep, 'protrusionlist.mat']);
          load([root, filesep, num2str(cells(1,i)), filesep, 'edge_vels\edge vel mapping_6', filesep, 'retractionlist.mat']);
     
   distance_arr_slk = [distance_arr_slk, distance]; 
   cell_size_slk = [cell_size_slk, avg_cell_area]; 
   
   num_protrusions_slk = [num_protrusions_slk, size(protrusions, 1)/((size(protvalsWindow,2))*(2/3))]; 
   num_retractions_slk = [num_retractions_slk,  size(retractions, 1)/((size(protvalsWindow,2))*(2/3))]; 
   
    
   prot_speed_slk = [prot_speed_slk, mean([protrusions.mspeed], 'all')]; 
   retr_speed_slk = [retr_speed_slk, mean([retractions.mspeed], 'all')]; 
   
  
   
 end   
 
 %%
  root = 'E:\seph backup\LOK SLK ERM KD\cropped\siERM';
  %cells = [1,2,3,4,5,6,8,9,10,11,14,15,16,18];  
  % T2 cells = [1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,19,20,21,23,24,26,28,29]; 
   cells = [1,3,4,5,6,7,9,10,11,12,13]; 

 distance_arr_erm= []; 
 cell_size_erm = []; 
 num_protrusions_erm = []; 
 num_retractions_erm = []; 
 
 prot_speed_erm= []; 
 retr_speed_erm = []; 
 
 prot_area_erm = []; 
 retr_area_erm  = []; 
 
 for i = 1:length(cells)
     
     load([root, filesep, num2str(cells(1,i)), filesep, 'edge_vels\edge vel mapping_6', filesep, 'Protrusion and FRET Values.mat']); 
          load([root, filesep, num2str(cells(1,i)), filesep, 'edge_vels\edge vel mapping_6', filesep, 'protrusionlist.mat']);
          load([root, filesep, num2str(cells(1,i)), filesep, 'edge_vels\edge vel mapping_6', filesep, 'retractionlist.mat']);
     
   distance_arr_erm = [distance_arr_erm, distance]; 
   cell_size_erm = [cell_size_erm, avg_cell_area]; 
   
   num_protrusions_erm= [num_protrusions_erm, size(protrusions, 1)/((size(protvalsWindow,2))*(2/3))]; 
   num_retractions_erm = [num_retractions_erm,  size(retractions, 1)/((size(protvalsWindow,2))*(2/3))]; 
   
    
   prot_speed_erm = [prot_speed_erm, mean([protrusions.mspeed], 'all')]; 
   retr_speed_erm = [retr_speed_erm, mean([retractions.mspeed], 'all')]; 
   
   prot_area_erm = [prot_area_erm, mean([protrusions.Area], 'all')]; 
   retr_area_erm = [retr_area_erm, mean([retractions.Area], 'all')]; 
   
 end   
 
 %% 
 
 number_protrusions = struct; 
  number_protrusions.CTRL = num_protrusions; 
 % number_protrusions.SLK = num_protrusions_slk; 
  number_protrusions.ERM = num_protrusions_erm; 
 f= figure; 

grid on; 
violinplot(number_protrusions); 

% yticks(-100:10:40); 
% ylim([-100, 40]); 

% 
  y = [ num_protrusions, num_protrusions_erm]; %num_protrusions_slk,
  group = repelem(1:2, 1, [numel(num_protrusions),numel(num_protrusions_erm)]); %numel(num_protrusions_slk),

  [p, anova, stats] = anova1(y, group); 
 
 [compare]= multcompare(stats); 
  
save([datadir,filesep,  'num_protrusions.mat'], 'number_protrusions'); %'compare'); 

saveas(f, ([datadir, filesep, 'num_protrusions.fig'])); 
%%
 number_retractions = struct; 
 number_retractions.CTRL = num_retractions; 
% number_retractions.SLK = num_retractions_slk; 
 number_retractions.ERM = num_retractions_erm; 
 f= figure; 

grid on; 
violinplot(number_retractions); 

% yticks(-100:10:40); 
 ylim([0 0.6]); 

% 
  y = [num_retractions,  num_retractions_erm]; %num_retractions_slk,
  group = repelem(1:2, 1, [numel(num_retractions),numel(num_retractions_erm)]); %),numel(num_retractions_slk),

  [p, anova, stats] = anova1(y, group); 
 
 [compare]= multcompare(stats); 
  
save([datadir,filesep,  'num_retractions.mat'], 'number_retractions'); %  'compare'); 

saveas(f, ([datadir, filesep, 'num_retractions.fig'])); 
%%  
 distance = struct; 
  distance.CTRL = distance_arr; 
 % distance.SLK = distance_arr_slk; 
 distance.ERM = distance_arr_erm; 
 f= figure; 

grid on; 
violinplot(distance); 

% yticks(-100:10:40); 
 ylim([0.5 , 5]); 

% 
  y = [distance_arr, distance_arr_erm]; % distance_arr_slk,
  group = repelem(1:2, 1, [numel(distance_arr),numel(distance_arr_erm)]); %numel(distance_arr_slk)

  [p, anova, stats] = anova1(y, group); 
 
 [compare]= multcompare(stats); 
  
save([datadir,filesep,  'distance.mat'], 'distance'); 

saveas(f, ([datadir, filesep, 'distance.fig'])); 
 
%% 
 retraction_speed = struct; 
  retraction_speed.CTRL = retr_speed; 
 % retraction_speed.SLK = retr_speed_slk; 
 retraction_speed.ERM = retr_speed_erm; 
 f= figure; 

grid on; 
violinplot(retraction_speed); 

% yticks(-100:10:40); 

%ylim([-6 -3 ]); 

% 
  y = [retr_speed,retr_speed_erm]; % , retr_speed_slk
  group = repelem(1:2, 1, [numel(retr_speed),numel(retr_speed_erm)]); %,numel(retr_speed_slk),

  [p, anova, stats] = anova1(y, group); 
 
 [compare]= multcompare(stats); 
  
save([datadir,filesep,  'retraction_speed.mat'], 'retraction_speed'); % 'compare'); 

saveas(f, ([datadir, filesep, 'retraction_speed.fig'])); 

%% 

retraction_area = struct; 
  retraction_area.CTRL = retr_area; 
  %Cell_area.SLK = retr_area_slk; 
 retraction_area.ERM = retr_area_erm; 
 f= figure; 

grid on; 
violinplot(retraction_area); 

% yticks(-100:10:40); 

%ylim([-6 -3 ]); 

% 
  y = [retr_area, retr_area_erm]; % retr_area_slk,
  group = repelem(1:2, 1, [numel(retr_area),numel(retr_area_erm)]); %numel(retr_area_slk),

  [p, anova, stats] = anova1(y, group); 
 
 [compare]= multcompare(stats); 
  
save([datadir,filesep,  'retraction_area.mat'], 'retraction_area'); % 'compare'); 

saveas(f, ([datadir, filesep, 'retraction_area.fig'])); 