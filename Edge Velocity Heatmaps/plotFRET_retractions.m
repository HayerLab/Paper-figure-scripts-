%script that loads retraction/protrusionlist mat files created from
%getEdgeVelStats_edits

% graphing averages for all cells -- large for loop
clear; clc; 
root='F:\Seph\data\data_200709 - Trial 1 Rac\cropped';


% see google sheets for these 
%cells=[6,7,9,10,12,13,14,15,16,17,18]; %trial 1
%cells = [2,3,5,6,7,8,9,11,12,13]; % trial 2
%cells =[2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19]; %trial 3
%cells=[2,3,4,5,7,8,9,11]; Rac trial 1 
%cells = [12,18,19,20,24,26,27,28,29,30,31,32,33]; %CDC42 t1
 %cells = [1,2,4,5,8,9,10,11]; 
%cells = [10,15,19]; 

%cells =[1,2,4,6,7,8,10,12,13,14,15,16,17,18,19]; CDC 42 T2
%cells = [1,4,5,6,7,8,9,11,12,13,14,15,17,18,21,22]; 
cells = [2,3,4,5,7,8,9,11]; 
cell_arr=cell(2,size(cells,2),1);

depths = [3,6,10,15,20,25]; 

for i = 1:6 
    depth = num2str(depths(1,i)); 

for loop=1:size(cells,2)
    fileKey=strcat('cell_',num2str(cells(1,loop)));
    
    load([root,filesep,fileKey,filesep,strcat('edge vel mapping_',depth),filesep,'protrusionlist.mat'],'protrusions');
    load([root,filesep,fileKey,filesep, strcat('edge vel mapping_',depth),filesep,'retractionlist.mat'],'retractions');
    load([root,filesep,fileKey,filesep,strcat('edge vel mapping_',depth),filesep,'Protrusion and FRET values.mat']);
    
     FRET_normed_retractions=zeros(size(retractions,1),1);
 
 for i=1:size(retractions,1)
     
   
     sum =0;
    
     counter=0;
     
         for j=retractions(i).framestart:retractions(i).frameend
            col =protvalsWindowF(:,j-1);
            logicalIndexes= (col >-5) & (col <5);  


            for k=1:180
                 if (logicalIndexes(k,1))
                     sum=sum+fretvalsF(k,j);
                     
                     counter=counter+1;
                 end     
             end   

mean_neutral=sum/counter; 
         end 

FRET_normed_retractions(i)=retractions(i).fretavg/mean_neutral;
 end 
     
 cell_arr{1,loop}=FRET_normed_retractions;    
 
 %protrusions
 FRET_normed_protrusions=zeros(size(protrusions,1),1);
 
 for i=1:size(protrusions,1)
     
   
     sum =0;
     counter=0;
     
         for j=protrusions(i).framestart:protrusions(i).frameend
            col =protvalsWindowF(:,j-1); 
                                          
            logicalIndexes= (col >-5) & (col <5);  


            for k=1:180
                 if (logicalIndexes(k,1))
                     sum=sum+fretvalsF(k,j);
                     counter=counter+1;
                 end     
             end   

mean_neutral=sum/counter; 
         end 

FRET_normed_protrusions(i)=protrusions(i).fretavg/mean_neutral;

 end
 
 cell_arr{2,loop}=FRET_normed_protrusions;
 
  f=figure; hold on; 
  
  ax(1)=subplot(1,2,1); 
notBoxPlot(FRET_normed_retractions);
xticks(1);
xticklabels('Retractions');

ax(2)= subplot(1,2,2);

notBoxPlot(FRET_normed_protrusions);
xticks(1);
xticklabels('Protrusions');
linkaxes(ax,'y');

saveas(f, strcat(root,'\',fileKey,filesep, 'edge vel mapping_',depth,'\FRET values graph.png')); 
    
end

%% graphs all protrusions and retractions at the same time 
close all;
    
    protrusionsfinal =[];
    retractionsfinal=[];
    
for i=1:size(cell_arr,2)
    protrusionsfinal=[protrusionsfinal cell_arr{2,i}'];
    
    retractionsfinal=[retractionsfinal cell_arr{1,i}'];
     
end 

f2=figure; 
hold on;

 ax(2)=subplot(1,2,2); 

 [R,stats]= notBoxPlot(retractionsfinal,'jitter',0.01);
   ylim([0.7 1.3]);
    xticks(1);
    xticklabels('Retractions');

    ax(1)= subplot(1,2,1);

   [P,stats2] = notBoxPlot(protrusionsfinal,'jitter',0.2);
   ylim([0.7 1.4]);
    xticks(1);
    xticklabels('Protrusions');
    linkaxes(ax,'y');
    
   
hold off; 
    
saveas(f2,strcat(root,'\graphs depth_', depth, '\overall averages Rho\','final graph.png'));

%% plots histogram of protrusions/retraction FRET vakues
close all; 
[N, edges] = histcounts(retractionsfinal,20); 
[M,edges2] = histcounts(protrusionsfinal,20);
f3=figure; 

max_M = max(M);
max_N=max(N);

%max = max(max_N, max_M); 
ax_front = axes;
   b_front=barh(edges(2:end),N);
   axis([-(max_M+3),(max_N)+3,0.5,1.5])
   set(b_front,'facecolor',[0.2,0.4,1])
   axis off; 
   
   ax_back=axes; 
   b_back=barh(edges2(2:end),-M);
   axis([-(max_M+3),max_N+3,0.5,1.5])
   set(b_back,'facecolor',[1,0.4,0.2])
  % xticks([-18:3:18])
  set(gca, 'xtick', [-84:6:84])
set(gca, 'xticklabel', [[84:-6:0],[6:6:84]])
grid on

   grid on; 
   axes(ax_front);
   
   saveas(f3,strcat(root,'\graphs depth_', depth, '\overall averages Rho\','final histo.svg'));
   
   
 %% creates individual histograms for figures 
 f5=figure; 
 
 hold on; 
 ylim([0.5 1.4]); 
 

 
 [line1, xi]=ksdensity(protrusionsfinal,'kernel','normpdf','npoints',size(protrusionsfinal,2)); 
 [line2, ai]=ksdensity(retractionsfinal,'kernel','normpdf','npoints',size(retractionsfinal,2)); 
 
 plot(line1,xi); 
 plot(line2,ai); 
 
%  direction = [0 1]; 
% rotate(f5, direction, 180); 
%  

% saveas(f5,'F:\Seph\research paper\Fig 1\F_v.svg'); 

 
 
  histogram(protrusionsfinal,'Normalization','probability', 'DisplayStyle', 'stairs');
  histogram(protrusionsfinal,'Normalization','probability', 'DisplayStyle', 'stairs');
   
   %% plots mean and standard deviation without the dots 
   close all; 
f4=figure; 
   
   hold on 
    
    xlim([0 6])
   yline(1,'--'); 
  % ytitle(Normalized FRET Intensity); 
  xticks([2 4])
  xticklabels({'Protrusions', 'Retractions'}); 
  
  C = [protrusionsfinal retractionsfinal]; 
grp = [zeros(1,size(protrusionsfinal,2)),ones(1,size(retractionsfinal,2))];
  boxplot(C,grp,'labels',{'Protrusions', 'Retractions'}); 
  title(strcat('Rac - Depth',' ', depth)); 
  
  ylim([0.7 1.3]);
    

   
  % errorbar(2,mean(protrusionsfinal),std(protrusionsfinal),'-sk','CapSize',15,'LineWidth',3); 
 %  errorbar(4,mean(retractionsfinal), std(retractionsfinal),'-sk', 'CapSize',15,'LineWidth',3); 
   
  [p_val]= ranksum(protrusionsfinal, retractionsfinal) % 'Vartype', 'unequal'); %look at the ttest2 matlab documentation, chose option with enqual variance her
   
  group = {[1,2]};
  sigstar(group, p_val)
   
   hold off; 
   
  saveas(f4,strcat(root,'\graphs depth_', depth, '\overall averages Rho\','error bar.fig'));
   saveas(f4,strcat(root,'\graphs depth_', depth, '\overall averages Rho\','error bar sig.png'));

   save(strcat(root,'\graphs depth_', depth, '\overall averages Rho\','statistics,'),'stats','stats2','p_val','protrusionsfinal','retractionsfinal'); 
end 