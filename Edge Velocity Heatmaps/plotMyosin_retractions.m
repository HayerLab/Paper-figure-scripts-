%script that loads retraction/protrusionlist mat files created from
%getEdgeVelStats_edits or something

% graphing averages for all cells -- large for loop
clear; clc; 
root='F:\Seph\data\data_201128 - Trial 1, 2 Rho Actin\cropped\Y27362';

%cells=[6,7,9,10,12,13,14,15,16,17,18]; %trial 1
%cells = [2,3,5,6,7,8,9,11,12,13]; % trial 2
%cells =[2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19]; %trial 3

%cells = [1,4,5,6,7,8,9,11,12,13,14,15,17,18,21,22];  % R, Actin T1/2
 cells =[5,6,7,8,9,10,11,13,15,16,17,18,19,20,21,22,23,24];



%cells = [1,2,3,4,5,7,8]; 
cell_arr=cell(2,size(cells,2),1);

for loop=1:size(cells,2)
    fileKey=strcat(num2str(cells(1,loop)));
    
    load([root,filesep,fileKey,filesep,'edge vel mapping_25',filesep,'protrusionlist.mat'],'protrusions');
    load([root,filesep,fileKey,filesep,'edge vel mapping_25',filesep,'retractionlist.mat'],'retractions');
    load([root,filesep,fileKey,filesep,'edge vel mapping_25',filesep,'Protrusion and FRET values.mat']);
    
     myosin_normed_retractions=zeros(size(retractions,1),1);
 
 for i=1:size(retractions,1)
     
   
     sum =0;
    
     counter=0;
     
         for j=retractions(i).framestart:retractions(i).frameend
       
            col =protvalsWindowF(:,j-1);
            
            logicalIndexes= (col >-5) & (col <5);  


            for k=1:180
                 if (logicalIndexes(k,1))
                     sum=sum+myosinF(k,j);
                     
                     counter=counter+1;
                 end     
             end   

mean_neutral=sum/counter; 
         end 

myosin_normed_retractions(i)=retractions(i).myosinavg/mean_neutral;
 end 
     
 cell_arr{1,loop}=myosin_normed_retractions;    
 
 %protrusions
 myosin_normed_protrusions=zeros(size(protrusions,1),1);
 
 for i=1:size(protrusions,1)
     
   
     sum =0;
     counter=0;
     
         for j=protrusions(i).framestart:protrusions(i).frameend
            col =protvalsWindowF(:,j-1); 
                                          
            logicalIndexes= (col >-5) & (col <5);  


            for k=1:180
                 if (logicalIndexes(k,1))
                     sum=sum+myosinF(k,j);
                     counter=counter+1;
                 end     
             end   

mean_neutral=sum/counter; 
         end 

myosin_normed_protrusions(i)=protrusions(i).myosinavg/mean_neutral;

 end
 
 cell_arr{2,loop}=myosin_normed_protrusions;
 
  f=figure; hold on; 
  
  ax(1)=subplot(1,2,1); 
notBoxPlot(myosin_normed_retractions);
xticks(1);
xticklabels('Actin-Retractions');

ax(2)= subplot(1,2,2);


notBoxPlot(myosin_normed_protrusions);
xticks(1);
xticklabels('Actin-Protrusions');
linkaxes(ax,'y');

saveas(f, strcat(root,'\',fileKey,'\edge vel mapping_25\actin values norm.png')); 
    
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
   [R,stats]= notBoxPlot(retractionsfinal,'jitter',0.2);
      ylim([-1 4]);
    xticks(1);
    xticklabels('M-Retractions');
   

    ax(2)= subplot(1,2,1);

   [P,stats2] = notBoxPlot(protrusionsfinal,'jitter',0.2);
   ylim([-1 4]);
    xticks(1);
    xticklabels('M-Protrusions');
    linkaxes(ax,'y');
    
   
hold off; 
    
saveas(f2,strcat(root,'\graphs depth 25\overall averages Actin\','final graph norm.svg'));

%% plots histogram of protrusions/retraction myosin vakues
close all; 
[N, edges] = histcounts(retractionsfinal,20); 
[M,edges2] = histcounts(protrusionsfinal,20);

max_num= max(N,M); 
max_M = max(M);
max_N=max(N);


f3=figure; 
ax_front = axes;
   b_front=barh(edges(2:end),N);
   axis([-(max_M+3),max_N+3,0,4])
   set(b_front,'facecolor',[0.2,0.4,1])
   axis off; 
   
   ax_back=axes; 
   b_back=barh(edges2(2:end),-M);
   axis([-(max_M+3),max_N+3,0,4])
   set(b_back,'facecolor',[1,0.4,0.2])
  % xticks([-18:3:18])
  set(gca, 'xtick', [-84:6:84])
%set(gca, 'xticklabel', [[90:-3:0],[3:3:90]])
set(gca, 'xticklabel', [[84:-12:0],[6:6:84]])

   grid on; 
   axes(ax_front);
   
   saveas(f3,strcat(root,'\graphs depth 25\overall averages Actin\','final histo.svg'));

 %% plots mean and standard deviation without the dots 
   close all; 
f4=figure; 
   
   hold on 
    ylim([0 2]);
    xlim([0 6])
   yline(1,'--'); 
  % ytitle(Normalized FRET Intensity); 
%   xticks([1 2])
%   xticklabels({'Protrusions', 'Retractions'}); 
    
C = [protrusionsfinal retractionsfinal]; 
grp = [zeros(1,size(protrusionsfinal,2)),ones(1,size(retractionsfinal,2))];
  boxplot(C,grp,'labels',{'Protrusions', 'Retractions'}); 
 % title('M - Depth 3'); 
   
%    errorbar(2,mean(protrusionsfinal),std(protrusionsfinal),'-sk','CapSize',15,'LineWidth',3); 
%    errorbar(4,mean(retractionsfinal), std(retractionsfinal),'-sk', 'CapSize',15,'LineWidth',3); 
   
  [p_val,boolean]= ranksum(protrusionsfinal, retractionsfinal) %look at the ttest2 matlab documentation, chose option with enqual variance her
   
  group = {[1,2]};
  sigstar(group, p_val)
   
   hold off; 
   
  saveas(f4,strcat(root,'\graphs depth 25\overall averages Actin\','error bar.fig'));
  saveas(f4,strcat(root,'\graphs depth 25\overall averages Actin\','error bar sig.png'));
  
  save(strcat(root,'\graphs depth 25\overall averages Actin\','statistics'),'stats','stats2','p_val','protrusionsfinal','retractionsfinal'); 
