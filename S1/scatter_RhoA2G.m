% used to generate scatter plot of 2 values of interest for each pixel in
% a cell. Data from RatioData.mat files work, for example. 

% chooses 2 frames at random for a given cell, creates scatter plot, and
% Pearson r coefficient of correlation

% Seph, Sept 14 2022
clear; clc; 
R_arr =[]; 
 root = 'D:\221209 - 40x 2x2 bin_RhoB_cyto\cropped';
 datadir = 'D:\221209 - 40x 2x2 bin_RhoB_cyto\cyto bias correlation-YFP-new';
 
if  ~exist(datadir)
    mkdir(datadir); 
end 

counter = 0;
R_total = 0; 

for cell = 1:9
   
    % specific cells that don't work from debris in frame 
    if cell ==3 || cell ==7
        continue; 
    end 
    
  
load([root, filesep, strcat(num2str(cell)), filesep, 'output-YFP cyto', filesep, 'RatioData.mat']); 

%generate random frames 
frame_1 =randi(size(imRatio, 2)); 
frame_2 = randi(size(imRatio, 2)); 
frame_3 = randi(size(imRatio, 2)); 

if frame_2 == frame_1
    frame_2 =  randi(size(imRatio, 2)); 
end 

if frame_2 == frame_3
    frame_2 =  randi(size(imRatio, 2)); 
end 
if frame_1 == frame_3
    frame_1 =  randi(size(imRatio, 2)); 
end 

frames = [frame_1, frame_2, frame_3]; 

for w = 1:3
    close all; 
counter = counter+1 

y = imRatio{1,frames(w)}; 
%x = im_mRuby{1, frames(w)}; 
x = im_YFP{1,frames(w)}; 
z = maskFinal{1,frames(w)}; 

%[imRatioPrelim, mRubyPrelim] = getCellMask_nucleus(x,y,z,datadir, counter) ; 


%y= imRatioPrelim; 
%x = mRubyPrelim; 
x = x(:); 
x =x(~isnan(x));
y = y(:); 
y = y(~isnan(y)); 

% %arbitrary cutoff to remove nucleus/high cytoplasm areas 
% x_cut = x(x<=5); 
% y_cut = y(x<=5); 

[r, p] = corrcoef(x, y); 
R_arr=[R_arr; r(1,2)]; 

f = figure; 
hold on; 
ylim([0.5, 1.5]); 
%xlim([0,2.5]); 
ylabel('RhoB');
xlabel('YFP signal'); 
c = ksdensity([x, y], [x, y]); 
scatter(x,y,2, c, 'fill');
%cb = colorbar(); 

str=sprintf('r= %1.2f',r(1,2));
T = text(max(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'right');
h = lsline;
linear_fit_data = polyfit(get(h,'xdata'),get(h,'ydata'),1); 

cur_frame = frames(w); 
deposit = ([datadir, filesep, num2str(counter)]); 
if  ~exist(deposit)
    mkdir(deposit); 
end 


save([deposit,filesep,  'site data.mat'],'cell','r','cur_frame', 'linear_fit_data');
saveas(f,strcat(deposit,'\','YFP vs. DORA RhoB.svg'))
saveas(f,strcat(deposit,'\','YFP vs. DORA RhoB.fig'))

end 
 
end 
 pd = fitdist(R_arr,'Normal')
   ci = paramci(pd)
save([datadir, filesep, 'R values.mat'], 'R_arr', 'pd','ci'); 
