%line plots to display cytopasm bias
root = 'C:\Users\gmarsh8\OneDrive - McGill University\research paper\results good_Feb2023\Supp Fig 1\YFP-FRET lineplots and perpixel corr examples\RhoB'; 
YFP_table=readtable([root, filesep,'Plot Values-YFP.csv']); 
FRET_table=readtable([root, filesep, 'Plot Values-FRET.csv']); 

YFP_arr = table2array(YFP_table); 
FRET_arr = table2array(FRET_table); 
YFP_max = max(YFP_arr(:,2)); 
FRET_max = max(FRET_arr(:,2)); 

%normalization
for i = 1:size(FRET_arr,1)
    FRET_arr(i,1) = FRET_arr(i,1)./FRET_arr(size(FRET_arr,1),1); 
     FRET_arr(i,2) = FRET_arr(i,2)./FRET_max; 
    YFP_arr(i,1) = YFP_arr(i,1)./YFP_arr(size(YFP_arr,1),1); 
    YFP_arr(i,2) = YFP_arr(i,2)./YFP_max; 
    
    
   
end 


f = figure; 
hold on; 
xlabel('distance edge to nucleus (norm.)'); 

yyaxis left; 

plot(YFP_arr(:,1), YFP_arr(:,2)); 
ylabel('YFP Norm.'); 

yyaxis right; 
ylim([0.8 1]); 
 
ylabel('FRET Norm.'); 

plot(FRET_arr(:,1), FRET_arr(:,2)); 

saveas(f, [root, filesep, 'lineplot.fig']); 
saveas(f, [root, filesep, 'lineplot.svg']); 