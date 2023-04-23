clc; clear; 
close all; 

root  = 'F:\Seph\data\data_200116 - Trial 3 Rho, Myosin'; 
event = 1; 
rawdir = ([root, filesep,'cropped',filesep, 'cell_8']); 
datadir =([root, filesep,'+Y FRET buildup - retractions', filesep, num2str(event)]); 

if  ~exist(datadir)
    mkdir(datadir)
end 

depths =[3,6,10,15,20,25]; 

FRET_depths = cell(2,6);

load([rawdir, filesep,'edge vel mapping_3',filesep,'Protrusion and FRET values.mat'],'protvalsWindowF', 'fretvalsF'); 

for row = 1: size(protvalsWindowF,1)
        for col = 1:size(protvalsWindowF,2)-1
          temporary(row,col)=(protvalsWindowF(row,col)+protvalsWindowF(row,col+1))/2;
        end
end

protval_map = temporary; 

for i=1:size(depths, 2)
   
    load([rawdir, filesep,strcat('edge vel mapping_',num2str(depths(1,i))),filesep,'Protrusion and FRET values.mat'],'fretvalsF')
    
    depth =  num2str(depths(1,i)); 
    
   FRET_depths{1,i}=depth;  
   FRET_depths{2,i}= fretvalsF(:,2:end-1);
    
end 


load('F:\Seph\code\supporting_functions\trackingcode\CMAP_blue_grey_yellow.mat');
w=figure;
imagesc(protval_map,[-13,13]);title('Edge Velocity');
colormap(w,cmap);


 
%% specify rectangle dimensions 
%close all; 
f = figure; 
 
imagesc(protvalsWindowF,[-12 12]);title('Edge Velocity');
colormap(f,cmap);
xticks([0 24 48 71 95 120 144])
 xticklabels({'0','10','20','30','40','50', '60'});  
hold on; 
coor_start = 100;
coor_end = 115;
time_start=51;
time_end =91; 

 
rectangle('Position',[time_start,coor_start,time_end-time_start,coor_end-coor_start],'LineWidth',2); 

hold off; 

f2 = figure; 


region = protval_map(coor_start:coor_end, time_start:time_end); 



% calculate average velocity in specified area 
for t = 1:size(region,2)
vel_avg(1,t) = mean(region(:,t)); 
end 


hold on; 
%x = [time_start:1:time_end]; 
x = [-20:1:20]; 



%find timept of max acceleration (pos or neg depending on analysis
acceleration=zeros(1,size(vel_avg,2)-1);
for i=1:size(vel_avg,2)-1
    acceleration(1,i)=[vel_avg(1,i+1)-vel_avg(1,i)]; 
end 

% changre this depending on if a protrusion or retraction 
for j = 1:size(vel_avg,2)
   zero_pt = (acceleration(1,j) <0  && vel_avg(1,j) >= 0 && vel_avg(1,j+1) <= 0); 
   if (zero_pt ==1)
       
       break;   
   end 

end 

if (abs(vel_avg(1,j+1)) < abs(vel_avg(1,j)))

    align_pt =j+1; 
    
else 
    align_pt = j; 

end 
xline(0); 
  xticks([-19.2 -16.8 -14.4 -12.0 -9.6 -7.2 -4.8 -2.4 0 2.4 4.8 7.2 9.6 12.0 14.4 16.8 19.2])
xticklabels({'-8' '' '-6' '' '-4' '' '-2' '' '0' '' '2' '' '4' '' '6' '' '8'}); 

yyaxis right;
ylabel('Dora-RhoB');
ylim([0.5 1.5]); 

yline(1, '--k'); 

FRET_temp=cell(2,6); 
for z =1:6
FRET_temp{1,z}=depths(1,z); 
FRET_temp{2,z} = FRET_depths{2,z}(coor_start:coor_end,time_start+align_pt-1-20:time_start+align_pt-1+20); 
FRET_temp{2,z}(isnan(FRET_temp{2,z}))=0; 
for w=1:size(FRET_temp{2,z},2)
FRET_avg(1,w) =mean(FRET_temp{2,z}(:,w));     
end 

plot(x,FRET_avg); 
 
end 

yyaxis left;
ylim([-15.384 15.384]); 
yticks([-15.384 -12.82 -10.256 -7.692 -5.128 -2.564 0 2.564 5.128 7.692 10.256 12.82 15.384]); 
yticklabels({'-12' '-10', '-8', '-6', '-4', '-2', '0', '2', '4', '6', '8', '10' '12'}); 

 
vel_arr = vel_avg(1,align_pt-20:align_pt+20);
plot(x,vel_arr); 
ylabel('Edge Velocity (um/min)'); 
xlabel('Time (min)'); 

 

%legend('Edge Vel','Peak Retraction','0.975 um','1.95 um', '3.25 um', '4.88 um', '6.50 um','8.13 um','Location','southoutside'); 

%saveas(f, [datadir,filesep,'Vel Map.png']); 
%saveas(f2,[datadir,filesep,'retraction']); 

%save(strcat(datadir,'\','retraction_statistics.mat'),'time_start','time_end','coor_start','coor_end','FRET_temp','vel_arr'); 
