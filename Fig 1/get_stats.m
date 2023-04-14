%% run stats on cell motility statistics 

root = 'C:\Users\gmarsh8\OneDrive - McGill University\research paper\results good_Feb2023\Fig 1\cell mvmt statistics'; 

load ([root, filesep, 'Rho T3 stats.mat']); 
stats_rho = stats; 

load ([root, filesep, 'RacT4 stats.mat']); 
stats_rac = stats; 


load ([root, filesep, 'Cdc42 T2 stats.mat']); 
stats_Cdc42 = stats; 

number_protrusions =[]; 
number_retractions = []; 
prot_size =[]; 
retr_size = []; 
prot_speed = []; 
retr_speed = []; 
prot_perimeter =[]; 
retr_perimeter = []; 
prot_timelength = []; 
retr_timelength = [];

number_protrusions =[number_protrusions, stats_rho.num_prots]; 
number_retractions = [number_retractions, stats_rho.num_retr]; 
prot_size =[prot_size,stats_rho.prot_avg_size]; 
retr_size = [retr_size, stats_rho.ret_avg_size]; 
prot_speed = [prot_speed, stats_rho.prot_speed]; 
retr_speed = [retr_speed, stats_rho.ret_speed]; 
prot_perimeter =[prot_perimeter, stats_rho.prot_perimeter]; 
retr_perimeter = [retr_perimeter, stats_rho.retr_perimeter]; 
prot_timelength = [prot_timelength, stats_rho.prot_timelength]; 
retr_timelength = [retr_timelength, stats_rho.retr_timelength];

number_protrusions =[number_protrusions, stats_rac.num_prots]; 
number_retractions = [number_retractions, stats_rac.num_retr]; 
prot_size =[prot_size,stats_rac.prot_avg_size]; 
retr_size = [retr_size, stats_rac.ret_avg_size]; 
prot_speed = [prot_speed, stats_rac.prot_speed]; 
retr_speed = [retr_speed, stats_rac.ret_speed]; 
prot_perimeter =[prot_perimeter, stats_rac.prot_perimeter]; 
retr_perimeter = [retr_perimeter, stats_rac.retr_perimeter]; 
prot_timelength = [prot_timelength, stats_rac.prot_timelength]; 
retr_timelength = [retr_timelength, stats_rac.retr_timelength];



number_protrusions =[number_protrusions, stats_Cdc42.num_prots]; 
number_retractions = [number_retractions, stats_Cdc42.num_retr]; 
prot_size =[prot_size,stats_Cdc42.prot_avg_size]; 
retr_size = [retr_size, stats_Cdc42.ret_avg_size]; 
prot_speed = [prot_speed, stats_Cdc42.prot_speed]; 
retr_speed = [retr_speed, stats_Cdc42.ret_speed]; 
prot_perimeter =[prot_perimeter, stats_Cdc42.prot_perimeter]; 
retr_perimeter = [retr_perimeter, stats_Cdc42.retr_perimeter]; 
prot_timelength = [prot_timelength, stats_Cdc42.prot_timelength]; 
retr_timelength = [retr_timelength, stats_Cdc42.retr_timelength];

[number_p_val]= ranksum(number_protrusions, number_retractions) 
[size_p_val]= ranksum(prot_size, retr_size)
[speed_p_val]= ranksum(prot_speed, abs(retr_speed))
[time_p_val]= ranksum(prot_timelength, retr_timelength)
[perimeter_p_val]= ranksum(prot_perimeter, retr_perimeter) 

 save('C:\Users\gmarsh8\OneDrive - McGill University\research paper\results good_Feb2023\Fig 1\cell mvmt statistics\compiled_stats.mat', 'number_protrusions','number_retractions', 'prot_size', 'retr_size','prot_speed', 'retr_speed', 'prot_perimeter', 'retr_perimeter', 'prot_timelength', 'retr_timelength'); 
 
 save('C:\Users\gmarsh8\OneDrive - McGill University\research paper\results good_Feb2023\Fig 1\cell mvmt statistics\compiled_p vals.mat','number_p_val', 'size_p_val', 'speed_p_val','time_p_val', 'perimeter_p_val'); 
 
 


