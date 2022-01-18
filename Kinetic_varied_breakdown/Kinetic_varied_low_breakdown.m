userpath='/project/sfinley_292/vtseruny/FinleyRotation/LCK_Ant_proj';
MPCTprofile=parallel.cluster.Local;
MPCTprofile.NumWorkers=24;
MPCTprofile.JobStorageLocation='/project/sfinley_292/vtseruny/FinleyRotation/LCK_Ant_proj';
parpool(MPCTprofile, 23);

parameters = [ 4.5, 1000, 250, 500, 500, 870, 540, 600, 300, 120, 25, 800, ...
               60, 600, 150, 150, 250, 20, 350, 1500, 1250, 0, 400, 1000, .1, ...
               0.023104906, 1000, 50, 40.26185694, 4.852880242, 6.373803376, ...
               .1, 416.8324872, 30.89383828, 137.0133105, 0.233721543, 0.1, ...
               2.026709313, 0.391443224, 1.310391485, 0.101597, 96.548, 331.68, ...
               270, 153.82, 484.9, 405.43, 0.1153, 360, 0.1, 500, 101, 28, 96, ...
               12.8, .1, 20, 0.05, 200, 200, 20, 0.1, 10, 0.1, 75, 0.1, 142, .1, ...
               110, .1, 15, 0.1, 55, 0.1, 250, 0.1, 65, .1, 16333, .1, 2400, 0.01, ...
               9666, 0.03, .1, 1266, 2.28, .1, 9666, .18, .1, 66, 6, 0.1, 186, 10, ...
               .1, 2000, .6, 12, 33, 200, 33, 60, 86, 300, 33, 10000, 66, 10000, ...
               66, 53, 60, 33, 1200, 33, 15200, 7, 2000, 12, 66, 66, 1, 0.0015, ...
               100, 100, 100, 10, 1.248432393, 0.094342447, 0.963404796, 67.84415559, ...
               1.459578998, 9.918819462, 44.15493881 ];

sim_size = 1e4;
rng(12121995);
ERK_times  = zeros(size(sim_size,1));
ERK_concs  = zeros(size(sim_size,1));
timepoints = 0:0.05:30;
timepoints = timepoints';
parameters(5) = 0;
parameters(42:49) = [96.548, 331.68, 270, 153.82, 484.9, 405.43, 0.1153, 360];
indices = [24 26 43 44 45 46 48 49 29 30 32 33 34 35 37 38 56 70 60 61 72 ...
           94 96 97 99 100 101 102 103 104 106 107 108 109 110 111 112 113 ...
           115 117 121 27 131 132 124 126 127 128];
mu = parameters(indices);
param_samples = zeros(sim_size, length(indices));

for ind = 1:length(indices)
    
    pdf = makedist('Normal','mu',mu(ind),'sigma',mu(ind)/3);
    pdf = truncate(pdf, 0, 2*mu(ind));
    param_values = random(pdf, sim_size, 1);
    param_samples(:,ind) = param_values;
    
end

%% Only the CD3z
parameters(5)  = 0;   %Add CD28
parameters(64) = 0.1; %Take away GADS
parameters(66) = 0.1; %Add Grb2
parfor ind=1:sim_size
   
   params     = parameters;
   params(indices) = param_samples(ind,:);
   try
       [~, tpts, ~, observables_out] = model_func(timepoints, params);
   catch
       sprintf('Integration error occured at ind = %f', ind)
       ERK_times(ind) = NaN;
       ERK_concs(ind) = NaN;
       continue
   end
   
   ERK_pp = observables_out(:,133);
   halftime_index_curr = sum(ERK_pp<125);
   ERK_times(ind) = timepoints(halftime_index_curr);
   ERK_concs(ind) = ERK_pp(end)/250;

end
ERK_times_CD3z = ERK_times;
ERK_concs_CD3z = ERK_concs;
"Finished with Cd3z..."

%% Only the Grb2 effect 
parameters(5)  = 500;  %Add CD28
parameters(64) = 0;    %Take away GADS
parameters(66) = 0.1;  %Add Grb2
parfor ind=1:sim_size
   
   params     = parameters;
   params(indices) = param_samples(ind,:);
   try
       [~, tpts, ~, observables_out] = model_func(timepoints, params);
   catch
       sprintf('Integration error occured at ind = %f', ind)
       ERK_times(ind) = NaN;
       ERK_concs(ind) = NaN;
       continue
   end
   
   ERK_pp = observables_out(:,133);
   halftime_index_curr = sum(ERK_pp<125);
   ERK_times(ind) = timepoints(halftime_index_curr);
   ERK_concs(ind) = ERK_pp(end)/250;

end
ERK_times_Grb2 = ERK_times;
ERK_concs_Grb2 = ERK_concs;
"Finished with Grb2..."

%% Only the GADS effect 
parameters(5)  = 500;  %Add CD28
parameters(64) = 0.1;  %Add GADS
parameters(66) = 0;    %Take away Grb2
parfor ind=1:sim_size
   
   params     = parameters;
   params(indices) = param_samples(ind,:);
   try
       [~, tpts, ~, observables_out] = model_func(timepoints, params);
   catch
       sprintf('Integration error occured at ind = %f', ind)
       ERK_times(ind) = NaN;
       ERK_concs(ind) = NaN;
       continue
   end
   
   ERK_pp = observables_out(:,133);
   halftime_index_curr = sum(ERK_pp<125);
   ERK_times(ind) = timepoints(halftime_index_curr);
   ERK_concs(ind) = ERK_pp(end)/250;

end
ERK_times_GADS = ERK_times;
ERK_concs_GADS = ERK_concs;
"Finished with GADS..."


%% Only the Kinetic effect 
parameters(5) = 500; %Add CD28
parameters(64)= 0;   %Take away GADS
parameters(66)= 0;   %Take away Grb2
mu(3:8) = [320, 270, 188, 352, 0.1153, 1170]; %Add Kinetic Effects
for ind = 3:8
    
    pdf = makedist('Normal','mu',mu(ind),'sigma',mu(ind)/3);
    pdf = truncate(pdf, 0, 2*mu(ind));
    param_values = random(pdf, sim_size, 1);
    param_samples(:,ind) = param_values;
    
end

parfor ind=1:sim_size
    
    params     = parameters;
    params(indices) = param_samples(ind,:);
    try
        [~, tpts, ~, observables_out] = model_func(timepoints, params);
    catch
        sprintf('Integration error occured at ind = %f', ind)
        ERK_times(ind) = NaN;
        ERK_concs(ind) = NaN;
        continue
    end
    
    ERK_pp = observables_out(:,133);
    halftime_index_curr = sum(ERK_pp<125);
    ERK_times(ind) = timepoints(halftime_index_curr);
    ERK_concs(ind) = ERK_pp(end)/250;
    
end
ERK_times_Kin = ERK_times;
ERK_concs_Kin = ERK_concs;
"Finished with all of Kin!!!"
save('Kinetic_varied_low_breakdown.mat','ERK_times_CD3z','ERK_concs_CD3z', ...
            'ERK_times_Grb2','ERK_concs_Grb2','ERK_times_GADS','ERK_concs_GADS', ...
               'ERK_times_Kin','ERK_concs_Kin')