MPCTprofile=parallel.cluster.Local;
MPCTprofile.NumWorkers=24;
parpool(MPCTprofile, 23);

parameters = [ 45, 1000, 250, 500, 500, 870, 540, 600, 300, 120, 25, 800, ...
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

sim_size = 1e5;
rng(12121995);
ERK_times  = zeros(size(sim_size,1));
ERK_concs  = zeros(size(sim_size,1));
timepoints = 0:0.05:30;
timepoints = timepoints';
parameters(42:49) = [167, 320, 270, 188, 352, 205, 0.1153, 1170];
parameters(61) = 65.20;

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
ERK_times_CD28 = ERK_times;
ERK_concs_CD28 = ERK_concs;
param_samples_CD28 = param_samples;
"Finished with half of Kinetic_varied..."

parameters(5) = 0;
mu(3:8) = [331.68, 270, 153.82, 484.9, 0.1153, 360];

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
ERK_times_CD3z = ERK_times;
ERK_concs_CD3z = ERK_concs;
param_samples_CD3z = param_samples;
"Finished with all of Kinetic_varied!!!"

try
	save('Kinetic_varied_1e5_high.mat', 'ERK_times_CD28', 'ERK_concs_CD28', 'param_samples_CD28')
	save('Kinetic_varied_1e5_high.mat', 'ERK_times_CD3z', 'ERK_concs_CD3z', 'param_samples_CD3z','-append')
catch
	"Issue When Saving"
	try 
		xlswrite("times_cd28.xlsx",ERK_times_CD28)
        xlswrite("concs_cd28.xlsx",ERK_concs_CD28)
        xlswrite("times_cd3z.xlsx",ERK_times_CD3z)
        xlswrite("concs_cd3z.xlsx",ERK_times_CD3z)
        xlswrite("param_samples_CD28.xlsx", param_samples_CD28)
        xlswrite("param_samples_CD3z.xlsx", param_samples_CD3z)
    catch
        "Issue When Writing Excel"
    end

end

