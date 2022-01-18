clear all
clc

%% Dose response script
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



ant_concs = logspace(-2, log10(45), 50);
timepoints = 0:0.05:30;
timepoints = timepoints';


% This is for CD28-CAR
parameters(42:49) = [167, 320, 270, 188, 352, 205, 0.1153, 1170];
for ind=1:length(ant_concs)
    
    parameters(1) = ant_concs(ind);
    [~, tpts, ~, observables_out] = model_func(timepoints, parameters);
    
    ERK_pp = observables_out(:,133);
    halftime_index_curr=sum(ERK_pp<125);
    ERK_times_CD28(ind)=timepoints(halftime_index_curr);
     
end

% This is for CD3z-CAR
parameters(5) = 0;
parameters(42:49) = [96.548, 331.68, 270, 153.82, 484.9, 405.43, 0.1153, 360];
for ind=1:length(ant_concs)
    
    parameters(1) = ant_concs(ind);
    [~, tpts, ~, observables_out] = model_func(timepoints, parameters);
    
    ERK_pp = observables_out(:,133);
    halftime_index_curr=sum(ERK_pp<125);
    ERK_times_CD3z(ind)=timepoints(halftime_index_curr);
     
end




figure(1)
semilogx(ant_concs, ERK_times_CD3z,'linewidth',7)
hold on 
semilogx(ant_concs, ERK_times_CD28,'linewidth',7)
ylim([0 32])
xlim([0 45])
yticks([0 10 20 30])
yticklabels(["0", "10", "20", "30"])
ylabel("Activation time, min", 'FontWeight', 'Bold','fontsize',26)
xlabel("Antigen concentration, molec. per \mu^2", 'FontWeight', 'Bold','fontsize',26)
set(gca,'XTickLabel', get(gca,'XTickLabel'),'fontsize',26,'FontWeight','bold')
set(gca,'YTickLabel', get(gca,'YTickLabel'),'fontsize',26,'FontWeight','bold')
