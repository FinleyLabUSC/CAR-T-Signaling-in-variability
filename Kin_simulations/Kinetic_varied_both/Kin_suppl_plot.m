load("Kinetic_varied\Kinetic_varied_high_1e5.mat")
CD28_def=ERK_times_CD28;
CD3z_def=ERK_times_CD3z;
load("Kinetic_varied_zap\Kinetic_varied_1e5_high.mat")
CD28_zap=ERK_times_CD28;
CD3z_zap=ERK_times_CD3z;
load("Kinetic_varied_lck\Kinetic_varied_1e5_high.mat")
CD28_lck=ERK_times_CD28;
CD3z_lck=ERK_times_CD3z;
load("Kinetic_varied_both\Kinetic_varied_1e5_high.mat")
CD28_both=ERK_times_CD28;
CD3z_both=ERK_times_CD3z;

figure(1)
histogram(CD28_def,0:0.5:30)
hold on
histogram(CD28_zap,0:0.5:30)
histogram(CD28_lck,0:0.5:30)
histogram(CD28_both,0:0.5:30)
xlabel("Activation time, min", 'FontWeight', 'Bold')
ylabel("Number of Cells", 'FontWeight', 'Bold')
yticks([0 5e3 1e4 1.5e4 2e4 2.5e4 3e4 3.5e4 4e4 4.5e4 5e4])
yticklabels(["0", "5000", "10000", "15000", "20000", "25000" ,"30000", "35000", "40000", "45000", "50000"])

figure(2)
histogram(CD3z_def,0:0.5:30)
hold on
histogram(CD3z_zap,0:0.5:30)
histogram(CD3z_lck,0:0.5:30)
histogram(CD3z_both,0:0.5:30)
xlabel("Activation time, min", 'FontWeight', 'Bold')
ylabel("Number of Cells", 'FontWeight', 'Bold')
yticks([0 5e3 1e4 1.5e4 2e4 2.5e4 3e4 3.5e4 4e4 4.5e4 5e4])
yticklabels(["0", "5000", "10000", "15000", "20000", "25000" ,"30000", "35000", "40000", "45000", "50000"])


