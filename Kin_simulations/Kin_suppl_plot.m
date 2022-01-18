load("Kinetic_varied_def\Kinetic_varied_low_1e5.mat")
CD28_def=ERK_times_CD28;
CD3z_def=ERK_times_CD3z;
load("Kinetic_varied_zap\Kinetic_varied_1e5_low.mat")
CD28_zap=ERK_times_CD28;
CD3z_zap=ERK_times_CD3z;
load("Kinetic_varied_lck\Kinetic_varied_1e5_low.mat")
CD28_lck=ERK_times_CD28;
CD3z_lck=ERK_times_CD3z;
load("Kinetic_varied_both\Kinetic_varied_1e5_low.mat")
CD28_both=ERK_times_CD28;
CD3z_both=ERK_times_CD3z;

figure(1)
histogram( CD28_def( (CD28_def<30) ),0:0.5:30, 'FaceColor', '#7E2F8E')
hold on
histogram( CD28_zap( (CD28_zap<30) ),0:0.5:30, 'FaceColor', '#EDB120')
histogram( CD28_lck( (CD28_lck<30) ),0:0.5:30, 'FaceColor', '#77AC30')
histogram( CD28_both( (CD28_both<30) ),0:0.5:30, 'FaceColor', '#A2142F')
xlabel("Activation time, min", 'FontWeight', 'Bold','FontSize',36)
ylabel("Number of Cells", 'FontWeight', 'Bold','FontSize',36)
yticks([0 5e3 1e4 1.5e4 2e4 2.5e4 3e4 3.5e4 4e4 4.5e4 5e4])
xlim([0 30])
ylim([0 4e4])
yticklabels(["0", "5000", "10000", "15000", "20000", "25000" ,"30000", "35000", "40000", "45000", "50000"])
set(gca,'XTickLabel', get(gca,'XTickLabel'),'fontsize',36,'FontWeight','bold')
set(gca,'YTickLabel', get(gca,'YTickLabel'),'fontsize',36,'FontWeight','bold')
legend(["Default Param.", "Optim. Kcat_{ZAP} ",  "Optim. Kcat_{LCKPU\_CD3z}", "Optim. Both"]);


figure(2)
histogram( CD3z_def( (CD3z_def<30) ),0:0.5:30, 'FaceColor', '#7E2F8E')
hold on
histogram( CD3z_zap( (CD3z_zap<30) ),0:0.5:30, 'FaceColor', '#EDB120')
histogram( CD3z_lck( (CD3z_lck<30) ),0:0.5:30, 'FaceColor', '#77AC30')
histogram( CD3z_both( (CD3z_both<30) ),0:0.5:30, 'FaceColor', '#A2142F')
xlabel("Activation time, min", 'FontWeight', 'Bold','FontSize',36)
ylabel("Number of Cells", 'FontWeight', 'Bold','FontSize',36)
xlim([0 30])
ylim([0 4e4])
yticks([0 5e3 1e4 1.5e4 2e4 2.5e4 3e4 3.5e4 4e4 4.5e4 5e4])
yticklabels(["0", "5000", "10000", "15000", "20000", "25000" ,"30000", "35000", "40000", "45000", "50000"])
set(gca,'XTickLabel', get(gca,'XTickLabel'),'fontsize',36,'FontWeight','bold')
set(gca,'YTickLabel', get(gca,'YTickLabel'),'fontsize',36,'FontWeight','bold')
legend(["Default Param.", "Optim. Kcat_{ZAP} ",  "Optim. Kcat_{LCKPU\_CD3z}", "Optim. Both"]);



figure(3)
x = categorical({'Active','Inactive'});
plt = bar(x, [ sum(CD28_def<30),  sum(CD28_def>=30);  
                  sum(CD28_zap<30), sum(CD28_zap>=30);
                    sum(CD28_lck<30), sum(CD28_lck>=30);
                        sum(CD28_both<30), sum(CD28_both>=30);]');
plt(1).FaceAlpha = 0.68;
plt(2).FaceAlpha = 0.68;
plt(3).FaceAlpha = 0.68;
plt(4).FaceAlpha = 0.68;
plt(1).FaceColor = '#7E2F8E';
plt(2).FaceColor = '#EDB120';
plt(3).FaceColor = '#77AC30';
plt(4).FaceColor = '#A2142F';
ylim([0 1.1e5])
yticks([0 2e4 4e4 6e4 8e4 10e4])
yticklabels(["0", "20000", "40000", "60000", "80000", "100000"])
ylabel("Number of Cells", 'FontWeight', 'Bold','FontSize',36)
set(gca,'YTickLabel', get(gca,'YTickLabel'),'fontsize',36,'FontWeight','bold')
set(gca,'XTickLabel', get(gca,'XTickLabel'),'fontsize',36,'FontWeight','bold')
legend(["Default Param.", "Optim. Kcat_{ZAP} ",  "Optim. Kcat_{LCKPU\_CD3z}", "Optim. Both"]);
xtickangle(90)

figure(4)
x = categorical({'Active','Inactive'});
plt = bar(x, [ sum(CD3z_def<30),  sum(CD3z_def>=30);  
                  sum(CD3z_zap<30), sum(CD3z_zap>=30);
                    sum(CD3z_lck<30), sum(CD3z_lck>=30);
                        sum(CD3z_both<30), sum(CD3z_both>=30);]');
plt(1).FaceAlpha = 0.68;
plt(2).FaceAlpha = 0.68;
plt(3).FaceAlpha = 0.68;
plt(4).FaceAlpha = 0.68;
plt(1).FaceColor = '#7E2F8E';
plt(2).FaceColor = '#EDB120';
plt(3).FaceColor = '#77AC30';
plt(4).FaceColor = '#A2142F';
ylim([0 1.1e5])
yticks([0 2e4 4e4 6e4 8e4 10e4])
yticklabels(["0", "20000", "40000", "60000", "80000", "100000"])
ylabel("Number of Cells", 'FontWeight', 'Bold','FontSize',36)
set(gca,'YTickLabel', get(gca,'YTickLabel'),'fontsize',36,'FontWeight','bold')
set(gca,'XTickLabel', get(gca,'XTickLabel'),'fontsize',36,'FontWeight','bold')
legend(["Default Param.", "Kcat_{ZAP} Optim.",  "Kcat_{LCKPU\_CD3z} Optim.", "Both Optim."]);
xtickangle(90)
legend(["Default Param.", "Optim. Kcat_{ZAP} ",  "Optim. Kcat_{LCKPU\_CD3z}", "Optim. Both"]);

