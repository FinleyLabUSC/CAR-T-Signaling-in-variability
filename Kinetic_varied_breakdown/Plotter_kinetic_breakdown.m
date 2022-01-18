close all
clear all
clc


%% Histogram figure

load("Kinetic_varied_high_breakdown.mat")

figure(1)
histogram(ERK_times_CD3z( (ERK_times_CD3z<30) ), 0:0.5:30)
hold on 
histogram(ERK_times_Grb2( (ERK_times_Grb2<30) ), 0:0.5:30)
ylim([0 3e3])
xlim([0 30])
yticks([0 500 1000 1500 2000 2500 3000])
yticklabels(["0", "500", "1000", "1500", "2000", "2500", "3000"])
ylabel("Number of Cells", 'FontWeight', 'Bold','fontsize',30)
xlabel("Activation time, min", 'FontWeight', 'Bold','fontsize',30)
set(gca,'XTickLabel', get(gca,'XTickLabel'),'fontsize',30,'FontWeight','bold')
set(gca,'YTickLabel', get(gca,'YTickLabel'),'fontsize',30,'FontWeight','bold')


figure
histogram(ERK_times_CD3z( (ERK_times_CD3z<30) ), 0:0.5:30)
hold on 
histogram(ERK_times_GADS( (ERK_times_GADS<30) ), 0:0.5:30)
ylim([0 3e3])
xlim([0 30])
yticks([0 500 1000 1500 2000 2500 3000])
yticklabels(["0", "500", "1000", "1500", "2000", "2500", "3000"])
ylabel("Number of Cells", 'FontWeight', 'Bold','fontsize',30)
xlabel("Activation time, min", 'FontWeight', 'Bold','fontsize',30)
set(gca,'XTickLabel', get(gca,'XTickLabel'),'fontsize',30,'FontWeight','bold')
set(gca,'YTickLabel', get(gca,'YTickLabel'),'fontsize',30,'FontWeight','bold')



figure
histogram(ERK_times_CD3z( (ERK_times_CD3z<30) ), 0:0.5:30)
hold on 
histogram(ERK_times_Kin( (ERK_times_Kin<30) ), 0:0.5:30)
ylim([0 3e3])
xlim([0 30])
yticks([0 500 1000 1500 2000 2500 3000])
yticklabels(["0", "500", "1000", "1500", "2000", "2500", "3000"])
ylabel("Number of Cells", 'FontWeight', 'Bold','fontsize',30)
xlabel("Activation time, min", 'FontWeight', 'Bold','fontsize',30)
set(gca,'XTickLabel', get(gca,'XTickLabel'),'fontsize',30,'FontWeight','bold')
set(gca,'YTickLabel', get(gca,'YTickLabel'),'fontsize',30,'FontWeight','bold')



load("Kinetic_varied_low_breakdown.mat")

figure
histogram(ERK_times_CD3z( (ERK_times_CD3z<30) ), 0:0.5:30)
hold on 
histogram(ERK_times_Grb2( (ERK_times_Grb2<30) ), 0:0.5:30)
ylim([0 3e3])
xlim([0 30])
yticks([0 500 1000 1500 2000 2500 3000])
yticklabels(["0", "500", "1000", "1500", "2000", "2500", "3000"])
ylabel("Number of Cells", 'FontWeight', 'Bold','fontsize',30)
xlabel("Activation time, min", 'FontWeight', 'Bold','fontsize',30)
set(gca,'XTickLabel', get(gca,'XTickLabel'),'fontsize',30,'FontWeight','bold')
set(gca,'YTickLabel', get(gca,'YTickLabel'),'fontsize',30,'FontWeight','bold')



figure
histogram(ERK_times_CD3z( (ERK_times_CD3z<30) ), 0:0.5:30)
hold on 
histogram(ERK_times_GADS( (ERK_times_GADS<30) ), 0:0.5:30)
ylim([0 3e3])
xlim([0 30])
yticks([0 500 1000 1500 2000 2500 3000])
yticklabels(["0", "500", "1000", "1500", "2000", "2500", "3000"])
ylabel("Number of Cells", 'FontWeight', 'Bold','fontsize',30)
xlabel("Activation time, min", 'FontWeight', 'Bold','fontsize',30)
set(gca,'XTickLabel', get(gca,'XTickLabel'),'fontsize',30,'FontWeight','bold')
set(gca,'YTickLabel', get(gca,'YTickLabel'),'fontsize',30,'FontWeight','bold')


figure
histogram(ERK_times_CD3z( (ERK_times_CD3z<30) ), 0:0.5:30)
hold on 
histogram(ERK_times_Kin( (ERK_times_Kin<30) ), 0:0.5:30)
ylim([0 3e3])
xlim([0 30])
yticks([0 500 1000 1500 2000 2500 3000])
yticklabels(["0", "500", "1000", "1500", "2000", "2500", "3000"])
ylabel("Number of Cells", 'FontWeight', 'Bold','fontsize',30)
xlabel("Activation time, min", 'FontWeight', 'Bold','fontsize',30)
set(gca,'XTickLabel', get(gca,'XTickLabel'),'fontsize',30,'FontWeight','bold')
set(gca,'YTickLabel', get(gca,'YTickLabel'),'fontsize',30,'FontWeight','bold')





%% Bar figure

load("Kinetic_varied_high_breakdown.mat")

figure
x = categorical({'Active','Inactive'});
plt = bar(x, [ sum(ERK_times_CD3z<30),  sum(ERK_times_CD3z>=30);  
                  sum(ERK_times_Grb2<30), sum(ERK_times_Grb2>=30)]');
plt(1).FaceAlpha = 0.68;
plt(2).FaceAlpha = 0.68;
ylim([0 1.1e4])
yticks([0 2e3 4e3 6e3 8e3 10e3])
yticklabels(["0", "2000", "4000", "6000", "8000", "10000"])
ylabel("Number of Cells", 'FontWeight', 'Bold','fontsize',30)
set(gca,'YTickLabel', get(gca,'YTickLabel'),'fontsize',30,'FontWeight','bold')
xtickangle(90)

figure
x = categorical({'Active','Inactive'});
plt = bar(x, [ sum(ERK_times_CD3z<30),  sum(ERK_times_CD3z>=30);  
                  sum(ERK_times_GADS<30), sum(ERK_times_GADS>=30)]');
plt(1).FaceAlpha = 0.68;
plt(2).FaceAlpha = 0.68;
ylim([0 1.1e4])
yticks([0 2e3 4e3 6e3 8e3 10e3])
yticklabels(["0", "2000", "4000", "6000", "8000", "10000"])
ylabel("Number of Cells", 'FontWeight', 'Bold','fontsize',30)
set(gca,'YTickLabel', get(gca,'YTickLabel'),'fontsize',30,'FontWeight','bold')
xtickangle(90)


figure
x = categorical({'Active','Inactive'});
plt = bar(x, [ sum(ERK_times_CD3z<30),  sum(ERK_times_CD3z>=30);  
                  sum(ERK_times_Kin<30), sum(ERK_times_Kin>=30)]');
plt(1).FaceAlpha = 0.68;
plt(2).FaceAlpha = 0.68;
ylim([0 1.1e4])
yticks([0 2e3 4e3 6e3 8e3 10e3])
yticklabels(["0", "2000", "4000", "6000", "8000", "10000"])
ylabel("Number of Cells", 'FontWeight', 'Bold','fontsize',30)
set(gca,'YTickLabel', get(gca,'YTickLabel'),'fontsize',30,'FontWeight','bold')
xtickangle(90)


load("Kinetic_varied_low_breakdown.mat")

figure
x = categorical({'Active','Inactive'});
plt = bar(x, [ sum(ERK_times_CD3z<30),  sum(ERK_times_CD3z>=30);  
                  sum(ERK_times_Grb2<30), sum(ERK_times_Grb2>=30)]');
plt(1).FaceAlpha = 0.68;
plt(2).FaceAlpha = 0.68;
ylim([0 1.1e4])
yticks([0 2e3 4e3 6e3 8e3 10e3])
yticklabels(["0", "2000", "4000", "6000", "8000", "10000"])
ylabel("Number of Cells", 'FontWeight', 'Bold','fontsize',30)
set(gca,'YTickLabel', get(gca,'YTickLabel'),'fontsize',30,'FontWeight','bold')
xtickangle(90)


figure
x = categorical({'Active','Inactive'});
plt = bar(x, [ sum(ERK_times_CD3z<30),  sum(ERK_times_CD3z>=30);  
                  sum(ERK_times_GADS<30), sum(ERK_times_GADS>=30)]');
plt(1).FaceAlpha = 0.68;
plt(2).FaceAlpha = 0.68;
ylim([0 1.1e4])
yticks([0 2e3 4e3 6e3 8e3 10e3])
yticklabels(["0", "2000", "4000", "6000", "8000", "10000"])
ylabel("Number of Cells", 'FontWeight', 'Bold','fontsize',30)
set(gca,'YTickLabel', get(gca,'YTickLabel'),'fontsize',30,'FontWeight','bold')
xtickangle(90)


figure
x = categorical({'Active','Inactive'});
plt = bar(x, [ sum(ERK_times_CD3z<30),  sum(ERK_times_CD3z>=30);  
                  sum(ERK_times_Kin<30), sum(ERK_times_Kin>=30)]');
plt(1).FaceAlpha = 0.68;
plt(2).FaceAlpha = 0.68;
ylim([0 1.1e4])
yticks([0 2e3 4e3 6e3 8e3 10e3])
yticklabels(["0", "2000", "4000", "6000", "8000", "10000"])
ylabel("Number of Cells", 'FontWeight', 'Bold','fontsize',30)
set(gca,'YTickLabel', get(gca,'YTickLabel'),'fontsize',30,'FontWeight','bold')
xtickangle(90)
