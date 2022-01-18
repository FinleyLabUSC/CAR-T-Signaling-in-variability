close all
clear all
clc

load("Ant_varied_results.mat")

figure(1)
histogram(ERK_times_CD3z( (ERK_times_CD3z<30) ), 0:0.5:30)
hold on 
histogram(ERK_times_CD28( (ERK_times_CD28<30) ), 0:0.5:30)
ylim([0 4e4])
xlim([0 30])
yticks([0 5000 10000 15000 20000 25000 30000 35000 40000])
yticklabels(["0", "5000", "10000", "15000", "20000", "25000", "30000", "35000", "40000"])
ylabel("Number of Cells", 'FontWeight', 'Bold','fontsize',28)
xlabel("Activation time, min", 'FontWeight', 'Bold','fontsize',28)
set(gca,'XTickLabel', get(gca,'XTickLabel'),'fontsize',28,'FontWeight','bold')
set(gca,'YTickLabel', get(gca,'YTickLabel'),'fontsize',28,'FontWeight','bold')


figure(2)
x = categorical({'Active','Inactive'});
plt = bar(x, [ sum(ERK_times_CD3z<30),  sum(ERK_times_CD3z>=30);  
                  sum(ERK_times_CD28<30), sum(ERK_times_CD28>=30)]');
plt(1).FaceAlpha = 0.68;
plt(2).FaceAlpha = 0.68;
ylim([0 1.1e5])
yticks([0 2e4 4e4 6e4 8e4 10e4])
yticklabels(["0", "20000", "40000", "60000", "80000", "100000"])
ylabel("Number of Cells", 'FontWeight', 'Bold','fontsize',28)
set(gca,'YTickLabel', get(gca,'YTickLabel'),'fontsize',28,'FontWeight','bold')
set(gca,'XTickLabel', get(gca,'XTickLabel'),'fontsize',28,'FontWeight','bold')
set(gca,'YTickLabel', get(gca,'YTickLabel'),'fontsize',28,'FontWeight','bold')
xtickangle(90)
