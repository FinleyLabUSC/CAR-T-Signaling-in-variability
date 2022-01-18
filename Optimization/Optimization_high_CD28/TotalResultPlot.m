close all
clear all
clc

p0 = [1170, 0.1, 20, 0.9634, 67.8442];

load('PSOresults_01.mat')
paramfit(:,:,1) = gp_aggregator(1:11,:);
paramfit(:,:,1) = log10(paramfit(:,:,1)./p0);

load('PSOresults_1.mat')
paramfit(:,:,2) = gp_aggregator(1:11,:);
paramfit(:,:,2) = log10(paramfit(:,:,2)./p0);

load('PSOresults_5.mat')
paramfit(:,:,3) = gp_aggregator(1:11,:);
paramfit(:,:,3) = log10(paramfit(:,:,3)./p0);

load('PSOresults_10.mat')
paramfit(:,:,4) = gp_aggregator(1:11,:);
paramfit(:,:,4) = log10(paramfit(:,:,4)./p0);

param_names = ["Kcat_{LCKPU\_CD3z}", "CSK_{on}", "Kcat_{ZAP}", "Kcat_{CD45\_LCK505}", "Kcat_{CD45\_A1}"];
x = [-1, 0, 1];
color = ['b', 'r', 'g', 'c', 'm'];

for idx=1:5
    
    scatter( repmat(-1+(idx-3)/11,1,11), paramfit(:,idx,1), 150, 'MarkerFaceColor', color(idx), ...
                                         'MarkerFaceAlpha', .2, 'MarkerEdgeColor','none')
    hold on
    
end

for idx=1:5
    
    scatter( repmat(0+(idx-3)/11,1,11), paramfit(:,idx,2), 150, 'MarkerFaceColor', color(idx), ...
                                          'MarkerFaceAlpha', .2,'MarkerEdgeColor','none')
    hold on
    
end


for idx=1:5
    
    scatter( repmat(1+(idx-3)/11,1,11), paramfit(:,idx,4), 150, 'MarkerFaceColor', color(idx),  ...
                                           'MarkerFaceAlpha', .2,'MarkerEdgeColor','none')
    hold on
    
end

ylim([0 2.5])
xlim([-1.5 +1.5])
xticks([-1 0 +1])
xticklabels(["-1", "0", "+1"])
set(gca,'XTickLabel', get(gca,'XTickLabel'),'fontsize',48,'FontWeight','bold')
set(gca,'XTickLabel', get(gca,'XTickLabel'),'fontsize',48,'FontWeight','bold')
ylabel("LOG_{10}(fold change)", 'FontWeight', 'Bold','fontsize',52)
xlabel("LOG_{10}(\lambda)", 'FontWeight', 'Bold','fontsize',52)
legend(param_names)

