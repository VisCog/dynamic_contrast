clear
%close all

vepFolder = [cd filesep 'fitdata_vep' filesep 'model_fits'];
psyFolder = [cd filesep 'fitdata_vep_psychophysics' filesep 'model_fits'];


vepFiles = dir([vepFolder filesep '*congruent.mat']);
psyFiles = dir([psyFolder filesep '*congruent.mat']);

% grouping variable
idx_psy_ns = find(cellfun(@(x) strcmpi(x(1:2),'NS'), {psyFiles.name}));
idx_psy_am = find(cellfun(@(x) strcmpi(x(1:2),'AM'), {psyFiles.name}));
idx_psy_bd = find(cellfun(@(x) strcmpi(x(1:2),'BD'), {psyFiles.name}));
grp_psy = nan(size(psyFiles));
%grp_psy(idx_psy_am) = 1;
%grp_psy(idx_psy_bd) = 2;
grp_psy(idx_psy_ns) = 3;

idx_vep_ns = find(cellfun(@(x) strcmpi(x(1:2),'NS'), {vepFiles.name}));
idx_vep_am = find(cellfun(@(x) strcmpi(x(1:2),'AM'), {vepFiles.name}));
idx_vep_bd = find(cellfun(@(x) strcmpi(x(1:2),'BD'), {vepFiles.name}));
grp_vep = nan(size(vepFiles));
%grp_vep(idx_psy_am) = 1;
%grp_vep(idx_psy_bd) = 2;
grp_vep(idx_psy_ns) = 3;

% quick alt to gathertable
all_psy = [];
for i = 1:length(psyFiles)
    load([psyFiles(i).folder filesep psyFiles(i).name]);
    all_psy = [all_psy p];
end
all_vep = [];
for i = 1:length(vepFiles)
    load([vepFiles(i).folder filesep vepFiles(i).name]);
    all_vep = [all_vep p];
end

ymaxval = [all_vep.meanModelErr all_psy.meanModelErr ...
    all_vep.maxModelErr all_psy.maxModelErr ...
    all_vep.minkModelErr all_psy.minkModelErr];
ymaxidx = [~isnan(grp_vep) ~isnan(grp_psy) ...
    ~isnan(grp_vep) ~isnan(grp_psy)...
    ~isnan(grp_vep) ~isnan(grp_psy)];
ymax = max(ymaxval(ymaxidx));
ymax = ceil(ymax*100)/100;

figure(1); clf; set(gcf, 'Name', 'Model error')

subplot(1,3,1); hold on; grid on;
quickPlot([all_vep.meanModelErr], [all_psy.meanModelErr], grp_vep, grp_psy, 0)
ylabel('Error: mean model')
title('mean model')
ylim([0 ymax])
set(gca, 'FontSize', 12)

subplot(1,3,2); hold on; grid on;
quickPlot([all_vep.maxModelErr], [all_psy.maxModelErr], grp_vep, grp_psy, 0)
ylim([0 ymax])
ylabel('Error: max model')
title('max model')
set(gca, 'FontSize', 12)

subplot(1,3,3); hold on; grid on;
quickPlot([all_vep.minkModelErr], [all_psy.minkModelErr], grp_vep, grp_psy, 0)
ylim([0 ymax])
ylabel('Error: minkowski model')
title('minkowski model')
set(gca, 'FontSize', 12)




figure(2); clf; set(gcf, 'Name', 'free parameters')
hold on; grid on;
quickPlot([all_vep.n], [all_psy.n], grp_vep, grp_psy, 0)
ylabel('minkowski parameter')
set(gca, 'FontSize', 12)
title('minkowski parameter')
% crit = 40; % outlier crit for plot
% outlr = [all_psy.n]>crit;
% ylim([0 crit]);
% tmp = [all_psy.n]; 
% text(2, crit+1, [num2str(sum(outlr)) ' values above ' num2str(crit) ' not shown'], 'HorizontalAlignment', 'center')
set(gca, 'YScale', 'log')

function quickPlot(vepvector, psyvector, vepgroup, psygroup, legendtf)
c = ['r', 'g', 'b'];
names = {'AM', 'BD', 'NS'};

tmp = vepgroup;
for g = unique(tmp(~isnan(tmp)))'
    scatter(1, vepvector(vepgroup == g), c(g), 'jitter', 0.2);
end

tmp = psygroup;
for g = unique(tmp(~isnan(tmp)))'
    scatter(2, psyvector(psygroup == g), c(g), 'jitter', 0.2);
end

xlim([0.5 2.5])
scatter([1 2], [mean(vepvector) mean(psyvector)], 50, 'k', 'filled');
text(1.1, mean(vepvector), [num2str(mean(vepvector),3)])
text(2.1, mean(psyvector), [num2str(mean(psyvector),3)])

xticks([1 2])
xticklabels({'VEP', 'PSYCHO'})
if legendtf == 1
    l1 = scatter(nan, nan, c(1), 'DisplayName', names{1});
    l2 = scatter(nan, nan, c(2), 'DisplayName', names{2});
    l3 = scatter(nan, nan, c(3), 'DisplayName', names{3});
    legend([l1 l2 l3])
end
end
