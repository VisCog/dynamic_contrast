clear
%close all

vepFolder = [cd filesep 'fitdata_vep' filesep 'model_fits'];
psyFolder = [cd filesep 'fitdata_vep_psychophysics' filesep 'model_fits'];


%vepFiles = dir([vepFolder filesep '*congruent.mat']);
%psyFiles = dir([psyFolder filesep '*congruent.mat']);
vepFiles = dir([vepFolder filesep 'NS_*congruent.mat']);
psyFiles = dir([psyFolder filesep 'NS_*congruent.mat']);


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
    clear p
end
all_vep = [];
for i = 1:length(vepFiles)
    load([vepFiles(i).folder filesep vepFiles(i).name]);
    all_vep = [all_vep p];
    clear p
end

% vepFiles_orth = dir([vepFolder filesep '*orthogonal.mat']);
% psyFiles_orth = dir([psyFolder filesep '*orthogonal.mat']);
% all_psy_orth = [];
% for i = 1:length(psyFiles_orth)
%     load([psyFiles_orth(i).folder filesep psyFiles_orth(i).name]);
%     all_psy_orth = [all_psy_orth p];
%     clear p
% end
% all_vep_orth = [];
% for i = 1:length(vepFiles_orth)
%     load([vepFiles_orth(i).folder filesep vepFiles_orth(i).name]);
%     all_vep_orth = [all_vep_orth p];
%     clear p
% end

ymaxval = [...
    all_vep.meanModelErr...
    all_psy.meanModelErr ...
    all_vep.maxModelErr...
    all_psy.maxModelErr ...
    all_vep.minkModelErr...
    all_psy.minkModelErr];

ymaxidx = [~isnan(grp_vep)...
    ~isnan(grp_psy) ...
    ~isnan(grp_vep) ...
    ~isnan(grp_psy)...
    ~isnan(grp_vep) ...
    ~isnan(grp_psy)];

ymax = max(ymaxval);
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
set(gca, 'YScale', 'log')

tmp = [all_vep.n];
vep_n = tmp(grp_vep==3);
disp(['Geometric mean for controls (n = ' num2str(sum(grp_vep==3)) '), VEP: ' num2str(geomean(vep_n),4)])
tmp = [all_psy.n];
psy_n = tmp(grp_psy==3);
disp(['Geometric mean for controls (n = ' num2str(sum(grp_psy==3)) '), psychophysics: ' num2str(geomean(psy_n),4)])


disp(['Median for controls (n = ' num2str(sum(grp_vep==3)) '), VEP: ' num2str(median(vep_n),4)])
disp(['Median for controls (n = ' num2str(sum(grp_psy==3)) '), psychophysics: ' num2str(median(psy_n),4)])


log_vep_n = log(vep_n);
log_psy_n = log(psy_n);
[h,p,ci,stats] = ttest(log_vep_n, log_psy_n);
disp(['Paired t-test on log-transformed values: t(' num2str(stats.df) ') = ' num2str(abs(stats.tstat),3) ', p = ' num2str(p,3)])
[p,h,stats] = signrank(vep_n, psy_n,'method','approximate');
disp(['Paired sign-rank test on original values: Z = ' num2str(abs(stats.zval),3) ', p = ' num2str(p,3)])



% mean/max weighting
figure(3); clf; set(gcf, 'Name', 'Mean/max model weighting')

subplot(1,2,1); hold on; grid on;
quickPlot([all_vep.mnmxwghtModelErr], [all_psy.mnmxwghtModelErr], grp_vep, grp_psy, 0)
ylabel('Error: mean/max weighting model')
title('model error')
%ylim([0 max([all_vep.mnmxwghtModelErr all_psy.mnmxwghtModelErr])])
set(gca, 'FontSize', 12)

subplot(1,2,2); hold on; grid on;
quickPlot([all_vep.w], [all_psy.w], grp_vep, grp_psy, 0)
ylabel('mean vs max weight')
set(gca, 'FontSize', 12)
title('model weighting')
%set(gca, 'YScale', 'log')



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
plotMean(vepvector, psyvector, vepgroup, psygroup, 'median')

xticks([1 2])
xticklabels({'VEP', 'PSYCHO'})
if legendtf == 1
    l1 = scatter(nan, nan, c(1), 'DisplayName', names{1});
    l2 = scatter(nan, nan, c(2), 'DisplayName', names{2});
    l3 = scatter(nan, nan, c(3), 'DisplayName', names{3});
    legend([l1 l2 l3])
end
end

function plotMean(vepvector, psyvector, vepgroup, psygroup, type)
% hard-coded as controls only
switch type
    case 'arithmetic'
        tmpvep = mean(vepvector(vepgroup == 3));
        tmppsy = mean(psyvector(psygroup == 3));

    case 'geometric'
        tmpvep = geomean(vepvector(vepgroup == 3));
        tmppsy = geomean(psyvector(psygroup == 3));

    case 'harmonic'
        tmpvep = harmmean(vepvector(vepgroup == 3));
        tmppsy = harmmean(psyvector(psygroup == 3));
    case 'median'
        tmpvep = median(vepvector(vepgroup == 3));
        tmppsy = median(psyvector(psygroup == 3));
end
scatter([1 2], [tmpvep tmppsy], 50, 'k', 'filled');
text(1.1, tmpvep, [num2str(tmpvep,3)])
text(2.1, tmppsy, [num2str(tmppsy,3)])
end
