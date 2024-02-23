clear

sGroup = 'AM';
sCondition = 'congruent';

%%% load files
vepFolder = [cd filesep 'fitdata_vep' filesep 'model_fits'];     % vep
psyFolder = [cd filesep 'fitdata_vep_psychophysics' filesep 'model_fits']; %psychophysics

vepFiles = dir([vepFolder filesep sGroup '_*' sCondition '.mat']);  % grab sGroup
psyFiles = dir([psyFolder filesep sGroup '_*' sCondition '.mat']);

% % grouping variable
% idx_psy_ns = find(cellfun(@(x) strcmpi(x(1:2),'NS'), {psyFiles.name}));
% idx_psy_am = find(cellfun(@(x) strcmpi(x(1:2),'AM'), {psyFiles.name}));
% idx_psy_bd = find(cellfun(@(x) strcmpi(x(1:2),'BD'), {psyFiles.name}));
% grp_psy = nan(size(psyFiles));
% %grp_psy(idx_psy_am) = 1;
% %grp_psy(idx_psy_bd) = 2;
% grp_psy(idx_psy_ns) = 3; % grabbing NS only right now
% 
% idx_vep_ns = find(cellfun(@(x) strcmpi(x(1:2),'NS'), {vepFiles.name}));
% idx_vep_am = find(cellfun(@(x) strcmpi(x(1:2),'AM'), {vepFiles.name}));
% idx_vep_bd = find(cellfun(@(x) strcmpi(x(1:2),'BD'), {vepFiles.name}));
% grp_vep = nan(size(vepFiles));
% %grp_vep(idx_psy_am) = 1;
% %grp_vep(idx_psy_bd) = 2;
% grp_vep(idx_psy_ns) = 3;

grp_psy = ones(size(psyFiles))*3; % legacy index
grp_vep = ones(size(vepFiles))*3;

% load fits and compile
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

% plot lims
ymaxval = [...
    all_vep.meanModel_err...
    all_psy.meanModel_err ...
    all_vep.maxModel_err...
    all_psy.maxModel_err];
%     all_vep.minkModelErr...
%     all_psy.minkModelErr];
ymax = max(ymaxval);
ymax = ceil(ymax*100)/100;


%% Fig 1: Model error
figure(1); clf; set(gcf, 'Name', ['Model error (' sGroup ': ' sCondition ')']);
nplots = 4;

%%%%% plot: error for mean model fit

subplot(1,nplots,1); hold on; grid on; 
quickPlot([all_vep.meanModel_err], [all_psy.meanModel_err], grp_vep, grp_psy, 0);
ylabel('MSE: mean model'); title('mean model');
ylim([0 ymax]); set(gca, 'FontSize', 12);

% label unusual values with sID in case want to investigate
[ovals, osids] = checkoutlier([all_vep.meanModel_err], {all_vep.sID}); 
for o = 1:length(ovals), text(1.1, ovals(o), osids{o}, 'Interpreter','none'), end
[ovals, osids] = checkoutlier([all_psy.meanModel_err], {all_psy.sID});
for o = 1:length(ovals), text(2.1, ovals(o), osids{o}, 'Interpreter','none'), end

%%%%% plot: error for max model fit

subplot(1,nplots,2); hold on; grid on;
quickPlot([all_vep.maxModel_err], [all_psy.maxModel_err], grp_vep, grp_psy, 0);
ylabel('MSE: max model'); title('max model');
ylim([0 ymax]); set(gca, 'FontSize', 12);

[ovals, osids] = checkoutlier([all_vep.maxModel_err], {all_vep.sID});
for o = 1:length(ovals), text(1.1, ovals(o), osids{o}, 'Interpreter','none'), end
[ovals, osids] = checkoutlier([all_psy.maxModel_err], {all_psy.sID});
for o = 1:length(ovals), text(2.1, ovals(o), osids{o}, 'Interpreter','none'), end

%%%%% plot 3: error for minkowski model fit
 
subplot(1,nplots,3); hold on; grid on;
% quickPlot([all_vep.minkModelErr], [all_psy.minkModelErr], grp_vep, grp_psy, 0);
% ylabel('MSE: minkowski model'); title('minkowski model')
% ylim([0 ymax]); set(gca, 'FontSize', 12);
% 
% [ovals, osids] = checkoutlier([all_vep.minkModelErr], {all_vep.sID});
% for o = 1:length(ovals), text(1.1, ovals(o), osids{o}, 'Interpreter','none'), end
% [ovals, osids] = checkoutlier([all_psy.minkModelErr], {all_psy.sID});
% for o = 1:length(ovals), text(2.1, ovals(o), osids{o}, 'Interpreter','none'), end

quickPlot([all_vep.meanModel_wghtd_err], [all_psy.meanModel_wghtd_err], grp_vep, grp_psy, 0);
ylabel('MSE: weighted mean model'); title('weighted (linear attn) mean model')
ylim([0 ymax]); set(gca, 'FontSize', 12);

[ovals, osids] = checkoutlier([all_vep.meanModel_wghtd_err], {all_vep.sID});
for o = 1:length(ovals), text(1.1, ovals(o), osids{o}, 'Interpreter','none'), end
[ovals, osids] = checkoutlier([all_psy.meanModel_wghtd_err], {all_psy.sID});
for o = 1:length(ovals), text(2.1, ovals(o), osids{o}, 'Interpreter','none'), end

%%%%% plot 4: error for mixture (mean/max) model fit

subplot(1,nplots,4); hold on; grid on;
% quickPlot([all_vep.mnmxwghtModelErr], [all_psy.mnmxwghtModelErr], grp_vep, grp_psy, 0);
% ylabel('MSE: mixture (mean/max) model'); title('mixture (mean/max) model')
% ylim([0 ymax]); set(gca, 'FontSize', 12);
% 
% [ovals, osids] = checkoutlier([all_vep.mnmxwghtModelErr], {all_vep.sID});
% for o = 1:length(ovals), text(1.1, ovals(o), osids{o}, 'Interpreter','none'), end
% [ovals, osids] = checkoutlier([all_psy.mnmxwghtModelErr], {all_psy.sID});
% for o = 1:length(ovals), text(2.1, ovals(o), osids{o}, 'Interpreter','none'), end

quickPlot([all_vep.maxModel_wghtd_err], [all_psy.maxModel_wghtd_err], grp_vep, grp_psy, 0);
ylabel('MSE: weighted max model'); title('weighted (linear attn) max model')
ylim([0 ymax]); set(gca, 'FontSize', 12);

[ovals, osids] = checkoutlier([all_vep.maxModel_wghtd_err], {all_vep.sID});
for o = 1:length(ovals), text(1.1, ovals(o), osids{o}, 'Interpreter','none'), end
[ovals, osids] = checkoutlier([all_psy.maxModel_wghtd_err], {all_psy.sID});
for o = 1:length(ovals), text(2.1, ovals(o), osids{o}, 'Interpreter','none'), end

%% Fig 2: Parameter fits

figure(2); clf; set(gcf, 'Name', ['Parameter fits (' sGroup ': ' sCondition ')']);
hold on; grid on;

%%%%% minkowski
subplot(1,2,1); hold on; grid on;
quickPlot([all_vep.n], [all_psy.n], grp_vep, grp_psy, 0);
ylabel('minkowski parameter'); title('minkowski model')
set(gca, 'FontSize', 12);
set(gca, 'YScale', 'log');
[ovals, osids] = checkoutlier([all_vep.n], {all_vep.sID});
for o = 1:length(ovals), text(1.1, ovals(o), osids{o}, 'Interpreter','none'), end
[ovals, osids] = checkoutlier([all_psy.n], {all_psy.sID});
for o = 1:length(ovals), text(2.1, ovals(o), osids{o}, 'Interpreter','none'), end
set(gca, 'YMinorGrid', 'off');
% print values to command line
disp('MINKOWSKI FITS: ')
tmp = [all_vep.n]; vep_n = tmp(grp_vep==3); clear tmp;
tmp = [all_psy.n]; psy_n = tmp(grp_vep==3); clear tmp;
printValsToCommandLine(vep_n, psy_n, 3);

%%%%% mean/max mixture model
subplot(1,2,2); hold on; grid on;
quickPlot([all_vep.w], [all_psy.w], grp_vep, grp_psy, 0);
ylabel('mean (0) vs max (1) mixture model'); title('mixture model')
set(gca, 'FontSize', 12);
[ovals, osids] = checkoutlier([all_vep.w], {all_vep.sID});
for o = 1:length(ovals), text(1.1, ovals(o), osids{o}, 'Interpreter','none'), end
[ovals, osids] = checkoutlier([all_psy.w], {all_psy.sID});
for o = 1:length(ovals), text(2.1, ovals(o), osids{o}, 'Interpreter','none'), end

disp('MIXTURE MODEL FITS: ')
tmp = [all_vep.w]; vep_w = tmp(grp_vep==3); clear tmp;
tmp = [all_psy.w]; psy_w = tmp(grp_vep==3); clear tmp;
printValsToCommandLine(vep_w, psy_w, 3);


%% weighted mean
figure(3); clf; set(gcf, 'Name', ['Weighted mean fits (' sGroup ': ' sCondition ')']);
hold on; grid on;
subplot(1,2,1); hold on; grid on;
quickPlot([all_vep.wa], [all_psy.wa], grp_vep, grp_psy, 0);
ylabel('w*L + (1-w)*R'); title('weighted average')
set(gca, 'FontSize', 12);
[ovals, osids] = checkoutlier([all_vep.wa], {all_vep.sID});
for o = 1:length(ovals), text(1.1, ovals(o), osids{o}, 'Interpreter','none'), end
[ovals, osids] = checkoutlier([all_psy.wa], {all_psy.sID});
for o = 1:length(ovals), text(2.1, ovals(o), osids{o}, 'Interpreter','none'), end
% print values to command line
disp('WEIGHTED AVG FITS: ')
tmp = [all_vep.wa]; vep_wa = tmp(grp_vep==3); clear tmp;
tmp = [all_psy.wa]; psy_wa = tmp(grp_vep==3); clear tmp;
printValsToCommandLine(vep_wa, psy_wa, 3);





%% Subfunctions for plotting
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

function [ovals, osids] = checkoutlier(values, sIDs)
oidx = isoutlier(values, 'quartiles', 'ThresholdFactor',2);
ovals = values(oidx);
osids = sIDs(oidx);
end

function printValsToCommandLine(vepVals, psyVals, roundTo)

% grab relevant numbers for abstract

disp([' ** VEP data (n = ' num2str(length(vepVals)) '):'])
disp(['    Mean: ' num2str(mean(vepVals),roundTo) ', Median: ' num2str(median(vepVals),roundTo) ...
    ', SD: ' num2str(std(vepVals),roundTo) ', SE: ' num2str(std(vepVals)/length(vepVals),roundTo) ])

disp([' ** Psychophysics data (n = ' num2str(length(psyVals)) '):'])
disp(['    Mean: ' num2str(mean(psyVals),roundTo) ', Median: ' num2str(median(psyVals),roundTo) ...
    ', SD: ' num2str(std(psyVals),roundTo) ', SE: ' num2str(std(psyVals)/length(psyVals),roundTo) ])

[h,p,ci,stats] = ttest(vepVals, psyVals);
disp(' *** Paired t-test:')
disp(['       t(' num2str(stats.df) ') = ' num2str(abs(stats.tstat),roundTo) ', p = ' num2str(p,roundTo)])

[p,h,stats] = signrank(vepVals, psyVals,'method','approximate');
disp(' *** Paired sign-rank test:')
disp(['       Z = ' num2str(abs(stats.zval),roundTo) ', p = ' num2str(p,roundTo)])

log_vep = log(vepVals); log_psy = log(psyVals);
[h,p,ci,stats] = ttest(log_vep, log_psy);
disp(' *** Paired t-test on log-transformed values:')
disp(['       t(' num2str(stats.df) ') = ' num2str(abs(stats.tstat),roundTo) ', p = ' num2str(p,roundTo)])

end