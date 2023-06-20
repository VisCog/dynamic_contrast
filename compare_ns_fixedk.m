% Comparing VEP vs psychophysics fits, NS subjects only, when k1=k2=1 and
% only the normalization values are fit (calibration is also done)
clear
close all

vepFolder = [cd filesep 'fitdata_vep' filesep 'model_fits_fixedk'];
psyFolder = [cd filesep 'fitdata_vep_psychophysics' filesep 'model_fits_fixedk'];


vepFiles = dir([vepFolder filesep '*congruent.mat']);
psyFiles = dir([psyFolder filesep '*congruent.mat']);

% quick fix for that 1 i can't dl right now
psyFiles = psyFiles(2:end);

% quick alt to gathertable
all_psy = [];
for i = 1:length(psyFiles)
    load([psyFiles(i).folder filesep psyFiles(i).name]);
    p.U2 = p.U(2);
    p.U3 = p.U(3);
    all_psy = [all_psy p];
end
all_vep = [];
for i = 1:length(vepFiles)
    load([vepFiles(i).folder filesep vepFiles(i).name]);
    p.U2 = p.U(2);
    p.U3 = p.U(3);
    all_vep = [all_vep p];
end
% catch if things have gone wrong
if ~all(strcmpi({all_vep.sID},{all_psy.sID}))
    disp('!!!!!!!!')
    disp('SUBJECT IDs DO NOT LINE UP')
    return
end

compareFields = {'U2', 'U3', 'offset', 'step2normalizationErr'};%, 'sigma','smax'};

matVEP = [];
matPSY = [];
for i=1:length(compareFields)
    tmpVEP = [all_vep.(compareFields{i})];
    tmpPSY = [all_psy.(compareFields{i})];
    [rv, pv] = corr(tmpVEP', tmpPSY');
    matVEP = [matVEP tmpVEP'];
    matPSY = [matPSY tmpPSY'];

    axVals = [0 ceil(max([tmpVEP tmpPSY]))];

    figure(i); clf;
    line(axVals + [-1 1], axVals + [-1 1], 'Color', [.5 .5 .5], 'LineStyle', '--');
    hold on; grid on;
    scatter(tmpVEP, tmpPSY, 'k', 'filled');
    xlim(axVals); ylim(axVals); axis square;
    refline(polyfit(tmpVEP, tmpPSY, 1));
    title([compareFields{i} ' comparison: r = ' num2str(rv,2) ', p = ' num2str(pv,2)])

    xlabel(['VEP data: ' compareFields{i}]);
    ylabel(['psychophys data: ' compareFields{i}]);
    set(gca, 'FontSize', 13)
end

% plot all means
mnVEP = mean(matVEP);
mnPSY = mean(matPSY);

figure(length(compareFields)+1);
clf;
subplot(1,2,1); hold on; grid on;  set(gca, 'FontSize', 16);
for i=1:length(compareFields)
    swarmchart(repmat(i,1,size(matVEP,1)), matVEP(:,i));
    scatter(i, mnVEP(i), 50, 'k', 'filled');
end
ylim([0 ceil(max([matVEP(:); matPSY(:)]))])
xlim([0 length(compareFields)+1])
xticks(1:length(compareFields))
xticklabels(compareFields)
title('VEP data')

subplot(1,2,2); hold on; grid on; set(gca, 'FontSize', 16);
for i=1:length(compareFields)
    swarmchart(repmat(i,1,size(matPSY,1)), matPSY(:,i));
    s = scatter(i, mnPSY(i), 50, 'k', 'filled');
end
ylim([0 ceil(max([matPSY(:); matPSY(:)]))])
xlim([0 length(compareFields)+1])
xticks(1:length(compareFields))
xticklabels(compareFields)
title('psychophys data')

legend(s, 'mean')
