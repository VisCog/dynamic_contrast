% Comparing VEP vs psychophysics fits, NS subjects only, on a set of
% alternative models

clear
close all

vepFolder = [cd filesep 'fitdata_vep' filesep 'model_fits_altmodels'];
psyFolder = [cd filesep 'fitdata_vep_psychophysics' filesep 'model_fits_altmodels'];


vepFiles = dir([vepFolder filesep '*congruent.mat']);
psyFiles = dir([psyFolder filesep '*congruent.mat']);

% select normally-sighted participants only:
idx = cellfun(@(x) strcmpi(x(1:2),'NS'), {vepFiles.name});
vepFiles = vepFiles(idx);
idx = cellfun(@(x) strcmpi(x(1:2),'NS'), {psyFiles.name});
psyFiles = psyFiles(idx);

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

% grab toy data
exFiles = dir([vepFolder filesep '*xample_congruent.mat']);
toy_data = [];
for i = 1:2
    load([exFiles(i).folder filesep exFiles(i).name]);
    toy_data = [toy_data p];
end
toynames = strrep({exFiles.name}, '_congruent.mat', '');

%%% calibrations
VEPcalibErrMean = [all_vep.errMean];
VEPcalibErrIndv = [all_vep.errInd];

PSYcalibErrMean = [all_psy.errMean];
PSYcalibErrIndv = [all_psy.errInd];

figure(1); clf; set(gcf, 'Name', 'Calibration error')
%mean
subplot(1,2,1); 
quickPlot(VEPcalibErrMean, PSYcalibErrMean, [toy_data.errMean], toynames, 0)
ylabel('Calibration err (mean)')
set(gca, 'FontSize', 12)
%indiv
subplot(1,2,2); 
quickPlot(VEPcalibErrIndv, PSYcalibErrIndv, [toy_data.errInd], toynames, 1)
ylabel('Calibration err (indiv. trials)')
set(gca, 'FontSize', 12)

figure(2); clf; set(gcf, 'Name', 'Model error')
subplot(1,2,1);
quickPlot([all_vep.simpleAverageErr], [all_psy.simpleAverageErr], [toy_data.simpleAverageErr], toynames, 0)
ylabel('Error: Simple Average Model')
set(gca, 'FontSize', 12)
subplot(1,2,2);
quickPlot([all_vep.simpleMaxErr], [all_psy.simpleMaxErr], [toy_data.simpleMaxErr], toynames, 1)
ylabel('Error: Simple Max Model')
set(gca, 'FontSize', 12)

figure(3), clf; set(gcf,'Name', 'Simple softmax')
subplot(1,2,1);
quickPlot([all_vep.simpleSmaxErr], [all_psy.simpleSmaxErr], [toy_data.simpleSmaxErr], toynames, 0)
ylabel('Error: Simple Smax')

set(gca, 'FontSize', 12)
subplot(1,2,2);
quickPlot([all_vep.smax], [all_psy.smax], [toy_data.smax], toynames, 1)
ylabel('smax free param')
set(gca, 'FontSize', 12)


function quickPlot(vepvector, psyvector,toyvector, toynames, legendtf)
scatter(1,vepvector, 'r', 'jitter', 0.2);
hold on; grid on;
scatter(2,psyvector, 'b', 'jitter', 0.2);
xlim([0.5 3.5])
scatter([1 2], [mean(vepvector) mean(psyvector)], 50, 'k', 'filled');
t1 = scatter(2.95, toyvector(1), 'DisplayName', toynames{1});
t2 = scatter(3.05, toyvector(2), 'DisplayName', toynames{2});
xticks([1 2 3])
xticklabels({'VEP', 'PSYCHO', 'TOYDATA'})
if legendtf == 1
    legend([t1 t2]);
end
end