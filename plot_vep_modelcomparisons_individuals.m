%% Plots model fits
% Makes 4 figures, one for each model (mean, max, minkowski, mixture)
% Each figure shows plots for all participant fits

datatype = 'vep';
%datatype = 'vep_psychophysics';

%% Load relevant data
all_data = [];

switch lower(datatype)

    case 'vep'
        vepFolder = [cd filesep 'fitdata_vep' filesep 'model_fits'];     % vep
        vepFiles = dir([vepFolder filesep 'NS_*congruent.mat']);  % grab NS only right now
        for i = 1:length(vepFiles)
            load([vepFiles(i).folder filesep vepFiles(i).name]);
            all_data = [all_data p];
            clear p
        end

    case 'vep_psychophysics'
        psyFolder = [cd filesep 'fitdata_vep_psychophysics' filesep 'model_fits']; %psychophysics
        psyFiles = dir([psyFolder filesep 'NS_*congruent.mat']);
        for i = 1:length(psyFiles)
            load([psyFiles(i).folder filesep psyFiles(i).name]);
            all_data = [all_data p];
            clear p
        end

    otherwise
        disp([datatype ' is an undefined datatype'])
end


%% plot setting
n = size(all_data,2);
numRows = ceil(n / 4);

%% Mean model
fig4 = figure(4); clf; hold on;
set(gcf, "Name", [datatype ' data: Mean model fits'])
for i = 1:n
    subplot(numRows, 4, i); hold on;
    oneSubjectPlot(all_data(i).predModel_meanModel, all_data(i).meanResponseScaled, all_data(i).sID)
    text(200, 0.1, ['err = ' num2str(all_data(i).meanModelErr,2)], 'color', 'b')
    ylabel('contrast');
end
text(800,.8, 'Mean prediction in red (L+R)/2', 'color', 'r')
text(800,.6, 'response in black', 'color', 'k')

%% Max model
fig5 = figure(5); clf; hold on;
set(gcf, "Name", [datatype ' data: Max model fits'])
for i = 1:n
    subplot(numRows, 4, i); hold on;
    oneSubjectPlot(all_data(i).predModel_maxModel, all_data(i).meanResponseScaled, all_data(i).sID)
    text(200, 0.1, ['err = ' num2str(all_data(i).maxModelErr,2)], 'color', 'b')
    ylabel('contrast');
end
text(800,.8, 'Max prediction in red max(L,R)', 'color', 'r')
text(800,.6, 'response in black', 'color', 'k')


%% Minkwoski
fig6 = figure(6); clf; hold on;
set(gcf, "Name", [datatype ' data: Minkowski fits'])
for i = 1:n
    subplot(numRows, 4, i); hold on;
    oneSubjectPlot(all_data(i).predModel_minkModel, all_data(i).meanResponseScaled, all_data(i).sID)
    text(200, 0.2, ['m = ' num2str(all_data(i).n,2)], 'color', 'r')
    text(200, 0.1, ['err = ' num2str(all_data(i).minkModelErr,2)], 'color', 'b')
    ylabel('contrast');
end
text(800,.8, 'Minkowski prediction in red ([L^m+R^m]/2)^(^1^/^m^)', 'color', 'r')
text(800,.6, 'response in black', 'color', 'k')

%% Mean/max weighting
fig7 = figure(7); clf; hold on;
set(gcf, "Name", [datatype ' data: mean/max mixture model fits'])
for i = 1:n
    subplot(numRows, 4, i); hold on;
    oneSubjectPlot(all_data(i).predModel_mnmxwghtModel, all_data(i).meanResponseScaled, all_data(i).sID)
    text(200, 0.2, ['w = ' num2str(all_data(i).w,2)], 'color', 'r')
    text(200, 0.1, ['err = ' num2str(all_data(i).mnmxwghtModelErr,2)], 'color', 'b')

    ylabel('contrast');
end
text(800,.8, 'Mean/max weighting in red (1-w)*mean + (w)*max', 'color', 'r')
text(800,.6, 'response in black', 'color', 'k')


%%

function oneSubjectPlot(predModelData, meanResponseData, sID)
plot(predModelData, 'r');
plot(meanResponseData, 'k');
title(sID, 'Interpreter', 'none');
xlim([1 length(meanResponseData)]); ylim([0 1]);
end