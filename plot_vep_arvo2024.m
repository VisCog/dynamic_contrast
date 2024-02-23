clear

sCondition = 'congruent';

%% load files

vepFolder = [cd filesep 'fitdata_vep' filesep 'model_fits'];     % vep
psyFolder = [cd filesep 'fitdata_vep_psychophysics' filesep 'model_fits']; %psychophysics

vepFiles_AM = dir([vepFolder filesep 'AM_*' sCondition '.mat']);  
psyFiles_AM = dir([psyFolder filesep 'AM_*' sCondition '.mat']);

vepFiles_NS = dir([vepFolder filesep 'NS_*' sCondition '.mat']);  
psyFiles_NS = dir([psyFolder filesep 'NS_*' sCondition '.mat']);

% compile
AM_vep = [];
for i = 1:length(vepFiles_AM)
    load([vepFiles_AM(i).folder filesep vepFiles_AM(i).name]);
    AM_vep = [AM_vep p];
    clear p
end
AM_vep = arrayfun(@(x) setfield(x, 'kae', min(x.k_arvo)), AM_vep);

AM_psy = [];
for i = 1:length(psyFiles_AM)
    load([psyFiles_AM(i).folder filesep psyFiles_AM(i).name]);
    AM_psy = [AM_psy p];
    clear p
end
AM_psy = arrayfun(@(x) setfield(x, 'kae', min(x.k_arvo)), AM_psy);

NS_vep = [];
for i = 1:length(vepFiles_NS)
    load([vepFiles_NS(i).folder filesep vepFiles_NS(i).name]);
    NS_vep = [NS_vep p];
    clear p
end
NS_vep = arrayfun(@(x) setfield(x, 'kae', min(x.k_arvo)), NS_vep);

NS_psy = [];
for i = 1:length(psyFiles_NS)
    load([psyFiles_NS(i).folder filesep psyFiles_NS(i).name]);
    NS_psy = [NS_psy p];
    clear p
end
NS_psy = arrayfun(@(x) setfield(x, 'kae', min(x.k_arvo)), NS_psy);


%%


figure(1); clf; hold on;
%am
scatter([1 2],[[AM_vep.kae]' [AM_psy.kae]'], 'r', 'filled') % am 
%scatter([1 2]+0.2, [mean([AM_vep.kae]) mean([AM_psy.kae])], 70, 'k', 'filled')
for i = 1:length(AM_vep)
    plot([1 2], [AM_vep(i).kae AM_psy(i).kae], 'Color', [0.7 0.7 0.7], 'LineStyle', '--'); 
end

%ns
scatter([3 4],[[NS_vep.kae]' [NS_psy.kae]'], 'b', 'filled') % NS 
%scatter([3 4]+0.2, [mean([NS_vep.kae]) mean([NS_psy.kae])], 70, 'k', 'filled')
for i = 1:length(NS_vep)
    plot([3 4], [NS_vep(i).kae NS_psy(i).kae], 'Color', [0.7 0.7 0.7], 'LineStyle', '--'); 
end

xlim([0.5 4.5])
set(gca, 'XTick', [1 2 3 4])
set(gca, 'XTickLabel', {'AM: VEP', 'AM: Behavior', 'NS: VEP', 'NS: Behavior'})
ylim([0 1])
ylabel('k_A_E')
title('linear attenuation estimated using mono trial phase only')
set(gca, 'FontSize', 16)

%%

% model error: meanS, maxS, wS, meanSwgtd, maxSwghtd, wSwghted
tmp_vep_err_NS = [[NS_vep.meanModel_err]' [NS_vep.maxModel_err]' [NS_vep.mnmx_err]'...
    [NS_vep.meanModel_wghtd_err]' [NS_vep.maxModel_wghtd_err]' [NS_vep.mnmx_wght_err]'];
tmp_vep_err_AM = [[AM_vep.meanModel_err]' [AM_vep.maxModel_err]' [AM_vep.mnmx_err]'...
    [AM_vep.meanModel_wghtd_err]' [AM_vep.maxModel_wghtd_err]' [AM_vep.mnmx_wght_err]'];

tmp_psy_err_NS = [[NS_psy.meanModel_err]' [NS_psy.maxModel_err]' [NS_psy.mnmx_err]'...
    [NS_psy.meanModel_wghtd_err]' [NS_psy.maxModel_wghtd_err]' [NS_psy.mnmx_wght_err]'];
tmp_psy_err_AM = [[AM_psy.meanModel_err]' [AM_psy.maxModel_err]' [AM_psy.mnmx_err]'...
    [AM_psy.meanModel_wghtd_err]' [AM_psy.maxModel_wghtd_err]' [AM_psy.mnmx_wght_err]'];


figure(2);clf;hold on;

subplot(2,1,1); hold on;
p1 = scatter([1:size(tmp_vep_err_NS,2)]+0.05, tmp_vep_err_NS, 'b', 'DisplayName','control');
p2 = scatter([1:size(tmp_vep_err_AM,2)]-0.05, tmp_vep_err_AM, 'r', 'DisplayName', 'amblyopia');

p3 = scatter([1:6]+0.15, (mean(tmp_vep_err_NS)), 50, 'b', 'filled', 'DisplayName','control M');
p4 = scatter([1:6]-0.15, (mean(tmp_vep_err_AM)), 50, 'r', 'filled', 'DisplayName','amblyopia M');

meanNS = mean(tmp_vep_err_NS);
meanAM = mean(tmp_vep_err_AM);

for i=1:6
text(i+0.15, meanNS(i), ['  ' num2str(meanNS(i),2)], 'HorizontalAlignment', 'left')
text(i-0.15, meanAM(i), [num2str(meanAM(i),2) '  '], 'HorizontalAlignment', 'right')
end
xticks([1:6])
xticklabels({'Mean (S)', 'Max (S)', 'mixture (S)',...
    'Mean (wghtd S)', 'Max (wghtd S)', 'mixture (wghtd S)'})

%legend([p1(1) p2(1) p3(1) p4(1)], 'Location', 'northeast','box', 'off')

set(gca,'FontSize', 14)
ylabel('MSE')
title('VEP model error')
xlim([0.75 6.25])

subplot(2,1,2); hold on;
p1 = scatter([1:size(tmp_psy_err_NS,2)]+0.05, tmp_psy_err_NS, 'b', 'DisplayName','control');
p2 = scatter([1:size(tmp_psy_err_AM,2)]-0.05, tmp_psy_err_AM, 'r', 'DisplayName', 'amblyopia');

p3 = scatter([1:6]+0.15, (mean(tmp_psy_err_NS)), 50, 'b', 'filled', 'DisplayName','control M');
p4 = scatter([1:6]-0.15, (mean(tmp_psy_err_AM)), 50, 'r', 'filled', 'DisplayName','amblyopia M');

meanNS = mean(tmp_psy_err_NS);
meanAM = mean(tmp_psy_err_AM);

for i=1:6
text(i+0.15, meanNS(i), ['  ' num2str(meanNS(i),2)], 'HorizontalAlignment', 'left')
text(i-0.15, meanAM(i), [num2str(meanAM(i),2) '  '], 'HorizontalAlignment', 'right')
end
xticks([1:6])
xticklabels({'Mean (S)', 'Max (S)', 'mixture (S)',...
    'Mean (wghtd S)', 'Max (wghtd S)', 'mixture (wghtd S)'})

legend([p1(1) p2(1) p3(1) p4(1)], 'Location', 'northeast')
xlim([0.75 6.25])

set(gca,'FontSize', 14)
ylabel('MSE')
title('Psychophysics model error')



%%


figure(3); clf; hold on;
%am
scatter([1 2],[[AM_vep.w_orig]' [AM_psy.w_orig]'], 'r', 'filled') % am 
for i = 1:length(AM_vep)
    plot([1 2], [AM_vep(i).w_orig AM_psy(i).w_orig], 'Color', [0.7 0.7 0.7], 'LineStyle', '--'); 
end

%ns
scatter([3 4],[[NS_vep.w_orig]' [NS_psy.w_orig]'], 'b', 'filled') % NS 
for i = 1:length(NS_vep)
    plot([3 4], [NS_vep(i).w_orig NS_psy(i).w_orig], 'Color', [0.7 0.7 0.7], 'LineStyle', '--'); 
end

%am
scatter([5 6],[[AM_vep.w_wght]' [AM_psy.w_wght]'], 'r', 'filled') % am 
for i = 1:length(AM_vep)
    plot([5 6], [AM_vep(i).w_wght AM_psy(i).w_wght], 'Color', [0.7 0.7 0.7], 'LineStyle', '--'); 
end

%ns
scatter([7 8],[[NS_vep.w_wght]' [NS_psy.w_wght]'], 'b', 'filled') % NS 
for i = 1:length(NS_vep)
    plot([7 8], [NS_vep(i).w_wght NS_psy(i).w_wght], 'Color', [0.7 0.7 0.7], 'LineStyle', '--'); 
end



xlim([0.5 8.5])
set(gca, 'XTick', [1.5 3.5 5.5 7.5])
set(gca, 'XTickLabel', {'AM (S)', 'NS (S)',...
    'AM (weighted S)', 'NS (weighted S)'})
ylabel('w')
title('parameter for mixture model w*mean + (1-w)*max')
set(gca, 'FontSize', 16)
grid on
text(1:8, ones(1,8)*-0.75, {'VEP', 'Psycho', 'VEP', 'Psycho', 'VEP', 'Psycho', 'VEP', 'Psycho'},...
    'HorizontalAlignment','center', 'FontSize', 16)

%%
figure(4); clf; hold on;

%am
scatter([1 2],[[AM_vep.kae]' [AM_psy.kae]'], 'r', 'filled') % am 
for i = 1:length(AM_vep)
    plot([1 2], [AM_vep(i).kae AM_psy(i).kae], 'Color', [0.7 0.7 0.7], 'LineStyle', '--'); 
end

%ns
scatter([3 4],[[NS_vep.kae]' [NS_psy.kae]'], 'b', 'filled') % NS 
for i = 1:length(NS_vep)
    plot([3 4], [NS_vep(i).kae NS_psy(i).kae], 'Color', [0.7 0.7 0.7], 'LineStyle', '--'); 
end


xlim([0.5 4.5])
ylim([-0.05 1.05])
set(gca, 'XTick', [1.5 3.5 5.5 7.5])
set(gca, 'XTickLabel', {'amblyopia', 'control'})
ylabel('k')
set(gca, 'FontSize', 16)
%grid on
text(1:4, ones(1,4)*-0.15, {'VEP', 'behavior', 'VEP', 'behavior'},...
    'HorizontalAlignment','center', 'FontSize', 16)