% plot_psychophysics_individualtrials
%
% plots an example joystick response on a single trial
%

% Plots "fellow" fast on top, "amblyopic" fast on bottom
clear;
%sID = 'NS_G_QP_15';
%sID = 'NS_G_AU_22';
%sID = 'AM_LE_L_WA_43';
%sID = 'AM_RE_L_RK_15';
sID = 'BD_G_AU_47';
%sID = 'BD_LE_G_ES_62';
%sID = 'AM_LE_L_CJ_38';
sID = 'AM_RE_IA_24'; %Kelly vep

datatype = 'vep_psychophysics';
condition = 'congruent';

rawDataDir = [cd filesep 'output_' datatype]; % where the raw data live
saveResultsDir = [cd filesep 'fitdata_' datatype filesep 'model_fits']; % where to save


file = dir([rawDataDir filesep sID '*_run*.mat']);
data = load([file.folder filesep file.name]);

loadThis = dir([saveResultsDir filesep sID '-' condition '.mat']);
load([loadThis.folder filesep loadThis.name]);


% for plotting individual trials:
p.costflag = 0;
[err, predModel, respJoy, respJoyCalib, ...
    t, S, stim] = b_s.getErr(p, data);

% plot means, separated into "faster eye = left" and "faster eye = right"
% data.config.fastEye: 0 = left, 1 = right
LEfastTrials = find(~logical(data.config.fastEye)); %#ok<*UNRCH>
REfastTrials = find(logical(data.config.fastEye));

% convert cells to mats
respJoyMAT = cell2mat(respJoy');
predModelMAT = cell2mat(predModel');

% divide into LE/RE faster
respLEfast = respJoyMAT(LEfastTrials,:);
respREfast = respJoyMAT(REfastTrials,:);

modelLEfast = predModelMAT(LEfastTrials,:);
modelREfast = predModelMAT(REfastTrials,:);

LEfastStim = S(LEfastTrials);
REfastStim = S(REfastTrials);

% for plotting..
%plotYmin = floor(min([min(respJoyCalibMAT(:)) min(respJoyCalibMAT(:))])*10)/10 - 0.1;
%plotYmax = ceil(max([max(respJoyCalibMAT(:)) max(respJoyCalibMAT(:))])*10)/10 + 0.1;
plotYmin = 0 - 0.15;
plotYmax = 1 + 0.15;

plotTrial = 3;


%% Figure

fs=20;%font size
fn='Arial';

fig1 = figure(1); set(gcf, 'Name', 'Example Individual Trials');
clf; hold on;
tiledlayout(2,2);


nexttile([1 2]); % a fellow eye "fast" trial
hold on; set(gca, 'FontSize', fs, 'FontName', fn);

if p.k(1) == 1 % left eye is fellow

    thisTrialStim = LEfastStim{plotTrial};

    pC1 = plot(data.t, thisTrialStim(:,1), 'Color', [51 76 133]/255, ...
        'LineStyle',':', 'LineWidth', 3);
    pC2 = plot(data.t, thisTrialStim(:,2), 'Color', [175 134 53]/255, ...
        'LineStyle',':', 'LineWidth', 3);

    p4 = plot(data.t, respLEfast(plotTrial,:),... % raw joystick position
        '-', 'LineWidth', 5, 'Color', 'black');
    % modeled result
    %p5 = plot(data.t+p.delay, (modelLEfast(plotTrial,:)-p.intercept)/p.slope,...
    %p5 = plot(data.t+p.delay, (modelLEfast(plotTrial,:))/p.slope,...
    %     '-', 'LineWidth', 5, 'Color', 'red');
    p5 = plot(data.t, (modelLEfast(plotTrial,:)-p.intercept),...
        '-', 'LineWidth', 5, 'Color', 'red');

else % else right eye is fellow

    thisTrialStim = REfastStim{plotTrial};

    pC1 = plot(data.t, thisTrialStim(:,1), 'Color', [51 76 133]/255, ...
        'LineStyle',':', 'LineWidth', 3);
    pC2 = plot(data.t, thisTrialStim(:,2), 'Color', [175 134 53]/255, ...
        'LineStyle',':', 'LineWidth', 3);

    p4 = plot(data.t, respREfast(plotTrial,:),... %
        '-', 'LineWidth', 5, 'Color', 'black', 'DisplayName', 'true slider position');
    % modeled result
    %p5 = plot(data.t+p.delay, (modelREfast(plotTrial,:)-p.intercept)/p.slope,...
    % p5 = plot(data.t+p.delay, (modelREfast(plotTrial,:))/p.slope,...
    %    '-', 'LineWidth', 5, 'Color', 'red', 'DisplayName', 'predicted slider position');
    p5 = plot(data.t, (modelREfast(plotTrial,:)-p.intercept),...
        '-', 'LineWidth', 5, 'Color', 'red', 'DisplayName', 'predicted slider position');
end

xlim([min(data.t) max(data.t)+.5])
xticks([0:6:48])

%xlabel('time (sec)')

%ylabel('contrast')
ylim([plotYmin plotYmax])
% janky way to rescale the plot instead of rescaling the response
%yticks(linspace(plotYmin+0.1, plotYmax-0.1,3))
yticks([0 0.5 1])
yticklabels({'0.0', '0.5', '1.0'})

nexttile([1 2]); % amblyopic eye "fast"
hold on; set(gca, 'FontSize', fs, 'FontName', fn);

if p.k(1) == 1 % left eye is fellow


    thisTrialStim = REfastStim{plotTrial};

    pC1 = plot(data.t, thisTrialStim(:,1), 'Color', [51 76 133]/255, ...
        'LineStyle',':', 'LineWidth', 3);
    pC2 = plot(data.t, thisTrialStim(:,2), 'Color', [175 134 53]/255, ...
        'LineStyle',':', 'LineWidth', 3);

    p4 = plot(data.t, respREfast(plotTrial,:),... % calibrated position
        '-', 'LineWidth', 5, 'Color', 'black', 'DisplayName', 'true slider position');
    % modeled result
    %p5 = plot(data.t+p.delay, ((modelREfast(plotTrial,:))-p.intercept)/p.slope,...
    %p5 = plot(data.t+p.delay, ((modelREfast(plotTrial,:)))/p.slope,...
    %   '-', 'LineWidth', 5, 'Color', 'red', 'DisplayName', 'predicted slider position');
    p5 = plot(data.t, (modelREfast(plotTrial,:)-p.intercept),...
        '-', 'LineWidth', 5, 'Color', 'red', 'DisplayName', 'predicted slider position');
else


    thisTrialStim = LEfastStim{plotTrial};

    pC1 = plot(data.t, thisTrialStim(:,1), 'Color', [51 76 133]/255, ...
        'LineStyle',':', 'LineWidth', 3);
    pC2 = plot(data.t, thisTrialStim(:,2), 'Color', [175 134 53]/255, ...
        'LineStyle',':', 'LineWidth', 3);

    p4 = plot(data.t, respLEfast(plotTrial,:),... % calibrated position
        '-', 'LineWidth', 5, 'Color', 'black');
    % modeled result
    %p5 = plot(data.t+p.delay, ((modelLEfast(plotTrial,:))-p.intercept)/p.slope,...
    %p5 = plot(data.t+p.delay, ((modelLEfast(plotTrial,:)))/p.slope,...
    %    '-', 'LineWidth', 5, 'Color', 'red');
    p5 = plot(data.t, (modelLEfast(plotTrial,:)-p.intercept),...
        '-', 'LineWidth', 5, 'Color', 'red');

end


xlim([min(data.t) max(data.t)+.5])
xticks([0:6:48])

ylim([plotYmin plotYmax])

%xlabel('time (sec)')

%ylabel('contrast')
ylim([plotYmin plotYmax])

yticks([0 0.5 1])
yticklabels({'0.0', '0.5', '1.0'})

set(gcf, 'color', 'w')
set(gcf,'Position',[0 0 13 5])

set(fig1,'Units','Inches');
pos = get(fig1,'Position');
set(fig1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

savname = [saveResultsDir filesep sID '-figplot-exampletrial.pdf'];
print(fig1,savname,'-dpdf','-r0')
