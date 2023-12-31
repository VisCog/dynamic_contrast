%% plot_psychophysics_modelpredictions

% Has 4 plots, one for each (mean/max/mink/mix)
%
% can plot the mean response for an individual compared to these models
% (see flag around line 105-110 to turn off/on)
%
% and/or can plot the overall group mean response (see flag around line
% 235-240 to rurn off/on) 


clear
%close all

%% load
%
%%% datatype:
%  'psychophyics' - motor response, original behavioral paper
%  'vep' - vep response, data from the vep+beh study
%  'vep_psychophysics' - motor response, data from the vep+beh study
datatype = 'vep_psychophysics';
%datatype = 'vep';

%%% condition:
%  'congruent' or 'orthogonal' - only used for vep+beh study behavior
%condition = [];
condition = 'congruent';

switch lower(datatype)
    case 'psychophysics'
        subjectList_psychophysics; % puts variable called sID in workspace
    case {'vep', 'vep_psychophysics'}
        subjectList_vep; % puts variable called sID in workspace
end

%% select subjects
% all_sID = sID;

% look at ns only right now
% select normally-sighted participants only:
idx = cellfun(@(x) strcmpi(x(1:2),'NS'), sID);
all_sID = sID(idx);


%% folders
rawDataDir = [cd filesep 'output_' datatype]; % output folder (original data)
modelFitDir = [cd filesep 'fitdata_' datatype filesep 'model_fits']; % where model fit p is saved
saveResultsDir = [cd filesep 'fitdata_' datatype filesep 'plots']; % where to save

%%

allSubjectResponse = nan(length(all_sID), 720);

for i = 1:length(all_sID) % replace this with the sID # to run only 1 person
    sID = all_sID{i};
    disp(sID)

    switch lower(datatype)
        case 'psychophysics'
            dataFile = dir([rawDataDir filesep sID '_run_*']);
            modelFile = dir([modelFitDir filesep sID '*regular.mat']);
        case {'vep_psychophysics'}
            dataFile = dir([rawDataDir filesep sID  '-motor-', condition '.mat']);
            modelFile = dir([modelFitDir filesep sID '_' condition '.mat']);
        case {'vep'}
            dataFile = dir([rawDataDir filesep sID  '-vep-', condition '.mat']);
            modelFile = dir([modelFitDir filesep sID '_' condition '.mat']);
        otherwise
            disp(['datatype not defined: ' datatype])
            return;
    end

    %% load data and params
    if isempty(dataFile)
        disp(' ~~~ no data for this condition ~~~')
        continue
    end
    data = load([dataFile.folder filesep dataFile.name]);
    if strcmpi(datatype(1:3), 'vep')
        switch lower(datatype)
            case 'vep'
                data = eval(['data.' condition 'Vep']);
            case 'vep_psychophysics'
                data = eval(['data.' condition 'Motor']);
        end
    end

    model = load([modelFile.folder filesep modelFile.name]);
    p = model.p; clear model;

    if strcmp(datatype, 'vep_psychophysics')
        % nan trials on which they were NOT using the joystick:
        data.experiment.LEcontrast(~data.conditionInfo.joyUsedIndex,:) = NaN;
        data.experiment.REcontrast(~data.conditionInfo.joyUsedIndex,:) = NaN;
        data.experiment.response(~data.conditionInfo.joyUsedIndex,:) = NaN;

    end

    allSubjectResponse(i,1:length(p.meanResponseScaled)) = p.meanResponseScaled;

    %% Plot - calculations

    plotYmin = 0 - 0.3;
    plotYmax = 1 + 0.3;

    fastSin = ((sin(2*pi*data.t/6)+1)/2);
    slowSin = ((sin(2*pi*data.t/8)+1)/2);

    t = 0 : p.dt : ((length(p.meanResponseScaled)-1)*p.dt);


    %% Plot - figure
    if 0 % flag for plotting all subjects

        fs=20;%font size
        fn='Arial';
        labelsOn = 1;

        fig1 = figure(1); set(gcf, 'Name', [sID ' (' datatype ')']);
        clf; hold on;
        tiledlayout(4,2);

        % Top panel
        nexttile([1 2]); % mean model
        hold on; set(gca, 'FontSize', fs, 'FontName', fn);

        pC1 = plot(data.t, fastSin, 'Color', [51 76 133]/255, ...   %  fast stim
            'LineStyle',':', 'LineWidth', 3);
        pC2 = plot(data.t, slowSin, 'Color', [175 134 53]/255, ...  %  slow stim
            'LineStyle',':', 'LineWidth', 3);
        pM = plot(t, p.predModel_meanModel,...           % model
            '-', 'LineWidth', 5, 'Color', 'red', 'DisplayName', 'mean model');
        pR = plot(t, p.meanResponseScaled,...           % joystick position (calibrated)
            '-', 'LineWidth', 5, 'Color', 'black', 'DisplayName', 'joystick');

        % axies
        xlim([min(data.t)-.5 round(data.t(end))+.5]);
        xticks([0:6:round(data.t(end))]);
        ylim([plotYmin plotYmax]);
        yticks([0 0.5 1]);
        yticklabels({'0.0', '0.5', '1.0'});
        if labelsOn == 1
            xlabel('time (sec)')
            ylabel('contrast')
        end
        %   title('mean model')
        title([sID ' - ' datatype], 'interpreter', 'none')
        legend([pM pR], 'Location', 'eastoutside');
        text(0, 1.2, ['model MSE: ' num2str(p.meanModelErr,4)])

        % Middle Panel
        nexttile([1 2]); % max model
        hold on; set(gca, 'FontSize', fs, 'FontName', fn);

        pC1 = plot(data.t, fastSin, 'Color', [51 76 133]/255, ...   %  fast stim
            'LineStyle',':', 'LineWidth', 3);
        pC2 = plot(data.t, slowSin, 'Color', [175 134 53]/255, ...  %  slow stim
            'LineStyle',':', 'LineWidth', 3);
        pM = plot(t, p.predModel_maxModel,...           % model
            '-', 'LineWidth', 5, 'Color', 'red', 'DisplayName', 'max model');
        pR = plot(t, p.meanResponseScaled,...           % joystick position (calibrated)
            '-', 'LineWidth', 5, 'Color', 'black', 'DisplayName', 'joystick');
        %   title('max model')

        % axies
        xlim([min(data.t)-.5 round(data.t(end))+.5]);
        xticks([0:6:round(data.t(end))]);
        ylim([plotYmin plotYmax]);
        yticks([0 0.5 1]);
        yticklabels({'0.0', '0.5', '1.0'});
        if labelsOn == 1
            xlabel('time (sec)')
            ylabel('contrast')
        end
        legend([pM pR], 'Location', 'eastoutside');
        text(0, 1.2, ['model MSE: ' num2str(p.maxModelErr,4)])

        % Bottom Panel
        nexttile([1 2]); % minkowski model
        hold on; set(gca, 'FontSize', fs, 'FontName', fn);

        pC1 = plot(data.t, fastSin, 'Color', [51 76 133]/255, ...   %  fast stim
            'LineStyle',':', 'LineWidth', 3);
        pC2 = plot(data.t, slowSin, 'Color', [175 134 53]/255, ...  %  slow stim
            'LineStyle',':', 'LineWidth', 3);
        pM = plot(t, p.predModel_minkModel,...           % model
            '-', 'LineWidth', 5, 'Color', 'red', 'DisplayName', 'minkowski model');
        pR = plot(t, p.meanResponseScaled,...           % joystick position (calibrated)
            '-', 'LineWidth', 5, 'Color', 'black', 'DisplayName', 'joystick');
        %   title('minkowski model')

        % axies
        xlim([min(data.t)-.5 round(data.t(end))+.5]);
        xticks([0:6:round(data.t(end))]);
        ylim([plotYmin plotYmax]);
        yticks([0 0.5 1]);
        yticklabels({'0.0', '0.5', '1.0'});
        if labelsOn == 1
            xlabel('time (sec)')
            ylabel('contrast')
        end
        legend([pM pR], 'Location', 'eastoutside');
        text(0, 1.2, ['model MSE: ' num2str(p.minkModelErr,4) '; p.n = ' num2str(p.n,2)])

        % New, bottom-est panel
        nexttile([1 2]); % min/max weight model
        hold on; set(gca, 'FontSize', fs, 'FontName', fn);

        pC1 = plot(data.t, fastSin, 'Color', [51 76 133]/255, ...   %  fast stim
            'LineStyle',':', 'LineWidth', 3);
        pC2 = plot(data.t, slowSin, 'Color', [175 134 53]/255, ...  %  slow stim
            'LineStyle',':', 'LineWidth', 3);
        pM = plot(t, p.predModel_mnmxwghtModel,...           % model
            '-', 'LineWidth', 5, 'Color', 'red', 'DisplayName', 'mean (0) max (1) mixture model');
        pR = plot(t, p.meanResponseScaled,...           % joystick position (calibrated)
            '-', 'LineWidth', 5, 'Color', 'black', 'DisplayName', 'joystick');
        %   title('minkowski model')

        % axies
        xlim([min(data.t)-.5 round(data.t(end))+.5]);
        xticks([0:6:round(data.t(end))]);
        ylim([plotYmin plotYmax]);
        yticks([0 0.5 1]);
        yticklabels({'0.0', '0.5', '1.0'});
        if labelsOn == 1
            xlabel('time (sec)')
            ylabel('contrast')
        end
        legend([pM pR], 'Location', 'eastoutside');
        text(0, 1.2, ['model MSE: ' num2str(p.minkModelErr,4) '; p.n = ' num2str(p.n,2)])

        %%% fig save
        savname = [saveResultsDir filesep sID '-figplot-trialmeans-' condition '.fig'];
        saveas(fig1,savname)
        input('Press enter for next plot')
    end
end

if 1 % flag for plotting group

    tmp = p;
    clear p
    p.costflag=0;

    p.n = 10; % chose a value to send in here for minkowski
    p.w = 0.8; % chose a value to send in here for min/max weight
    [mink_err, ~, mink_pred] = b_s.minkowski(p, [fastSin;slowSin]', mean(allSubjectResponse)');
    [mix_err, ~, mix_pred] = b_s.meanmax_weighted(p, [fastSin;slowSin]', mean(allSubjectResponse)');

    mean_prediction = mean([fastSin;slowSin]);
    max_prediction = max([fastSin;slowSin]);

    t=data.t;

    fs=20;%font size
    fn='Arial';
    labelsOn = 1;

    fig1 = figure(1); set(gcf, 'Name', ['Mean of all participants (n = ' num2str(length(all_sID)) ', ' datatype ' data)']);
    clf; hold on;
    tiledlayout(4,2);

    % Top panel
    nexttile([1 2]); % mean model
    hold on; set(gca, 'FontSize', fs, 'FontName', fn);

    pC1 = plot(t, fastSin, 'Color', [51 76 133]/255, ...   %  fast stim
        'LineStyle',':', 'LineWidth', 3);
    pC2 = plot(t, slowSin, 'Color', [175 134 53]/255, ...  %  slow stim
        'LineStyle',':', 'LineWidth', 3);
    pM = plot(t, mean_prediction,...           % model
        '-', 'LineWidth', 5, 'Color', 'red', 'DisplayName', 'mean model');
    pR = plot(t, mean(allSubjectResponse),...           % joystick position (calibrated)
        '-', 'LineWidth', 5, 'Color', 'black', 'DisplayName', 'joystick');

    % axies
    xlim([min(t)-.5 round(t(end))+.5]);
    xticks([0:6:round(t(end))]);
    ylim([plotYmin plotYmax]);
    yticks([0 0.5 1]);
    yticklabels({'0.0', '0.5', '1.0'});
    if labelsOn == 1
        xlabel('time (sec)')
        ylabel('contrast')
    end
    title(['Mean of all participants (n = ' num2str(length(all_sID)) ', ' datatype ' data)'], 'interpreter', 'none')
    legend([pM pR], 'Location', 'eastoutside');

    % Middle Panel
    nexttile([1 2]); % max model
    hold on; set(gca, 'FontSize', fs, 'FontName', fn);

    pC1 = plot(t, fastSin, 'Color', [51 76 133]/255, ...   %  fast stim
        'LineStyle',':', 'LineWidth', 3);
    pC2 = plot(t, slowSin, 'Color', [175 134 53]/255, ...  %  slow stim
        'LineStyle',':', 'LineWidth', 3);
    pM = plot(t, max_prediction,...           % model
        '-', 'LineWidth', 5, 'Color', 'red', 'DisplayName', 'max model');
    pR = plot(t, mean(allSubjectResponse),...           % joystick position (calibrated)
        '-', 'LineWidth', 5, 'Color', 'black', 'DisplayName', 'joystick');

    % axies
    xlim([min(t)-.5 round(t(end))+.5]);
    xticks([0:6:round(t(end))]);
    ylim([plotYmin plotYmax]);
    yticks([0 0.5 1]);
    yticklabels({'0.0', '0.5', '1.0'});
    if labelsOn == 1
        xlabel('time (sec)')
        ylabel('contrast')
    end
    legend([pM pR], 'Location', 'eastoutside');

    % Bottom Panel
    nexttile([1 2]); % minkowski model
    hold on; set(gca, 'FontSize', fs, 'FontName', fn);

    pC1 = plot(t, fastSin, 'Color', [51 76 133]/255, ...   %  fast stim
        'LineStyle',':', 'LineWidth', 3);
    pC2 = plot(t, slowSin, 'Color', [175 134 53]/255, ...  %  slow stim
        'LineStyle',':', 'LineWidth', 3);
    pM = plot(t, mink_pred,...           % model
        '-', 'LineWidth', 5, 'Color', 'red', 'DisplayName', 'minkowski model');
    pR = plot(t, mean(allSubjectResponse),...           % joystick position (calibrated)
        '-', 'LineWidth', 5, 'Color', 'black', 'DisplayName', 'joystick');

    % axies
    xlim([min(t)-.5 round(t(end))+.5]);
    xticks([0:6:round(t(end))]);
    ylim([plotYmin plotYmax]);
    yticks([0 0.5 1]);
    yticklabels({'0.0', '0.5', '1.0'});
    if labelsOn == 1
        xlabel('time (sec)')
        ylabel('contrast')
    end
    legend([pM pR], 'Location', 'eastoutside');
    text(0, 1.2, ['minkowski parameter (p.n) = ' num2str(p.n)], 'FontSize', 12)

 nexttile([1 2]); % new even more bottom panel for mean/max mix
    hold on; set(gca, 'FontSize', fs, 'FontName', fn);

    pC1 = plot(t, fastSin, 'Color', [51 76 133]/255, ...   %  fast stim
        'LineStyle',':', 'LineWidth', 3);
    pC2 = plot(t, slowSin, 'Color', [175 134 53]/255, ...  %  slow stim
        'LineStyle',':', 'LineWidth', 3);
    pM = plot(t, mix_pred,...           % model
        '-', 'LineWidth', 5, 'Color', 'red', 'DisplayName', 'mean (0) max (1) mixture model');
    pR = plot(t, mean(allSubjectResponse),...           % joystick position (calibrated)
        '-', 'LineWidth', 5, 'Color', 'black', 'DisplayName', 'joystick');

    % axies
    xlim([min(t)-.5 round(t(end))+.5]);
    xticks([0:6:round(t(end))]);
    ylim([plotYmin plotYmax]);
    yticks([0 0.5 1]);
    yticklabels({'0.0', '0.5', '1.0'});
    if labelsOn == 1
        xlabel('time (sec)')
        ylabel('contrast')
    end
    legend([pM pR], 'Location', 'eastoutside');
    text(0, 1.2, ['Mixed model weight (p.w) = ' num2str(p.w)], 'FontSize', 12)


    


end