%% plot_psychophysics_meantrial
%
% plots the mean response for an individual over all trials
% Plots "fellow" fast on top, "amblyopic" fast on bottom
%
% Fall 2023: Consider this archival for now - not currently used for
% anything 


clear

%% datatype:
%  'vep' - vep response, data from the vep+beh study
%  'vep_psychophysics' - motor response, data from the vep+beh study
%datatype = 'vep_psychophysics';
datatype = 'vep';
subjectList_vep; % puts variable called sID in workspace
all_sID = sID; clear sID;
%datatype = 'psychophysics';
% condition:
%  'congruent' or 'orthogonal' - only used for vep+beh study behavior
%condition = [];
condition = 'congruent';


%% folders

rawDataDir = [cd filesep 'output_' datatype]; % output folder (original data)
modelFitDir = [cd filesep 'fitdata_' datatype filesep 'model_fits']; % where model fit p is saved
saveResultsDir = [cd filesep 'fitdata_' datatype filesep 'plots']; % where to save
for i = 1%:length(all_sID) % replace this with the sID # to run only 1 person
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
    [err, predModel, respJoy, respJoyCalib, ...
        t, S, stim] = b_s.getErr(p, data);

    if strcmp(datatype, 'vep_psychophysics')
        % nan trials on which they were NOT using the joystick:
        data.experiment.LEcontrast(~data.conditionInfo.joyUsedIndex,:) = NaN;
        data.experiment.REcontrast(~data.conditionInfo.joyUsedIndex,:) = NaN;
        data.experiment.response(~data.conditionInfo.joyUsedIndex,:) = NaN;
    end


    %% Plot - calculations

    % %     % plot means, separated into "faster eye = left" and "faster eye = right"
    % %     % data.config.fastEye: 0 = left, 1 = right
    % %     LEfastTrials = find(~logical(data.config.fastEye)); %#ok<*UNRCH>
    % %     REfastTrials = find(logical(data.config.fastEye));

    % Create indices to select trial portions where each eye is presenting
    % zero contrast information (the "monocular" trial portions)
    % use diff() function to avoid grabbing zeros that are not part of the
    % dropped cycle (i.e. that hit zero for 1 timepoint during regular cycles)
    diffLE = diff(data.experiment.LEcontrast,1,2);
    diffRE = diff(data.experiment.REcontrast,1,2);
    zeroContrastInLE = [zeros(size(diffLE,1),1) (diffLE == 0)];
    zeroContrastInRE = [zeros(size(diffRE,1),1) (diffRE == 0)];

    % index for all monocular trials, left OR right
    monoTrialsIndex = zeroContrastInLE | zeroContrastInRE;

    % convert cells to mats
    respJoyMAT = cell2mat(respJoy');
    respJoyCalibMAT = cell2mat(respJoyCalib');
    predModelMAT = cell2mat(predModel');


    % NaN out all the monocular portions (do not include in averages)
    respJoyMAT(monoTrialsIndex) = NaN;
    respJoyCalibMAT(monoTrialsIndex) = NaN;
    predModelMAT(monoTrialsIndex) = NaN;

    % %     % separate into LE/RE faster
    % %     respLEfast = respJoyMAT(LEfastTrials,:);
    % %     respREfast = respJoyMAT(REfastTrials,:);
    % %
    % %     calibLEfast = respJoyCalibMAT(LEfastTrials,:);
    % %     calibREfast = respJoyCalibMAT(REfastTrials,:);
    % %
    % %     modelLEfast = predModelMAT(LEfastTrials,:);
    % %     modelREfast = predModelMAT(REfastTrials,:);
% 
%     % if change
%     mn_respREfast = nanmean(calibREfast);
%     mn_modelREfast = nanmean(modelREfast);
%     mn_modelREfast = mn_modelREfast-(nanmean(mn_modelREfast)  - nanmean(respREfast));
%     mn_modelREfast = scaleif(mn_modelREfast, 0, 1);
%     mn_respREfast = scaleif(mn_respREfast, 0, 1);
% 
%     mn_respLEfast = nanmean(calibLEfast);
%     mn_modelLEfast = nanmean(modelLEfast);
%     mn_modelLEfast = mn_modelLEfast- (nanmean(mn_modelLEfast)  - nanmean(respLEfast));
%     mn_modelLEfast = scaleif(mn_modelLEfast, 0, 1);
%     mn_respLEfast = scaleif(mn_respLEfast, 0, 1);

    fastSin = ((sin(2*pi*data.t/6)+1)/2);
    slowSin = ((sin(2*pi*data.t/8)+1)/2);
    mn_maxmodel =max([fastSin;slowSin]);
    mn_meanmodel =mean([fastSin;slowSin]);
    p.n = 100;
    mn_resp = nanmean(respJoyCalibMAT, 1)';
    mn_resp = scaleif(mn_resp, 0, 1);
    err = b_s.minkowski(p,[fastSin;slowSin]', mn_resp);

    p = fit('b_s.minkowski', p,{'n'}, [fastSin;slowSin]',mn_resp);
    [err_mink, p, mink_model] = b_s.minkowski(p,[fastSin;slowSin]', mn_resp);

    % get errors
    mn_err = sum((mn_meanmodel-mn_resp).^2);
    max_err = sum((mn_maxmodel-mn_resp).^2);

    plotYmin = 0 - 0.3;
    plotYmax = 1 + 0.3;

    %% Plot - figure

    fs=20;%font size
    fn='Arial';
    labelsOn = 1;

    fig1 = figure(1); set(gcf, 'Name', 'Trial means');
    clf; hold on;
    tiledlayout(2,2);

    % Top panel
    nexttile([1 2]); % fellow eye "fast"
    hold on; set(gca, 'FontSize', fs, 'FontName', fn);


    pC1 = plot(data.t, fastSin, 'Color', [51 76 133]/255, ...   % fellow fast
        'LineStyle',':', 'LineWidth', 3);
    pC2 = plot(data.t, slowSin, 'Color', [175 134 53]/255, ...  % amb slow
        'LineStyle',':', 'LineWidth', 3);
    p4 = plot(data.t-p.delay, mn_respLEfast,...           % joystick position
        '-', 'LineWidth', 5, 'Color', 'black');
    p5 = plot(data.t, mn_modelLEfast,...    % modeled result
        '-', 'LineWidth', 5, 'Color', 'red');
    p_mn =  plot(data.t, mn_meanmodel,...    % modeled result
        '-', 'LineWidth', 3, 'Color',[ 0 .5 0]);
    p_max =  plot(data.t, mn_maxmodel,...    % modeled result
        '--', 'LineWidth', 3, 'Color',[ 0 .5 0]);
    p_max =  plot(data.t, mink_model,...    % modeled result
        '--', 'LineWidth', 3, 'Color',[ 1 .5 0]);


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

    % Bottom Panel
    nexttile([1 2]); % amblyopic eye "fast"
    hold on; set(gca, 'FontSize', fs, 'FontName', fn);
    pC1 = plot(data.t, slowSin, 'Color', [51 76 133]/255, ...   % fellow slow
        'LineStyle',':', 'LineWidth', 3);
    pC2 = plot(data.t, fastSin, 'Color', [175 134 53]/255, ...  % amb fast
        'LineStyle',':', 'LineWidth', 3);
    p4 = plot(data.t-p.delay, mn_respREfast,...           % joystick position
        '-', 'LineWidth', 5, 'Color', 'black');
    p5 = plot(data.t, mn_modelREfast,...    % modeled result
        '-', 'LineWidth', 3, 'Color', 'red');
    p_mn =  plot(data.t, mn_meanmodel,...    % modeled result
        '-', 'LineWidth', 3, 'Color',  [ 0 .5 0]);
    p_max =  plot(data.t, mn_maxmodel,...    % modeled result
        '--', 'LineWidth', 3, 'Color', [ 0 .5 0]);
    text(20, 1 , ['mn err = ', num2str(round(mn_err, 2))]);
    text(20, .8 , ['max err = ', num2str(round(max_err, 2))]);
    text(20, .6 , ['mink val = ', num2str(round(abs(p.n), 3))]);


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

    % Save
    if strcmpi(datatype, 'psychophysics')
        xIn = 11;
    else
        xIn = 7;
    end

    set(fig1,'Units','Inches');
    set(gcf, 'color', 'w');
    set(gcf,'Position',[0 0 xIn 5]);
    pos = get(fig1,'Position');
    set(fig1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

    %%% pdf save
    % savname = [saveResultsDir filesep sID '-figplot-trialmeans.pdf'];
    %print(fig1,savname,'-dpdf','-r0')
    %%% fig save
    savname = [saveResultsDir filesep sID '-figplot-trialmeans.fig'];
    saveas(fig1,savname)
    pause
end