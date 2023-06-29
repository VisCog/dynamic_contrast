%% fit_vep_dynamic_contrast
%
% This script fits  the VEP data (all trials, not separated by
% whether the joystick was used or not) or the psychophysics data (only
% trials where the joystick was used)
% Participant IDs are listed in current directory's file called subjectList_vep.m
%
% <old words deleted> This header needs to be updated


%%


clear
close all
datatype = 'vep'; %'vep' or 'vep_psychophysics'
subjectList_vep; % puts variable called sID in workspace

% June 2023: This script is set up to fit normally-sighted participants
% with equal attenuation (k1 = k2 = 1) if set analysistype = 'ns'
analysistype = 'ns';

if strcmpi(analysistype, 'ns')
    % select normally-sighted participants only:
    idx = cellfun(@(x) strcmpi(x(1:2),'NS'), sID);
    sID = sID(idx);
end

if strcmp(datatype, 'vep')
    rawDataDir = [cd filesep 'output_vep']; % where the raw data live
    saveResultsDir = [cd filesep 'fitdata_vep'  filesep 'model_fits'];



elseif strcmp(datatype, 'vep_psychophysics')
    rawDataDir = [cd filesep 'output_vep_psychophysics']; % where the raw data live
    saveResultsDir = [cd filesep 'fitdata_vep_psychophysics'  filesep 'model_fits'];
end

savePlotOn = 0; % if 1, saves plots
pauseForPlots = 0; % if 1, waits for you to press enter after each plot

if strcmpi(analysistype, 'ns')
    disp('%%%%%%%%%%%%%% You''re running the NS-only, fix k1=k2=1 analysis')
    saveResultsDir = strrep(saveResultsDir, 'model_fits', 'model_fits_fixedk');
    disp('%%%%%%%%%%%%%% If that''s unexpected, change analysistype in the code')
end

%% Individual subjects

for i = 1:length(sID) % replace this with the sID # to run only 1 person
    % we clear p at the end of each individual so these need to be inside
    % the loop
    p.sID = sID{i};
    p.abs = 1;
    p.clean_range  = 0.3; % only calibrate or fit data where there's this must range in the data, p.clean_range = 0, uses all data
    p.condition = 'congruent';

    %%% quick fix for not having this dl'd yet
    if all(sID{i} == 'NS_AA_17'), continue; end

    disp('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    disp(['Working on ' sID{i} ' (' num2str(i) ' of ' num2str(length(sID)) ')'])
    disp([' fitting  condition: ', p.condition]);

    % load raw data:
    file = dir([rawDataDir filesep sID{i} '*-',p.condition, '.mat']);
    load([file.folder filesep file.name]);
    % puts a struct called congruentVep or orthogonalVep into the workplace
    if strcmp(datatype, 'vep')
        data = eval([p.condition, 'Vep']);
    else
        data =  eval([p.condition, 'Motor']);
    end

    figure(10); clf;
    set(gcf, 'Name', [sID{i} '     datatype: ' datatype]);
    imagesc(data.experiment.response); colormap(gray);
    yticks(1:size(data.experiment.response,1));
    if savePlotOn == 1
        % save plot
        saveas(gcf, [saveResultsDir filesep p.sID '-RasterPlot.fig']);
    end
    if pauseForPlots == 1
        input('updated raster plot - press enter to continue')
    end

    %%%%% select only the trials on which they were using the joystick %%%%%
    % Important difference from before: only half of trials contain a
    % response (the other half they were told only to watch the stimulus).
    % the variable data.conditionInfo.joyUsedIndex indicates on which
    % trials participants used the joystick (regardless of whether they
    % were told to or not), use that to select trials for this analysis.
    if strcmp(datatype, 'vep_psychophysics')
        data.experiment.BEcontrastStart = data.experiment.BEcontrastStart(data.conditionInfo.joyUsedIndex,:);
        data.experiment.binoResponseStart = data.experiment.binoResponseStart(data.conditionInfo.joyUsedIndex,:);
        data.experiment.BEcontrastEnd = data.experiment.BEcontrastEnd(data.conditionInfo.joyUsedIndex,:);
        data.experiment.binoResponseEnd = data.experiment.binoResponseEnd(data.conditionInfo.joyUsedIndex,:);
        data.experiment.LEcontrast = data.experiment.LEcontrast(data.conditionInfo.joyUsedIndex,:);
        data.experiment.REcontrast = data.experiment.REcontrast(data.conditionInfo.joyUsedIndex,:);
        data.experiment.response = data.experiment.response(data.conditionInfo.joyUsedIndex,:);
    end

    % A difference between the EEG and previously-collected bino data:
    % in EEG the binocularly-presented contrast is split between the beginning
    % and end of each trial (i.e. 7 seconds at the beginning and 7 seconds at
    % the end).
    % data.experiment.binoResponseStart and .BEcontrastStart hold the data for
    % the start of the trial, and binoResponseEnd / BEcontrastend hold data for
    % the end.

    % b_s expects the data to be inside data.experiment.binoResponse:
    concatContrast  = [data.experiment.BEcontrastStart data.experiment.BEcontrastEnd]; % ntrials x duration
    concatResponse  = [data.experiment.binoResponseStart data.experiment.binoResponseEnd];

    % put everything in  -0.5 to 0.5 units
    data.experiment.binoS= concatContrast(1, :)-0.5;
    data.experiment.binoResponse = concatResponse - 0.5;
    nan_indx = round(length(data.experiment.binoS)/2-3:length(data.experiment.binoS)/2+3);
    data.experiment.binoResponse(:, nan_indx) = NaN; % block out switchy bit when concatenating trial

    [data.experiment.binoResponse, p.n_good]  = b_s.cleanData(data.experiment.binoResponse, p);

    disp(['# Calib trials = ', num2str(size( data.experiment.binoResponse, 1)),  ...
        ' Good Calib trials = ', num2str(p.n_good)]);

    p.startT = 1;

    % default to slope of 1, intercept is mean (usually around 0)
    p.slope = 1;
    p.intercept = -2*nanmean(data.experiment.binoResponse(:));

    p.junk = 0;
    binoMean.gvals = [];
    p.dt = diff(data.binocular.t(1:2));
    data.binocular.t = 0:p.dt:p.dt*size(data.experiment.binoResponse, 2)-p.dt;

    if strcmp(datatype, 'vep')
        freeList = {'slope', 'intercept', 'delay'};
        p. joystickfunction = 'delay + scale';
        p.delay = 0.15; p.penalizeDelay = 3; % delay  penalization in seconds
    elseif strcmp(datatype, 'vep_psychophysics')
        freeList = {'slope', 'intercept', 'delay'};
        p.delay =2; p.penalizeDelay = 3; % delay  penalization in seconds, based on previous psychophysical data
        p. joystickfunction = 'delay + scale';
    end
    p.costflag = 1; p = fit('b_s.getErrBinoMean', p, freeList, data);

    % Calculate the error, both for mean joystick position and for
    % indivdual trials
    p.costflag = 0;
    [~, p.errMean, data] = b_s.getErrBinoMean(p,data);
    [~, p.errInd, ~] = b_s.getErrBinoInd(p,data);

    % display fits and error
    disp(['   .. delay: ' num2str(round(p.delay,4)) '    intercept: ' num2str(round(p.intercept,4)) '    slope: ' num2str(round(p.slope,3))])
    disp(['   .. MSE:   for mean response: ' num2str(round(p.errMean,4)) '    for individual trials: ' num2str(round(p.errInd,4)) ]);
    figure(1); clf;  set(gcf, 'Name', 'Calibration');
    b_s.plotJoystickCalibration(data);
    title([p.sID ': slope ' num2str(p.slope,3) ', intcpt ' num2str(p.intercept,2) ', delay ' num2str(p.delay,2)],...
        'Interpreter', 'none');
    if savePlotOn == 1
        % save plot
        saveas(gcf, [saveResultsDir filesep p.sID '-Joystick-Calibration.fig']);
    end
    if pauseForPlots == 1
        input('updated calibration plot - press enter to continue')
    end

    % double check do we need to reshift? thought it was in here ..

    %%  post-calibration data, shift to correct units

    [ data.experiment.response, n_good]  = b_s.cleanData(data.experiment.response, p);
    disp(['#  trials = ', num2str(size( data.experiment.response, 1)),  ...
        ' Good  trials = ', num2str(n_good)]);

    %% Calibration done - move on to main model

    % shift units
    data.experiment.response = data.experiment.response-0.5;
    data.experiment.LEcontrast = data.experiment.LEcontrast-0.5;
    data.experiment.REcontrast = data.experiment.REcontrast-0.5;

    %%%%%%%%%%%%
    %% Step 1 %% Fit monocular trial portions, obtain attenuation estimates
    %%%%%%%%%%%%
    disp(' .. fitting attenuation on dropped cycles (monocular data)')

    % Create indices to select trial portions where each eye is presenting
    % zero contrast information (the "monocular" trial portions)
    % use diff() function to avoid grabbing zeros that are not part of the
    % dropped cycle (i.e. that are during regular cycles)
    diffLE = diff(data.experiment.LEcontrast,1,2);
    diffRE = diff(data.experiment.REcontrast,1,2);
    zeroContrastInLE = [zeros(size(diffLE,1),1) (diffLE == 0)];
    zeroContrastInRE = [zeros(size(diffRE,1),1) (diffRE == 0)];

    % index for all monocular trials, left OR right
    monoTrialsIndex = zeroContrastInLE | zeroContrastInRE;

    % Duplicate the participant's data struct, but replace the dichoptic
    % trial portions with NaN so they are not included in attenuation fit
    monoData = data;
    monoData.experiment.LEcontrast(~monoTrialsIndex) = NaN;
    monoData.experiment.REcontrast(~monoTrialsIndex) = NaN;
    monoData.experiment.response(~monoTrialsIndex) = NaN;

    % Run the fit - allow k to vary, fix Us to 0 and sigma to 1
    p.model ='b_s.softmax';
    % don't let certain parameters go below zero
    p.p = [1,1]; p.tau = NaN; p.m = [ 1 1];
    p.k = [1,1]; p.U = [0,0,0,0];
    p.sigma = 2; p.smax = 1; p.offset = 0;

    %freeList = {'k', 'offset'};
    % KM - not sure what offset is doing here, ask IF if I need updated b_s


    if strcmpi(analysistype, 'ns')
        % doing the normal-sighted analysis where k1=k2=1
        % skip the fit for this step

        freeList = {'offset'};
        p.costflag = 1; p = fit('b_s.getErr', p, freeList, monoData);
        disp(['   .. offset: ' num2str(round(p.offset,3)) ])


    else % do the regular old analysis
        freeList = {'k', 'offset'};
        p.costflag = 1; p = fit('b_s.getErr', p, freeList, monoData);
        disp(['   .. initial k left: ' num2str(round(p.k(1),3)) '   initial k right: ' num2str(round(p.k(2),3))])
        % Normalize relative weights
        p.k = p.k / (max(p.k));
    end

    % Grab model error
    p.costflag = 0;  [p.step1attenuationErr,~,~,~,~,~,~] = b_s.getErr(p, data);
    disp(['   .. normed k left: ' num2str(round(p.k(1),3)) '    normed k right: ' num2str(round(p.k(2),3))])
    disp(['   .. model MSE: ' num2str(round(p.step1attenuationErr, 4)) ])


    %%%%%%%%%%%%
    %% Step 2 %% Fit dichoptic trial portions, obtain suppression estimates
    %%%%%%%%%%%%
    disp(' .. fitting normalization on dichoptic data')
    % as above - replace monocular trial portions with NaN so they are not
    % included in the fit
    dichData = data;
    dichData.experiment.LEcontrast(monoTrialsIndex) = NaN;
    dichData.experiment.REcontrast(monoTrialsIndex) = NaN;
    dichData.experiment.response(monoTrialsIndex) = NaN;

    % run the fit - use k values from above, estimate U2 U3
    p.U = [0 0 0 0]; freeList = {'U(2)','U(3)'};
    p.costflag = 1;  p = fit('b_s.getErr',p, freeList, dichData);

    if p.abs == 1
        % previously only in the equation - leads to some
        % saving out as negative, previously dealt with that in
        % gatherTable but now gatherTable needs more
        % flexibility - do abs before saving individual model fit
        p.U = abs(p.U);
    end

    % grab error for this model fit, using ALL the data, and not including costs for bad param values
    p.costflag = 0;  [p.step2normalizationErr,predModel_softmax,~,data.dich,~,~,~] = b_s.getErr(p, data);
    % display outputs
    disp(['   .. U2 (right influence on left eye response): ' num2str(round(p.U(2),4))]);
    disp(['   .. U3 (left influence on right eye response): ' num2str(round(p.U(3),4))]);
    disp(['   .. sigma: ' num2str(round(p.sigma,4)) ])
    disp(['   .. model MSE: ' num2str(round(p.step2normalizationErr, 4)) ]);

    % test smax here .. ?

    %%%%%%%%%%%%
    %% Step 3 %% Re-fit k on all data
    %%%%%%%%%%%%
    if strcmpi(analysistype, 'ns')
        % doing the normal-sighted analysis where k1=k2=1
        % skip the fit for this step

    else % do the regular old analysis

        disp(' .. re-fitting k on all data')

        p.kBeforeRefit = p.k;
        freeList = {'k'};
        p.costflag = 1; p = fit('b_s.getErr', p, freeList, data);
        disp(['   .. initial k left: ' num2str(round(p.k(1),3)) '   initial k right: ' num2str(round(p.k(2),3))])

        % Normalize relative weights
        p.k = p.k / (max(p.k));
        disp(['   .. normed k left: ' num2str(round(p.k(1),3)) '    normed k right: ' num2str(round(p.k(2),3))])
        p.costflag = 0;
        [p.softmaxErr,predModel_softmax,~,data,~,~,~] = b_s.getErr(p, data);
        disp(['   .. FINAL model MSE (all data): ' num2str(round(p.softmaxErr, 4)) ])
    end
    %
    %     if strcmpi(analysistype, 'ns')
    %
    %         % try smax!
    %         freeList = {'smax'};
    %         p.costflag = 1; p = fit('b_s.getErr', p, freeList, data);
    %         disp(['   .. smax: ' num2str(p.smax,3)]);
    %         p.costflag = 0;
    %         [p.softmaxErr,predModel_softmax,~,data,~,~,~] = b_s.getErr(p, data);
    %         disp(['   .. FINAL model MSE (all data): ' num2str(round(p.softmaxErr, 4)) ])
    %     end


    %     %% grab cross-validated error
    %     tmp_p = p;
    %     % using all the data, not just dichoptic
    %     tmp_p.costflag = 0; kfoldErr = b_s.cross_calibrate(tmp_p,data, freeList);
    %     tmp_p.kfoldErr  = mean(kfoldErr);    tmp_p.kfoldStd = std(kfoldErr);
    %     p.crossvalerr = tmp_p;

    save([saveResultsDir filesep sID{i}, '_', p.condition], 'p');
    clear p
end