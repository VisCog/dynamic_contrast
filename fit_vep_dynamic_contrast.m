%% fit_vep_dynamic_contrast
%
% This script fits  the VEP data (all trials, not separated by
% whether the joystick was used or not) or the psychophysics data (only
% trials where the joystick was used)
% Participant IDs are listed in current directory's file called subjectList_vep.m
%
% <old words deleted> This header needs to be updated but the code is
% fairly well commented at this point


%%
clear
close all


%%%%% analysis details: %%%%%%
analysistype    = '';
% analysistype options:
%   ''          leave blank for standard 
%   'ns'        set k1 = k2 = 1, fit Us, normal-sighted only


datatype    = 'vep_psychophysics';  %'vep' or 'vep_psychophysics'
condition   = 'congruent';
savePlotOn      = 1;        % if 1, saves plots
pauseForPlots   = 0;        % if 1, waits for you to press enter after each plot before continuing

%%%%% calibration settings/defaults %%%%%
clean_range = 0.3; % calibration: use trials with this min range
startT      = 1;
slope       = 1;
intercept   = 'mean';
calibFreeList = {'slope', 'intercept', 'delay'};
joyfunction = 'delay + scale';
delay       = 0.5;  % intial guess for delay b/w stimulus and motor (sec)
penalizeDly = 3;    % penalize delays longer than this (sec)

%%%%% model-fitting settings/defaults %%%%%
modelStr    = 'b_s.softmax'; % default to softmax model
useAbs =    1;   % use absolute value formulae
pp =        [1,1];
ptau =      NaN;
pm =        [1 1];
psmax =     1;
monoFreeList = {'k'}; % free parameters for monocular fitting stage (overwritten later if analysistype='ns')
poffset =   0;
pk =        [1 1];
dichFreeList = {'U(2)','U(3)', 'sigma'}; % free parameters for dichoptic fitting stage
pU =        [0,0,0,0];
psigma =    1;




%% set up

subjectList_vep; % puts variable called sID in workspace

switch lower(datatype)
    case 'vep'
        rawDataDir = [cd filesep 'output_vep']; % where the raw data live
        saveResultsDir = [cd filesep 'fitdata_vep'  filesep 'model_fits'];
    case 'vep_psychophysics'
        rawDataDir = [cd filesep 'output_vep_psychophysics']; % where the raw data live
        saveResultsDir = [cd filesep 'fitdata_vep_psychophysics'  filesep 'model_fits'];
    otherwise
        disp([' datatype ' datatype ' not defined (typo?)'])
end

if strcmpi(analysistype, 'ns')
    disp('%%%%%%%%%%%%%% You''re running the NS-only, fix k1=k2=1 analysis')
    saveResultsDir = strrep(saveResultsDir, 'model_fits', 'model_fits_fixedk');
    disp('%%%%%%%%%%%%%% If that''s unexpected, change analysistype in the code')
    % select normally-sighted participants only:
    idx = cellfun(@(x) strcmpi(x(1:2),'NS'), sID);
    sID = sID(idx);
end

%% Individual subjects

for i = 1:length(sID) % replace this with the sID # to run only 1 person
    % we clear p at the end of each individual so these need to be inside
    % the loop
    p.sID = sID{i};
    p.abs = useAbs;
    p.clean_range  = clean_range; % only calibrate or fit data where there's this must range in the data, p.clean_range = 0, uses all data
    p.condition = condition;

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
    data.experiment.binoS        = concatContrast(1,:) - 0.5;
    data.experiment.binoResponse = concatResponse - 0.5;
    nan_indx = round(length(data.experiment.binoS)/2-3:length(data.experiment.binoS)/2+3);
    data.experiment.binoResponse(:, nan_indx) = NaN; % block out switchy bit when concatenating trial

    [data.experiment.binoResponse, p.n_good]  = b_s.cleanData(data.experiment.binoResponse, p);

    disp(['# Calib trials = ', num2str(size( data.experiment.binoResponse, 1)),  ...
        ' Good Calib trials = ', num2str(p.n_good)]);

    % Calibration defaults
    p.startT = startT;
    p.slope = slope;
    if ischar(intercept)
        switch lower(intercept)
            case 'mean'
                p.intercept = -2*nanmean(data.experiment.binoResponse(:)); %#ok<*NANMEAN>
            otherwise, disp(['intercept case ' intercept ' undefined']);% break
        end
    elseif isnumeric(intercept)
        p.intercept = intercept;
    else, disp('intercept undefined'); %break
    end
    p.junk = 0;
    binoMean.gvals = [];
    p.dt = diff(data.binocular.t(1:2));
    data.binocular.t = 0:p.dt:p.dt*size(data.experiment.binoResponse, 2)-p.dt;
    p.joystickfunction = joyfunction;
    p.delay = delay;
    p.penalizeDelay = penalizeDly; % delay  penalization in seconds

    % fit:
    p.costflag = 1; p = fit('b_s.getErrBinoMean', p, calibFreeList, data);

    % Calculate error, for mean joystick position and for individual trials
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

    %% Calibration done - move on to main model
    [ data.experiment.response, n_good]  = b_s.cleanData(data.experiment.response, p);
    disp(['#  trials = ', num2str(size( data.experiment.response, 1)),  ...
        ' Good  trials = ', num2str(n_good)]);

    % shift units
    data.experiment.response = data.experiment.response-0.5;
    data.experiment.LEcontrast = data.experiment.LEcontrast-0.5;
    data.experiment.REcontrast = data.experiment.REcontrast-0.5;


    %%%%%%%%%%%%
    %% Step 1 %% Fit monocular trial portions
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
    p.model = modelStr;
    p.p = pp;
    p.tau = ptau;
    p.m = pm;
    p.k = pk;
    p.U = pU;
    p.sigma = psigma;
    p.smax = psmax;
    p.offset = poffset;

    switch analysistype
        case 'ns'
            % doing the normal-sighted analysis where k1=k2=1
            % skip the fit for this step, do offset only
            freeList = {'offset'};
            p.costflag = 1; p = fit('b_s.getErr', p, freeList, monoData);
            disp(['   .. offset: ' num2str(round(p.offset,3)) ])

        otherwise % do the regular old analysis
            p.costflag = 1; p = fit('b_s.getErr', p, monoFreeList, monoData);
            disp(['   .. initial k left: ' num2str(round(p.k(1),3)) '   initial k right: ' num2str(round(p.k(2),3))])
            % Normalize relative weights
            p.k = p.k / (max(p.k));
    end

    % Grab model error (all data)
    p.costflag = 0;  [p.step1attenuationErr,~,~,~,~,~,~] = b_s.getErr(p, data);
    disp(['   .. normed k left: ' num2str(round(p.k(1),3)) '    normed k right: ' num2str(round(p.k(2),3))])
    disp(['   .. offset: ' num2str(round(p.offset,3)) ])
    disp(['   .. model MSE: ' num2str(round(p.step1attenuationErr, 4)) ])


    %%%%%%%%%%%%
    %% Step 2 %% Fit dichoptic trial portions
    %%%%%%%%%%%%
    disp(' .. fitting normalization on dichoptic data')
    % as above - replace monocular trial portions with NaN so they are not
    % included in the fit
    dichData = data;
    dichData.experiment.LEcontrast(monoTrialsIndex) = NaN;
    dichData.experiment.REcontrast(monoTrialsIndex) = NaN;
    dichData.experiment.response(monoTrialsIndex) = NaN;

    % run the fit
    p.costflag = 1;  p = fit('b_s.getErr',p, dichFreeList, dichData);

    if p.abs == 1
        % previously only in the equation - leads to some
        % saving out as negative, previously dealt with that in
        % gatherTable but now gatherTable needs more
        % flexibility - do abs before saving individual model fit
        p.U = abs(p.U);
    end

    % grab error (all data)
    p.costflag = 0;  [p.step2normalizationErr,~,~,~,~,~,~] = b_s.getErr(p, data);
    % display outputs
    disp(['   .. U2 (right influence on left eye response): ' num2str(round(p.U(2),4))]);
    disp(['   .. U3 (left influence on right eye response): ' num2str(round(p.U(3),4))]);
    disp(['   .. sigma: ' num2str(round(p.sigma,4)) ])
    disp(['   .. model MSE: ' num2str(round(p.step2normalizationErr, 4)) ]);

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
        [p.softmaxErr,predModel_softmax,~,~,~,~,~] = b_s.getErr(p, data);
        p.predModel_softmax = cell2mat(predModel_softmax');
        disp(['   .. FINAL model MSE (all data): ' num2str(round(p.softmaxErr, 4)) ])
    end


%         %% grab cross-validated error
%         % ask IF if that freelist is correct for cross-val error (it's
%         % whatever freeList was set to above on the last step)
%         tmp_p = p;
%         % using all the data, not just dichoptic
%         tmp_p.costflag = 0; kfoldErr = b_s.cross_calibrate(tmp_p,data, freeList);
%         tmp_p.kfoldErr  = mean(kfoldErr);    tmp_p.kfoldStd = std(kfoldErr);
%         p.crossvalerr = tmp_p;

    %% mean/max/minkowski models
    % three models: mean(L,R)  max(L,R), and single-parameter minkowski
    % equiation.
    %
    % unlike above where we analyze individual trials, here we compare the
    % participant's response to each of these three models by looking at
    % their MEAN response over all trials. Additionally, the mean response
    % is scaled to circumvent regression to the mean problems.

    % some of this repeats steps from above but want to keep it modular in
    % case things need to differ:
    minkData = data;

    % we don't need to split into left fast/right fast trials to obtain the
    % mean since we treat the eyes as equally-weighted, however we do need
    % to NaN out the monocular portions
    diffLE = diff(minkData.experiment.LEcontrast,1,2);
    diffRE = diff(minkData.experiment.REcontrast,1,2);
    zeroContrastInLE = [zeros(size(diffLE,1),1) (diffLE == 0)];
    zeroContrastInRE = [zeros(size(diffRE,1),1) (diffRE == 0)];

    % index for all monocular trials, left OR right
    monoTrialsIndex = zeroContrastInLE | zeroContrastInRE;

    % get rid of them for the averaging:
    minkData.experiment.LEcontrast(monoTrialsIndex) = NaN;
    minkData.experiment.REcontrast(monoTrialsIndex) = NaN;
    minkData.experiment.response(monoTrialsIndex) = NaN;

    % Use the getErr function to obtain the calibrated joystick position
    % for this data (we don't need the model prediction for this part)
    [~, ~, ~, respJoyCalib, ~, ~, ~] = b_s.getErr(p, minkData);


    respJoyCalib = cell2mat(respJoyCalib'); %reshape to trials x timepoints
    mn_resp = nanmean(respJoyCalib);

    % generate the stimuli timecourses:
    fastSin = ((sin(2*pi*data.t/6)+1)/2);
    slowSin = ((sin(2*pi*data.t/8)+1)/2);

    % after running through calibration there will be nans at the end due
    % to delay - these mess up the minkowski fit function - chop off NaNs
    nan_idx = isnan(mn_resp);
    mn_resp = mn_resp(~nan_idx);
    fastSin = fastSin(~nan_idx);
    slowSin = slowSin(~nan_idx);

    mn_resp_scaled = rescale(mn_resp, 0, 1);% range between 0 and 1


    S = [fastSin;slowSin]';

    % generate the predicted model values:
    meanModelPrediction = mean(S,2)';
    maxModelPrediction = max(S, [], 2)';

    % generate fit/error values
    p.n = 1;
    p.costflag = 1; p = fit('b_s.minkowski', p, {'n'}, S, mn_resp_scaled'); % fit minkowski
    p.costflag = 0;
    [mink_err, ~, mink_pred] = b_s.minkowski(p, S, mn_resp_scaled');


    % save these details for later in same format as above
    p.meanResponseScaled = mn_resp_scaled;
    p.predModel_meanModel = meanModelPrediction;
    p.meanModelErr = sum((meanModelPrediction-mn_resp_scaled).^2)/length(mn_resp_scaled);
    p.predModel_maxModel = maxModelPrediction;
    p.maxModelErr = sum((maxModelPrediction-mn_resp_scaled).^2)/length(mn_resp_scaled);
    p.predModel_minkModel = mink_pred';
    p.minkModelErr = mink_err;

    % display outputs
    disp(' .. simpler model fits')

    disp(['   .. mean model MSE: ' num2str(round(p.meanModelErr,4))]);
    disp(['   .. max model MSE: ' num2str(round(p.maxModelErr,4))]);
    disp(['   .. Minkowski parameter: ' num2str(round(p.n,2))]);
    disp(['   .. Minkowski model MSE: ' num2str(round(p.minkModelErr,4))]);

    %% save it

    save([saveResultsDir filesep sID{i}, '_', p.condition], 'p');
    clear p
end