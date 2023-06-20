% fit_psychophysics_dynamic_contrast
%
% This script fits all the things (data for Scientific Reports submission)
% Participant IDs are listed in current directory's file called
% subjectList_psychophysics.m
%
% This loads each one up, fits it, and saves fit parameters inside
% current directory/fitdata_psychophyics/model_fits/ as [subjectID.mat]
%
% Current implementation: calibrate joystick, fit monocular trial portions
% for estimate of attenuation (k), then fit dichoptic (exclude monocular)
% trial portions for estimate of suppression/inter-ocular normalization and
% sigma. NEW: re-fit k to all data for re-adjustment


%% Paths
% addpath(genpath(cd))

clear
rawDataDir = [cd filesep 'output_psychophysics']; % where the raw data live
saveResultsDir = [cd filesep 'fitdata_psychophysics' filesep 'model_fits'];

%% Load subject IDs

% load participants we're fitting:
subjectList_psychophysics; % puts variable called sID in workspace


%% Individual subjects
ct_k = 0;

for i = 1:length(sID) % replace this with the sID # to run only 1 person
    
    % we clear p at the end of each person so this needs to be inside the loop
    p.mid_range_flag = 0;   if p.mid_range_flag==1;  disp('p.mid_range_flag = 1 is this correct?'); end
    p.clean_range  = 0.5;   % only calibrate or fit data where there's this must range in the data, p.clean_range = 0, uses all data
    p.abs = 1;              % use absolute-value version of equations -- e.g. abs(p.U) -- to prevent negative numbers
    p.Kfold = 4;

    p.sID = sID{i};
    p.when = datestr(now);

    disp('++++++++++++++++++++++++++++++++++++++++++++++');
    disp(['Working on ' sID{i} ' (' num2str(i) ' of ' num2str(length(sID)) ')'])

    % load raw data:
    file = dir([rawDataDir filesep sID{i} '*_run*.mat']);
    data = load([file.folder filesep file.name]);

    % ---- reduced data methods here ------
    reducedTrialAnalysis = 0;
    if reducedTrialAnalysis == 1 % reduced trial analysis
        disp(' Conducting reduced trial length analysis')
        % the trials are 48 seconds long and the mono trial portion happens in
        % either the first 24-sec or the second 24-sec portion. We cut trial length
        % in half by determining which of the trial halves contains the mono
        % portion and use that half only, while excluding the trial half that
        % contains dichoptic presenetation only (this is such a simple concept but
        % hard to describe in words)

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

        % index for - does the first half / second half of the trial contain mono?
        nHalfSamples = size(monoTrialsIndex,2)/2;
        nSamples = size(monoTrialsIndex,2);
        iFirstHalf = any(monoTrialsIndex(1:size(monoTrialsIndex,1), 1:nHalfSamples), 2);
        iSecondHalf = any(monoTrialsIndex(1:size(monoTrialsIndex,1), nHalfSamples+1:nSamples), 2);

        iBoth = iFirstHalf == iSecondHalf;

        % trial indices:
        iFirstHalf(iBoth) = 0;
        iSecondHalf(iBoth) = 0;
        reducedData = data;

        reducedData.experiment.response = [];
        reducedData.experiment.REcontrast = [];
        reducedData.experiment.LEcontrast = [];

        reducedData.experiment.response = data.experiment.response(iFirstHalf, 1:nHalfSamples);
        reducedData.experiment.REcontrast = data.experiment.REcontrast(iFirstHalf, 1:nHalfSamples);
        reducedData.experiment.LEcontrast = data.experiment.LEcontrast(iFirstHalf, 1:nHalfSamples);

        reducedData.experiment.response = [reducedData.experiment.response;...
            data.experiment.response(iSecondHalf, nHalfSamples+1:nSamples)];
        reducedData.experiment.REcontrast = [reducedData.experiment.REcontrast;...
            data.experiment.REcontrast(iSecondHalf, nHalfSamples+1:nSamples)];
        reducedData.experiment.LEcontrast = [reducedData.experiment.LEcontrast;...
            data.experiment.LEcontrast(iSecondHalf, nHalfSamples+1:nSamples)];

        % ---- end reduced data method, continue using new structs ------
        data = reducedData;
    end


    %%%%%%%%%%%%
    %% Step 0 %% Fit joystick
    %%%%%%%%%%%%

    disp(' .. calibrating joystick position')

    % Fits are being done on the average position (across all 28 runs) for the
    % 14-second binocular period prior to the dichoptic onset.
    p.joystickfunction = 'delay + scale';
    p.startT = 2;
    p.delay = 1; p.penalizeDelay = 2; % penalize delays larger than this, in seconds
    p.intercept = 0.3; p.junk = 0;
    p.slope = 1; binoMean.gvals = [];
    p.dt = diff(data.binocular.t(1:2));
    data.experiment.binoS = (sin(2*pi*data.binocular.t/data.binocular.period)+1)/2 - 0.5; % stimulus time course, -0.5 to 0.5

    %cut out first 4s and move to  -0.5 to 0.5 scale
    data.experiment.binoResponse(:,data.binocular.t < 4) = NaN;
    data.experiment.binoResponse = data.experiment.binoResponse-0.5;
    [data.experiment.binoResponse, p.n_good]  = b_s.cleanData(data.experiment.binoResponse, p);
    p.intercept = -2*nanmean(data.experiment.binoResponse(:));
    p.costflag = 1; p = fit('b_s.getErrBinoMean', p, {'delay', 'slope', 'intercept'}, data);

    if strcmpi('AM_RE_G_RF_20', sID{i})
        % this person gets weirdly high slope which is the crux of later
        % fit problems (vU alues of 10^5+!!) using default 1 leads
        % to reasonable values later, and calibration MSE isn't different,
        % so opting for this
        p.slope = 1;
    end

    % Calculate the error, both for mean joystick position and for
    % indivdual trials
    p.costflag = 0;  [~, p.calibErr, data] = b_s.getErrBinoInd(p,data);

    % display fits and error
    disp(['   .. delay: ' num2str(round(p.delay,2)) ...
        '    intercept: ' num2str(round(p.intercept,2)) ...
        '    slope: ' num2str(round(p.slope,2)) ...
        '    MSE: ' num2str(round(p.calibErr,4))])
    disp(['   .. n_good: ' num2str(p.n_good)]);

    %figure(1); clf;  set(gcf, 'Name', 'Calibration');
    %p.costflag = 0; [err, err_noCost, data] = b_s.getErrBinoMean(p,data);
    %b_s.plotJoystickCalibration(data);
    %[err, err_noCost, data] = b_s.getErrBinoInd(p,data);


    %%  post-calibration data, shift to correct units
    data.experiment.response = data.experiment.response-0.5;
    data.experiment.LEcontrast = data.experiment.LEcontrast-0.5;
    data.experiment.REcontrast = data.experiment.REcontrast-0.5;
    [data.experiment.response, p.n_good]  = b_s.cleanData(data.experiment.response, p);

    % Notes on joystick calibration
    % the original binocular joystick response is in 0-1 units, but because
    % sometimes people don't get to 1, the fits allow units above 1

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

    % move to the right units

    % Run the fit - allow k to vary, fix Us to 0 and sigma to 1
    p.model ='b_s.softmax';
    p.m = [1,1];  p.U = [0,0,0,0];   p.sigma = 1;   p.smax = 1; p.tau = NaN;
    p.k = [1,1];
    freeList = {'k'};% don't let certain parameters go below zero
    p.costflag = 1; p = fit('b_s.getErr', p, freeList, monoData);

    if p.abs == 1
        % previously only in the equation - leads to some saving out as
        % negative - do abs before saving individual model fit so doesn't
        % have to be accounted for later (previously in gatherTable)
        p.k = abs(p.k);
    end

    p.costflag = 0;  [errLR,predModel,~,~,~,~,~] =    b_s.getErr(p, monoData);
    disp(['   .. initial k left: ' num2str(round(p.k(1),4)) '   initial k right: ' ,  num2str(round(p.k(2),4))])

    % Normalize relative weights
    p.k = p.k/(max(p.k));
    disp(['   .. normed k left: ' num2str(round(p.k(1),4)) '    normed k right: ' num2str(round(p.k(2),4))])

    % Grab model error
    p.costflag = 0;  [p.step1attenuationErr,~,~,~,~,~,~] = b_s.getErr(p, data);
    disp(['   .. model MSE: ' num2str(round(p.step1attenuationErr, 4)) ])

    %%%%%%%%%%%%
    %% Step 2 %% Fit dichoptic trial portions using softmax model, obtain suppression estimates
    %%%%%%%%%%%%
    disp(' .. fitting normalization on dichoptic data')

    % as above - replace monocular trial portions with NaN so they are not
    % included in the fit
    dichData = data;
    dichData.experiment.LEcontrast(monoTrialsIndex) = NaN;
    dichData.experiment.REcontrast(monoTrialsIndex) = NaN;
    dichData.experiment.response(monoTrialsIndex) = NaN;

    % run the fit - use k values from above, estimate U2 U3 and sigma

    p.usegridsearch = 1; % set to 1 for gridsearch, 0 for no gridsearch (saved in p for posterity)

    p.costflag = 1;
    if p.usegridsearch == 1 
        freeList = {'U(2)','U(3)','sigma'};
        p.U =[0 0 0 0];
        p.sigma = 1;
        gridParams = {'U(2)','U(3)','sigma'}; % list of parameters to grid
        gridList = {[0:0.25:3], [0:0.25:3], [0:0.1:1]};   % list of grid vectors for each parameter in gridParams
        [pBest,errBest] = b_s.gridsearch('b_s.getErr',p,gridParams,gridList,dichData);
        p = pBest;
        p.costflag = 0; p = fit('b_s.getErr',p, freeList, dichData);

    else
        freeList = {'U(2)','U(3)','sigma'};
        p.U = [0 0 0 0]; 
        p.sigma = 1;
        p = fit('b_s.getErr',p, freeList, dichData);
    end

    if p.abs == 1 % equation used absolute values - convert
        p.U = abs(p.U);
        p.sigma = abs(p.sigma);
    end

    % grab error for this model fit, using ALL the data, and not including costs for bad param values
    p.costflag = 0;  [p.step2normalizationErr,predModel_softmax,~,data.dich,~,~,~] = b_s.getErr(p, data);

    % display outputs
    disp(['   .. U2 (right influence on left eye response): ' num2str(round(p.U(2),4))]);
    disp(['   .. U3 (left influence on right eye response): ' num2str(round(p.U(3),4))]);
    disp(['   .. sigma: ' num2str(round(p.sigma,4)) ])
    disp(['   .. model MSE: ' num2str(round(p.step2normalizationErr, 4)) ])

    %%%%%%%%%%%%
    %% NEW! Step 3 %% fit ALL data to re-adjust k
    %%%%%%%%%%%%
    disp(' .. re-fitting attenuation on all data')

    p.kBeforeRefit = p.k;
    freeList = {'k'};
    p.costflag = 1; p = fit('b_s.getErr', p, freeList, data);

    if p.abs == 1% equation used absolute values - convert
        p.k = abs(p.k);
    end

    % Normalize relative weights
    p.k = p.k / (max(p.k));

    % grab error for this model fit, using ALL the data, and not including costs for bad param values
    p.costflag = 0;  [p.softmaxErr,predModel_softmax,~,data.dich,~,~,~] = b_s.getErr(p, data);

    % display outputs
    disp(['   .. re-fit k left: ' num2str(round(p.k(1),4)) '    re-fit k right: ' num2str(round(p.k(2),4))])
    disp(['   .. model MSE: ' num2str(round(p.softmaxErr, 4)) ])

    %% Alternative models section
    alternateModels = 1;
    if alternateModels == 1
        %% grab cross-validated error
        disp('  ... cross-validated error')
        orig_p = p;
        tmp_p = orig_p;
        % using all the data, not just dichoptic
        tmp_p.costflag = 0; kfoldErr = b_s.cross_calibrate(tmp_p, data, freeList);
        tmp_p.kfoldErr  = mean(kfoldErr);    tmp_p.kfoldStd = std(kfoldErr);
        p.softmax = tmp_p;

        %% alternative models
        % begin with models that don't have any free parameters that need
        % fitting, what's returned are not  kfold values, since kfold makes no
        % sense without free parameters

        disp('  ... simple models')

        [p.rivalry.Err, ~] = b_s.rivalry(p, data);
        [p.simpleAverage.Err, ~] = b_s.simpleAverage(p, data);
        [p.simpleMax.Err, ~] = b_s.simpleMax(p, data);
        [p.meanJoystick.Err, ~] = b_s.meanJoystick(p, data);
        [p.DualMeanJoystick.Err, ~] = b_s.DualMeanJoystick(p, data);

        % variable tau
        disp('  ... a + n with tau')

        tmp_p = orig_p;   tmp_p.model = 'b_s.softmax'; % this time with free tau
        freeList = {'U(2)','U(3)','sigma', 'tau'};     % using all the data, not just dichoptic
        tmp_p.costflag = 1; tmp_p.tau = 50; tmp_p = fit('b_s.getErr',tmp_p, freeList, data);
        tmp_p.costflag = 0; kfoldErr = b_s.cross_calibrate(tmp_p,data, freeList);    % grab cross-validated error
        tmp_p.kfoldErr  = mean(kfoldErr);    tmp_p.foldStd = std(kfoldErr);
        p.softmax_tau = tmp_p;

        % weighted average
        disp('  ... weighted avg')
        tmp_p = orig_p;  tmp_p.model = 'b_s.weightedAverage';
        tmp_p.wa = 0.5;  freeList = {'wa'};
        tmp_p.costflag = 1; tmp_p= fit('b_s.getErr',tmp_p, freeList, data);
        tmp_p.costflag = 0;  kfoldErr = b_s.cross_calibrate(tmp_p, data, freeList);
        tmp_p.kfoldErr  = mean(kfoldErr);    tmp_p.kfoldStd = std(kfoldErr);
        p.weightedAverage = tmp_p;

        disp('  ... D & S 2006')
        tmp_p = orig_p;  tmp_p.model = 'b_s.ds2006'; % Ding & Sperling 2006
        tmp_p.e = [0.5 0.5 0.5 0.5]; tmp_p.tau = 50;
        freeList = {'e'};
        tmp_p.costflag = 1; tmp_p = fit('b_s.getErr',tmp_p, freeList, data);
        tmp_p.costflag = 0;   kfoldErr = b_s.cross_calibrate(tmp_p,data, freeList);
        tmp_p.kfoldErr  = mean(kfoldErr);    tmp_p.kfoldStd = std(kfoldErr);
        p.ds2006 = tmp_p;

        disp('  ... D & S 2006 with tau')
        tmp_p = orig_p;  tmp_p.model = 'b_s.ds2006'; % Ding & Sperling 2006, allowing tau to vary
        tmp_p.e = [0.5 0.5 0.5 0.5]; tmp_p.tau = 50;
        freeList = {'e', 'tau'};
        tmp_p.costflag = 1; tmp_p = fit('b_s.getErr',tmp_p, freeList, data);
        tmp_p.costflag = 0;   kfoldErr = b_s.cross_calibrate(tmp_p,data, freeList);
        tmp_p.kfoldErr  = mean(kfoldErr);    tmp_p.kfoldStd = std(kfoldErr);
        p.ds2006_tau = tmp_p;

        disp('  ... BMG 2007')
        tmp_p = orig_p; tmp_p.model = 'b_s.bmg2007'; %Baker, Meese and Georgeson 2007
        tmp_p.m = 1.28;   tmp_p.pq = [7.99 6.59]; tmp_p.Z = 0.076;
        tmp_p.S = .6;  tmp_p.w= [1,1];% starting parameters based on Table 1
        freeList = {'w', 'S'};
        tmp_p.costflag = 1; tmp_p = fit('b_s.getErr',tmp_p, freeList, data);
        tmp_p.costflag = 0;   kfoldErr = b_s.cross_calibrate(tmp_p,data, freeList);
        tmp_p.S = abs(tmp_p.S); tmp_p.w = abs(tmp_p.w); % forced to be abs in the function
        tmp_p.kfoldErr  = mean(kfoldErr);    tmp_p.kfoldStd = std(kfoldErr);
        p.bmg2007 = tmp_p;

    end



    % save this
    if reducedTrialAnalysis == 1
        save([saveResultsDir filesep sID{i} '-reduced'], 'p');
    elseif p.mid_range_flag == 1
        save([saveResultsDir filesep sID{i} '-midrange'], 'p');
    elseif alternateModels == 1
        save([saveResultsDir filesep sID{i} '-altmodels'], 'p');
    else
        save([saveResultsDir filesep sID{i} '-regular'], 'p');
    end

    if min(p.k) == max(p.k)
        ct_k = ct_k+1;
    end

    % clean save file size
    clear binoMean data dichData diffLE diffRE doTwice err err_noCost errLR file ...
        freeList kfoldErr monoData monoTrialsIndex n_good orig_p p predModel ...
        predModel_softmax tmp_p zeroContrastInLE zeroContrastInRE ...
        p_fitk1 p_fitk2 iBoth iFirstHalf iSecondHalf reducedData ...
        tBoth tFirstHalf tSecondHalf
end


