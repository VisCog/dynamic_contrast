%% fit_altmodels_vep
%
% This script fits data for the VEP experiment (either the vep data itself
% or the corresponding psychophysical data) using a variety of models that
% are not our regular model.
%
%%% NOTE mean/max/minkowski modeling has been merged with our main fit
%%% script, fit_vep_dynamic_contrast. I'm keeping this script for
%%% historical reasons in case I need to double-check anything but consider
%%% it obsolete and I'll delete it soon once I confirm no more information
%%% is needed out of here.

clear
close all

%% analysis settings

%%% datatype:
datatype = 'vep_psychophysics'; %'vep' or 'vep_psychophysics'

%%% other settings
useAbs      = 1;
clean_range = 0.3;
condition   = 'congruent';

%%% default calibration settings
startT      = 1;
slope       = 1;
intercept   = 'mean';
delay       = 0.5;
calibfunction = 'delay + scale';
calibfreelist = {'slope', 'intercept', 'delay'};



%%% select/make folders
if strcmp(datatype, 'vep')  %%%% electrophysiological data
    rawDataDir = [cd filesep 'output_vep'];
    saveResultsDir = [cd filesep 'fitdata_vep'  filesep 'model_fits_altmodels'];

elseif strcmp(datatype, 'vep_psychophysics')  %%%% psychophysical data
    rawDataDir = [cd filesep 'output_vep_psychophysics'];
    saveResultsDir = [cd filesep 'fitdata_vep_psychophysics'  filesep 'model_fits_altmodels'];
end

plotResultsDir = strrep(saveResultsDir, 'model_fits', 'plots');

if ~isfolder(saveResultsDir)
    mkdir(saveResultsDir);
end
if ~isfolder(plotResultsDir)
    mkdir(plotResultsDir);
end

%% go
subjectList_vep; % puts variable called sID in workspace

% add mean and max example (only exist in VEP data)
sID = [sID 'meanExample', 'maxExample'];

disp( ' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp([' %% datatype: ' datatype])
disp([' %% condition: ' condition])
disp( ' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
i=1

for i = 1:length(sID)

    disp(['   Working on ' sID{i} ' (' num2str(i) ' of ' num2str(length(sID)) ')'])
    p.sID           = sID{i};
    p.abs           = useAbs;
    p.clean_range   = clean_range;
    p.condition     = condition;

    % load raw data:
    file = dir([rawDataDir filesep sID{i} '*-',p.condition, '.mat']);
    if isempty(file), disp('(doesn''t exist)'), continue, end
    load([file.folder filesep file.name]);
    % puts a struct called congruentVep or orthogonalVep
    %  (or congruentMotor, orthogonalMotor) into the workplace

    if strcmp(datatype, 'vep'), data = eval([p.condition, 'Vep']);
    elseif strcmp(datatype, 'vep_psychophysics'), data = eval([p.condition, 'Motor']);
    end

    %%%%% if psychophysics data, select only trials on which they responded
    if strcmp(datatype, 'vep_psychophysics')
        data.experiment.BEcontrastStart = data.experiment.BEcontrastStart(data.conditionInfo.joyUsedIndex,:);
        data.experiment.binoResponseStart = data.experiment.binoResponseStart(data.conditionInfo.joyUsedIndex,:);
        data.experiment.BEcontrastEnd = data.experiment.BEcontrastEnd(data.conditionInfo.joyUsedIndex,:);
        data.experiment.binoResponseEnd = data.experiment.binoResponseEnd(data.conditionInfo.joyUsedIndex,:);
        data.experiment.LEcontrast = data.experiment.LEcontrast(data.conditionInfo.joyUsedIndex,:);
        data.experiment.REcontrast = data.experiment.REcontrast(data.conditionInfo.joyUsedIndex,:);
        data.experiment.response = data.experiment.response(data.conditionInfo.joyUsedIndex,:);
    end

    % b_s expects the binocular data to be inside data.experiment.binoResponse:
    concatContrast  = [data.experiment.BEcontrastStart data.experiment.BEcontrastEnd]; % ntrials x duration
    concatResponse  = [data.experiment.binoResponseStart data.experiment.binoResponseEnd];
    data.experiment.binoS = concatContrast(1, :)-0.5;     % put everything in  -0.5 to 0.5 units
    data.experiment.binoResponse = concatResponse - 0.5;
    nan_indx = round(length(data.experiment.binoS)/2-3:length(data.experiment.binoS)/2+3);
    data.experiment.binoResponse(:, nan_indx) = NaN; % block out switchy bit when concatenating trial

    % apply clean_range criterion
    [data.experiment.binoResponse, p.n_good]  = b_s.cleanData(data.experiment.binoResponse, p);
    disp(['    # Calib trials total = ', num2str(size( data.experiment.binoResponse, 1)),  ...
        ', # used = ', num2str(p.n_good)]);

    %%%%%% calibration
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
    p.delay = delay;
    p.junk = 0;
    binoMean.gvals = [];
    p.dt = diff(data.binocular.t(1:2));
    data.binocular.t = 0:p.dt:p.dt*size(data.experiment.binoResponse, 2)-p.dt;
    p.joystickfunction = calibfunction;
    p.penalizeDelay = 3;

    % fit
    p.costflag = 1; p = fit('b_s.getErrBinoMean', p, calibfreelist, data);
    % error
    p.costflag = 0;
    [~, p.errMean, data] = b_s.getErrBinoMean(p,data);
    [~, p.errInd, ~] = b_s.getErrBinoInd(p,data);
    disp(['   .. delay: ' num2str(p.delay,4) '    intercept: ' num2str(p.intercept,3) '    slope: ' num2str(p.slope,3)])
    disp(['   .. MSE:   for mean: ' num2str(p.errMean,5) '    for individual trials: ' num2str(p.errInd,5) ]);

    %%%%%% fit
    data.experiment.response = data.experiment.response-0.5; % shift units
    data.experiment.LEcontrast = data.experiment.LEcontrast-0.5;
    data.experiment.REcontrast = data.experiment.REcontrast-0.5;

    % non-fitting models - get err
    [p.simpleAverageErr, meanModelPrediction] = b_s.simpleAverage(p, data);
    [p.simpleMaxErr, maxModelPrediction] = b_s.simpleMax(p, data);

    %     % simple softmaxmodel
    %     p.model = 'b_s.simpleSmax';
    %     p.smax = 1; freeList = 'smax';
    %     p.costflag = 1; p = fit('b_s.getErr', p, freeList, data);
    %     p.costflag = 0; [p.simpleSmaxErr,~,~,~,~,~,~] = b_s.getErr(p, data);

    disp('   Errors:')
    disp(['   - ' num2str(p.simpleAverageErr,3) ' mean model'])
    disp(['   - ' num2str(p.simpleMaxErr,3) ' max model'])
    %     disp(['   - ' num2str(p.simpleSmaxErr,3) ' simple max model (smax = ' num2str(p.smax,2) ')'])

    save([saveResultsDir filesep sID{i}, '_', p.condition], 'p');



    clear p
end