clear

%datatype = 'vep';
datatype = 'vep_psychophysics';

group = 'AM';


%% Load relevant data
subjectList_vep;

idx = cellfun(@(x) strcmpi(x(1:2),group), sID); % select group
sID = sID(idx);

i = 1
sCond = 'congruent';

switch lower(datatype)

    case 'vep'
        vepDataFile = dir([cd filesep 'output_vep' filesep sID{i} '-vep-' sCond '.mat']);
        load([vepDataFile.folder filesep vepDataFile.name]); % loads congruentVep
        data = congruentVep; clear congruentVep vepDataFile;
        

        vepFitFile = dir([cd filesep 'fitdata_vep' filesep 'model_fits' filesep sID{i} '_' sCond '.mat']);     % vep
        load([vepFitFile.folder filesep vepFitFile.name]); % loads p
        modelfit = p; clear p vepFitFile;

    case 'vep_psychophysics'
        psyDataFile = dir([cd filesep 'output_vep_psychophysics' filesep sID{i} '-motor-' sCond '.mat']);
        load([psyDataFile.folder filesep psyDataFile.name]); % loads congruentMotor;
        data = congruentMotor; clear congruentMotor psyDataFile;
        
        psyFitFile = dir([cd filesep 'fitdata_vep_psychophysics' filesep 'model_fits' filesep sID{i} '_' sCond '.mat']); %psychophysics
        load([psyFitFile.folder filesep psyFitFile.name]); % loads p
        modelfit = p; clear p psyFitFile;

    otherwise
        disp([datatype ' is an undefined datatype'])
end

% reduce data to only joystick-used trials
    if strcmp(datatype, 'vep_psychophysics')
        data.experiment.BEcontrastStart = data.experiment.BEcontrastStart(data.conditionInfo.joyUsedIndex,:);
        data.experiment.binoResponseStart = data.experiment.binoResponseStart(data.conditionInfo.joyUsedIndex,:);
        data.experiment.BEcontrastEnd = data.experiment.BEcontrastEnd(data.conditionInfo.joyUsedIndex,:);
        data.experiment.binoResponseEnd = data.experiment.binoResponseEnd(data.conditionInfo.joyUsedIndex,:);
        data.experiment.LEcontrast = data.experiment.LEcontrast(data.conditionInfo.joyUsedIndex,:);
        data.experiment.REcontrast = data.experiment.REcontrast(data.conditionInfo.joyUsedIndex,:);
        data.experiment.response = data.experiment.response(data.conditionInfo.joyUsedIndex,:);
        % update fast eye trials
        data.config.fastEye = data.config.fastEye(data.conditionInfo.joyUsedIndex,:);
    end

% calibrate response (no model fitting is done otherwise)
[err,predModel,respJoy,respJoyCalib,t,stim,n] = b_s.getErr(modelfit, data);

% vectors 
contrastLE = reshape([data.experiment.LEcontrast],1,[])';
contrastRE = reshape([data.experiment.REcontrast],1,[])';
response = reshape(cell2mat(respJoyCalib'),1,[])';

% unclear to me if it's a good idea to zero-center, skip for now
vars = table(contrastLE, contrastRE, contrastLE.*contrastRE, response);
mod = fitlm(vars, 'linear', 'RobustOpts', 'on', 'Intercept', true);



figure(2); clf; hold on;
scatter3(contrastRE,contrastLE, response, 'filled', 'ColorVariable', 'response');
zlabel('response'); ylabel('LE contrast'); xlabel('RE contrast')