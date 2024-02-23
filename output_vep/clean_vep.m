% clean_vep.m
%
% This script takes data (based on subjectID) from
% ../output_vep/output_vep_uncleaned, looks at the
% corresponding data in output-motor-cleaned to determine which trials
% participants used or didn't use the joystick (.joystickActuallyUsed), and
% saves that information to the struct which gets saved in
% output-vep-cleaned.
%
% Additionally, the data are scaled such that EACH TRIAL's min and max are
% 0 and 1, respectively.
%
% In each of the participant folders in there, there will be:
%  Data files
%  - subjectID-vep-congruent.mat
%  - subjectID-vep-orthogonal.mat
%
%  Load these for a structure called congruentVep or orthogonalVep
%  (everyone will have congruent data; people who came back for a second
%  EEG visit will have orthogonal data). The structs have the same fields.
%
%  Some important fields in these structs:
%
%     .experiment
%                .joystickReportOn = participant was told to use (1) or not
%                 use (0) the joystick to make a response (regardless of
%                 whether they actually did or not)
%                .joystickActuallyUsed = participant used (1) or did not
%                use (0) the joystick during that trial
%                .BEcontrastStart, .BEcontrastEnd = the contrast presented
%                 to both eyes during the first 3.5 and the last 3.5
%                 seconds of the trial
%                .binoResponseStart, .binoResponseEnd = the response made
%                 to the contrast presented during the start/end of trial
%                .REcontrast, .LEcontrast = the contrast presented to the
%                 left and right eyes during a trial
%                .response = the motor (joystick) response during the trial
%
%    .conditionInfo
%                  .condition = 'congruent' or 'orthogonal'
%                  .joyOnIndex = the same as experiment.joystickReportOn
%                  .joyUsedIndex = the same information as
%                   experiment.joystickActuallyUsed; use this to select
%                   which trials (==true) to send to the fit function
%                  .notes = strings of notes regarding where
%                   joyUsedIndex ~= joyOnIndex


clear; close all;
dataDir = [cd filesep 'output_vep_uncleaned'];

%% go thru subjects

% load participants we're fitting:
%run([cd '/../subjectList_vep.m']); % puts variable called sID in workspace

%sID = {'NS_HJ_24', 'AM_RE_KR_19', 'AM_LE_RT_26'}
sID = {'NS_HJ_24'}

allChan = 0; % flag to load the all-channels 2ar file instead of Oz-only


for i = 1:length(sID)

    disp(['Working on ' sID{i} ' (' num2str(i) ' of ' num2str(length(sID)) ')'])

    % reset flag
    existOrthog = 0;

    %% load

    if allChan == 1
        files = dir([dataDir filesep sID{i} '*hapvepv2ar.mat']);
    else
        files = dir([dataDir filesep sID{i} '*hapvepv.mat']);
    end

    switch size(files,1)

        case 0 % shouldn't occur
            disp('No files detected, this ought not to happen')
            return;

        case 1 % one file, confirm it's congruent (not orthog)
            existOrthog = 0;

            if isempty(strfind(files.name, 'orthog'))
                % all good
                data = load([files.folder filesep files.name]);
            elseif ~isempty(strfind(files.name, 'orthog'))
                disp('This person has an orthogonal file but not congruent - should not happen')
                return;
            end

        case 2 % two files, one should be congruent and one should be orthog
            existOrthog = 1;

            % for alpha reasons, file 1 should be cong and file 2 orthog
            % check this for now and fix later if turns out to be untrue..
            if isempty(strfind(files(1).name, 'orthog'))
                % all good
                data = load([files(1).folder filesep files(1).name]);
            elseif ~isempty(strfind(files(1).name, 'orthog'))
                disp('unexpectedly orthogonal file, look into this')
                return;
            end

            % quick check to make sure orthog is in second position
            if isempty(strfind(files(2).name, 'orthog'))
                disp('unexpectedly not an orthogonal file, look into this')
                return;
            end

        otherwise
            disp('More than 2 files detected, should not happen')
            return;
    end

    % start with congruent
    disp(' .. processing congruent run')
    try
        motorData = load([dataDir filesep '..' filesep '..' filesep 'output_vep_psychophysics' ...
            filesep sID{i} '-motor-congruent.mat']);
    catch ME
        if (strcmp(ME.identifier,'MATLAB:load:couldNotReadFile'))
            msg = ['output_vep_psychophysics' filesep sID{i} '-motor-congruent.mat' newline ...
                'does not exist. You need to process the motor data before' newline ...
                'processing the VEP data so that the trials can be correctly' newline ...
                'sorted into "joystick used/not used" trials. ' newline 'Did you do that?'];
            causeException = MException('MATLAB:myCode:couldNotReadFile',msg);
            ME = addCause(ME,causeException);
        end
        rethrow(ME)
    end

    % copy over motor file info (contains a bit extra info) and replace the
    % joystick position with the VEP
    congruentVep = motorData.congruentMotor;

    % clear the motor data to avoid duplication
    congruentVep.experiment.binoResponseStart = [];
    congruentVep.experiment.response = [];
    congruentVep.experiment.binoResponseEnd = [];

    % rescale everything 0-1 (trial-by-trial)
    A = [data.experiment.binoResponseStart data.experiment.response data.experiment.binoResponseEnd];
    rowmin = min(A,[],2);
    rowmax = max(A,[],2);
    A_scaled = rescale(A,'InputMin',rowmin,'InputMax',rowmax);

    if allChan == 1 % 3d
        tempStart = A_scaled(:,1:size(data.experiment.binoResponseStart,2),:);
        tempResponse = A_scaled(:, size(tempStart,2)+1 : (size(tempStart,2)) + size(data.experiment.response,2),:);
        tempEnd = A_scaled(:, size(tempStart, 2)+size(tempResponse, 2)+1 : end,:);
    else
        tempStart = A_scaled(:,1:size(data.experiment.binoResponseStart,2));
        tempResponse = A_scaled(:, size(tempStart,2)+1 : (size(tempStart,2)) + size(data.experiment.response,2));
        tempEnd = A_scaled(:, size(tempStart, 2)+size(tempResponse, 2)+1 : end);
    end

    congruentVep.experiment.response = tempResponse;
    congruentVep.experiment.binoResponseStart = tempStart;
    congruentVep.experiment.binoResponseEnd = tempEnd;

    % make a note
    if allChan == 1
        congruentVep.experiment.chanlocs = data.chanlocs;
        str = {'experiment.response contains VEP data, all channels, average reference'};
    else
        str = {'experiment.response contains VEP data'};
    end

    if ~isfield(congruentVep.conditionInfo, 'notes')
        congruentVep.conditionInfo.notes = str;
    else
        congruentVep.conditionInfo.notes{size(congruentVep.conditionInfo.notes, 1)+1, :} = str;
    end

    clear data motorData

    if allChan == 1
        save([cd filesep sID{i} '-vep-congruent-allchan.mat'], 'congruentVep');
    else
        save([cd filesep sID{i} '-vep-congruent.mat'], 'congruentVep');
    end

    % orthogonal, if it exists
    if allChan ~= 1 % skip orthog for now
        if existOrthog == 1
            disp(' .. processing orthogonal run')
            data = load([files(2).folder filesep files(2).name]);

            motorData = load([dataDir filesep '..' filesep '..' filesep 'output_vep_psychophysics' ...
                filesep sID{i} '-motor-orthogonal.mat']);
            orthogonalVep = motorData.orthogonalMotor;

            % rescale everything 0-1 (trial-by-trial)
            A = [data.experiment.binoResponseStart data.experiment.response data.experiment.binoResponseEnd];
            rowmin = min(A,[],2);
            rowmax = max(A,[],2);
            A_scaled = rescale(A,'InputMin',rowmin,'InputMax',rowmax);
            tempStart = A_scaled(:,1:size(data.experiment.binoResponseStart,2));
            tempResponse = A_scaled(:, size(tempStart,2)+1 : (size(tempStart,2)) + size(data.experiment.response,2));
            tempEnd = A_scaled(:, size(tempStart, 2)+size(tempResponse, 2)+1 : end);


            orthogonalVep.experiment.response = tempResponse;
            orthogonalVep.experiment.binoResponseStart = tempStart;
            orthogonalVep.experiment.binoResponseEnd = tempEnd;

            % make a note
            str = {'experiment.response contains VEP data'};
            if ~isfield(orthogonalVep.conditionInfo, 'notes')
                orthogonalVep.conditionInfo.notes = str;
            else
                orthogonalVep.conditionInfo.notes{size(orthogonalVep.conditionInfo.notes, 1)+1, :} = str;
            end

            clear data motorData

            save([cd filesep sID{i} '-vep-orthogonal.mat'], 'orthogonalVep');

        else
            disp(' .. no orthogonal run')

        end
    end

end

clear; close all;