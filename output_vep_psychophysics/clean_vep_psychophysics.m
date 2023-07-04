% clean_motor.m
%
% We need to do 2 cleaning operations: (1) determine that a participant
% did/did not use the joystick as they were instructed (i.e. confirming
% they followed directions), and remove any "flat" trials where they didn't
% use the joystick as this will increase the error in fit. and (2) very
% rarely there is at least 1 case where a participant might have 2 files
% (i.e. started, stopped, and then started a new run). Need to combine
% those.
%
% This script takes data from ../output_vep_psychophyics/
% output_vep_psychophysics_uncleaned, processes it, then saves it in
% ../output_vep_psychophyics/.
%
% For each participant there will be:
%  Data files
%  - subjectID-motor-congruent.mat
%  - subjectID-motor-orthogonal.mat
%
%  Load these for a structure called congruentMotor or orthogonalMotor
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
%
% The two figures are for inspection to make sure that this script has
% successfully classified peoples' joystick trials.
%
% These files contain two plots for quality control inspection. Top plot is
% all trials flagged as joystick used, and bottom is all trials flagged as
% joystick not used. Look at these plots to make sure that the trials were
% properly sorted, and adjust the script if necessary. Most commonly there
% is some weird spike in the flat-line data that causes the script to
% classify this trial as "joystick used" when it really wasn't. See special
% cases below for dealing with them.

clear; close all;
plotOn = 1;
skipCong = 0; % flag for skipping this condition
skipOrth = 0;

dataDir = [cd filesep 'output_vep_psychophysics_uncleaned'];

%% Load subject IDs

% load participants we're fitting:
%run([cd '/../subjectList_vep.m']); % puts variable called sID in workspace
% or individual subjects for when new people added:
% sID = {'AM_LE_BA_15','AM_LE_LF_16'};
sID = {'AM_RE_AK_25'}

for i = 1:length(sID)

    disp(['Working on ' sID{i} ' (' num2str(i) ' of ' num2str(length(sID)) ')'])

    %% load - congruent condition
    if ~skipCong
        disp(' .. processing congruent run')
        % notes
        % AM_LE_LF_16: trial 43 should have, didn't (confirmed)
        % AM_LE_QK_24: trial 2, 4 didn't need to, did (confirmed)
        % AM_RE_KC_24: trial 17 should have, didn't (confirmed); trial 16,18
        %              didnt need to but did (confirmed); weird electrical
        %              interference (?) on trial 28 (confirmed); needs crit=1
        % AM_RE_XV_19: trial 34 has single-sample blips (confirmed)
        % BD_LE_HU_63: trial 2 didnt need to, did (confirmed)
        % NS_YJ_28: trial 12,44 didn't need to, did (confirmed)

        file.congru = dir([dataDir filesep sID{i} '*_EEG_*.mat']);
        switch size(file.congru,1)
            case 0 % shouldn't occur since everyone does the congruent first
                disp('!!!!!! No congruent condition detected, this ought not to happen !!!!!!!')
                return;

            case 1 % one file, nothing went wrong
                data = load([file.congru.folder filesep file.congru.name]);

            otherwise
                disp('More than 1 file detected, code something to concatenate runs')
                return;
        end

        data.conditionInfo.condition = 'congruent'; % init field - congruent condition

        % this is where they were told to use (or not) the joystick
        joyOnIndex = logical(data.experiment.joystickReportOn');

        crit = 2;
        % using > 2 instead of ~= 0 removes trials with a tiny bump/waver
        % (validated by checking a few participants to confirm)

        % special case: 2 is not a good crit for AM_RE_KC_24, flags incorrectly
        if strcmpi(sID{i},'AM_RE_KC_24')
            crit = 1;
        end
        joyUsedIndex = sum(abs(diff([data.experiment.binoResponseStart ...
            data.experiment.response data.experiment.binoResponseEnd],1,2)),2) > crit;
        % joyUsedIndex: whether they actually used (or not) the joystick, such
        % that 1 = used it and 0 = didn't use it
        trialID = (1:length(data.experiment.joystickReportOn))';


        % % % Dealing with special cases % % %
        %
        if strcmpi(sID{i},'AM_RE_XV_19')
            % - AM_RE_XV_19's trial 34 has this very strange blip
            %   must be electrial (joystick isn't moving but has single sample
            %   that is wildly different)
            t = 34;
            data.experiment.binoResponseStart(t, 74) = mean(data.experiment.binoResponseStart(t,[73 75]));
            data.experiment.response(t, 675) = mean(data.experiment.response(t,[674 676]));
            joyUsedIndex(t) = 0;
            clear t
            %
        elseif strcmpi(sID{i},'AM_RE_KC_24')
            % - AM_RE_KC_24's trial 28 has oscillation (??) between 0.18
            %   0.19, I suspect it was just positioned oddly - this flags
            %   it as a trial on which the lever was moved, but it wasn't
            t = 28;
            repl = mean([data.experiment.binoResponseStart(t,:) data.experiment.response(t,:) data.experiment.binoResponseEnd(t,:)]);
            data.experiment.binoResponseStart(t, :) = repl;
            data.experiment.response(t, :) = repl;
            data.experiment.binoResponseEnd(t, :) = repl;
            joyUsedIndex(t) = 0;
            clear t repl
        end

        % grab the trials where they should have used it and check that they
        % did (i.e., indexed trials should be 1):
        if ~all(joyOnIndex==joyUsedIndex)
            % make a note
            str = {'Participant did not obey all condition instructions, use variable joyUsedIndex==1 to determine when they were moving the joystick'};
            if ~isfield(data.conditionInfo, 'notes')
                data.conditionInfo.notes = str;
            else
                data.conditionInfo.notes{size(data.conditionInfo.notes, 1)+1, :} = str;
            end
        else
            str = {'Participant obeyed all condition instructions, joyOnIndex will == joyUsedIndex'};
            if ~isfield(data.conditionInfo, 'notes')
                data.conditionInfo.notes = str;
            else
                data.conditionInfo.notes{size(data.conditionInfo.notes, 1)+1, :} = str;
            end
        end

        if any(joyOnIndex(~joyUsedIndex))
            % if true, the subject didn't move the joystick on a trial where
            % they should have
            fixTrials = trialID(joyOnIndex & joyOnIndex ~= joyUsedIndex);

            % make a note
            str = {['Participant didn''t use joystick as instructed for trials: ' num2str(fixTrials')]};
            if ~isfield(data.conditionInfo, 'notes')
                data.conditionInfo.notes = str;
            else
                data.conditionInfo.notes{size(data.conditionInfo.notes, 1)+1, :} = str;
            end
            clear fixTrials
        end

        if any(~joyOnIndex(joyUsedIndex))
            % if true, the subject didn't move the joystick on a trial where
            % they should have
            fixTrials = trialID(~joyOnIndex & joyOnIndex ~= joyUsedIndex);

            % make a note
            str = {['Participant instructed not to use the joystick but did, for trials: ' num2str(fixTrials')]};
            if ~isfield(data.conditionInfo, 'notes')
                data.conditionInfo.notes = str;
            else
                data.conditionInfo.notes{size(data.conditionInfo.notes, 1)+1, :} = str;
            end
            clear fixTrials
        end

        % Save all of the condition info
        data.experiment.joystickActuallyUsed = joyUsedIndex;
        data.conditionInfo.joyOnIndex = joyOnIndex; % repetitive but useful later
        data.conditionInfo.joyUsedIndex = joyUsedIndex;


        % quick visual plot - check this to confirm trials were sorted
        % correctly
        if plotOn==1
            fig1=figure(1); clf;

            subplot(2,1,1); % show trials where joystick was used
            t = trialID(data.conditionInfo.joyUsedIndex);
            p = plot([data.experiment.binoResponseStart(t,:) data.experiment.response(t,:) data.experiment.binoResponseEnd(t,:)]');
            title([sID{i} ', congruent: Motor response on trials where joystick marked as "used"'])
            legend(p, num2str(t), 'Location', 'EastOutside');
            subplot(2,1,2); % show trials where joystick was not used
            t = trialID(~data.conditionInfo.joyUsedIndex);
            p = plot([data.experiment.binoResponseStart(t,:) data.experiment.response(t,:) data.experiment.binoResponseEnd(t,:)]');
            title([sID{i} ', congruent: Trials where joystick marked as "do not use"'])
            legend(p, num2str(t), 'Location', 'EastOutside');
            set(fig1, 'Position', [0 0 800 800]);
            input('Press enter to close plot and continue', 's')
           close(fig1);
        end

        clear crit file joyOnIndex joyUsedIndex str trialID

        % Remove fields unnecessary for analysis to keep things clean
        data = rmfield(data, {'bigbox', 'gammaTable', 'holdScreenSec', 'keyCode', ...
            'keyIsDown', 'keybIndices', 'keybNames', 'keyboardnum', 'nextKey', ...
            'response', 'responseBinoEnd', 'responseBinoStart', 'runNum', ...
            'runsComplete', 'screen_clut', 'seconds', 'smallbox', 'theScreenArea', ...
            'homedir', 'savedir'});

        congruentMotor = data;
        clear data


        save([cd filesep sID{i} '-motor-congruent.mat'], 'congruentMotor');


    end % skipflag

    %% load - orthogonal condition
    if ~skipOrth
        disp(' .. processing orthogonal run')
        % notes
        % AM_RE_KC_24 used the joystick on every trial! (confirmed)
        % AM_LE_QK_24, trial 51: used joystick, shouldn't have (confirmed,
        % fixed by changing crit = 1.99 instead of crit = 2)
        % AM_LE_BA_15: 48 didn't move but should have, 47 did move but
        % shouldn't have (confirmed)


        file.orthog = dir([dataDir filesep sID{i} '*_EEGorthog_*.mat']);
        switch size(file.orthog,1)
            case 0 % shouldn't occur since everyone does the congruent first
                disp('No orthogonal condition detected, skipping')
                continue;

            case 1 % one file, nothing went wrong
                data = load([file.orthog.folder filesep file.orthog.name]);

            otherwise
                disp([num2str(size(file.orthog,1)) ' files detected, concatenating runs'])

                data = load([file.orthog(1).folder filesep file.orthog(1).name]);

                % flag this because the config may indicate a different
                % number/order of trials than was actually run
                data.conditionInfo.notes = {'Multiple runs were concatenated in this file - the config won''t be strictly correct; multiple offsetL/R values reflect each concatenated run'};

                data.experiment.joystickReportOn(end) = [];
                for f=2:size(file.orthog,1)
                    tmp = load([file.orthog(f).folder filesep file.orthog(f).name]);
                    tmp.experiment.joystickReportOn(end) = [];
                    data.experiment.joystickReportOn = [data.experiment.joystickReportOn tmp.experiment.joystickReportOn];
                    data.experiment.binoResponseStart = [data.experiment.binoResponseStart; tmp.experiment.binoResponseStart];
                    data.experiment.binoResponseEnd = [data.experiment.binoResponseEnd; tmp.experiment.binoResponseEnd];
                    data.experiment.BEcontrastStart = [data.experiment.BEcontrastStart; tmp.experiment.BEcontrastStart];
                    data.experiment.BEcontrastEnd = [data.experiment.BEcontrastEnd; tmp.experiment.BEcontrastEnd];
                    data.experiment.response = [data.experiment.response; tmp.experiment.response];
                    data.experiment.LEcontrast = [data.experiment.LEcontrast; tmp.experiment.LEcontrast];
                    data.experiment.REcontrast = [data.experiment.REcontrast; tmp.experiment.REcontrast];
                    data.offsetL = [data.offsetL; tmp.offsetL];
                    data.offsetR = [data.offsetR; tmp.offsetR];
                    clear tmp
                end
        end


        data.conditionInfo.condition = 'orthogonal'; % init field - 

        % this is where they were told to use (or not) the joystick
        joyOnIndex = logical(data.experiment.joystickReportOn');

        crit = 1.99;


        joyUsedIndex = sum(abs(diff([data.experiment.binoResponseStart ...
            data.experiment.response data.experiment.binoResponseEnd],1,2)),2) > crit;
        % joyUsedIndex: whether they actually used (or not) the joystick, such
        % that 1 = used it and 0 = didn't use it
        trialID = (1:length(data.experiment.joystickReportOn))';


        % grab the trials where they should have used it and check that they
        % did (i.e., indexed trials should be 1):
        if ~all(joyOnIndex==joyUsedIndex)
            % make a note
            str = {'Participant did not obey all condition instructions, use variable joyUsedIndex==1 to determine when they were moving the joystick'};
            if ~isfield(data.conditionInfo, 'notes')
                data.conditionInfo.notes = str;
            else
                data.conditionInfo.notes{size(data.conditionInfo.notes, 1)+1, :} = str;
            end
        else
            str = {'Participant obeyed all condition instructions, joyOnIndex will == joyUsedIndex'};
            if ~isfield(data.conditionInfo, 'notes')
                data.conditionInfo.notes = str;
            else
                data.conditionInfo.notes{size(data.conditionInfo.notes, 1)+1, :} = str;
            end
        end


        if any(joyOnIndex(~joyUsedIndex))
            % if true, the subject didn't move the joystick on a trial where
            % they should have
            fixTrials = trialID(joyOnIndex & joyOnIndex ~= joyUsedIndex);

            % make a note
            str = {['Participant didn''t use joystick as instructed for trials: ' num2str(fixTrials')]};
            if ~isfield(data.conditionInfo, 'notes')
                data.conditionInfo.notes = str;
            else
                data.conditionInfo.notes{size(data.conditionInfo.notes, 1)+1, :} = str;
            end
            clear fixTrials
        end

        if any(~joyOnIndex(joyUsedIndex))
            % if true, the subject did move when they shouldn't have
            fixTrials = trialID(~joyOnIndex & joyOnIndex ~= joyUsedIndex);

            % make a note
            str = {['Participant instructed not to use the joystick but did, for trials: ' num2str(fixTrials')]};
            if ~isfield(data.conditionInfo, 'notes')
                data.conditionInfo.notes = str;
            else
                data.conditionInfo.notes{size(data.conditionInfo.notes, 1)+1, :} = str;
            end
            clear fixTrials
        end

        % Save all of the condition info
        data.experiment.joystickActuallyUsed = joyUsedIndex;
        data.conditionInfo.joyOnIndex = joyOnIndex; % repetitive but useful later
        data.conditionInfo.joyUsedIndex = joyUsedIndex;


        % quick visual plot - check this to confirm trials were sorted
        % correctly
        if plotOn ==1
            fig1=figure(1); clf;

            subplot(2,1,1); % show trials where joystick was used
            t = trialID(data.conditionInfo.joyUsedIndex);
            p = plot([data.experiment.binoResponseStart(t,:) data.experiment.response(t,:) data.experiment.binoResponseEnd(t,:)]');
            title([sID{i} ', orthogonal: Motor response on trials where joystick marked as "used"'])
            legend(p, num2str(t), 'Location', 'EastOutside');
            subplot(2,1,2); % show trials where joystick was not used
            t = trialID(~data.conditionInfo.joyUsedIndex);
            p = plot([data.experiment.binoResponseStart(t,:) data.experiment.response(t,:) data.experiment.binoResponseEnd(t,:)]');
            title([sID{i} ', orthogonal: Trials where joystick marked as "do not use"'])
            legend(p, num2str(t), 'Location', 'EastOutside');
            set(fig1, 'Position', [0 0 800 800]);
            input('Press enter to close plot and continue', 's')
            close(fig1);
        end
        clear crit file joyOnIndex joyUsedIndex str trialID

        % Remove fields unnecessary for analysis to keep things clean
        data = rmfield(data, {'bigbox', 'gammaTable', 'holdScreenSec', 'keyCode', ...
            'keyIsDown', 'keybIndices', 'keybNames', 'keyboardnum', 'nextKey', ...
            'response', 'responseBinoEnd', 'responseBinoStart', 'runNum', ...
            'runsComplete', 'screen_clut', 'seconds', 'smallbox', 'theScreenArea', ...
            'homedir', 'savedir'});

        orthogonalMotor = data;
        clear data


        save([cd filesep sID{i} '-motor-orthogonal.mat'], 'orthogonalMotor');

    end
end

clear; close all;