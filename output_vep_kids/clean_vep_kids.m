% clean_vep_kids.m
%
% This script was adapted from dynamic_contrast/output_vep/clean_vep.m to
% work on the newly collected kids data, so if some comments seem wrong/out
% of date that is why. The difference is that kids were tasked with a
% detection task unrelated to the contrast, rather than the contrast-
% estimation task. 
%
% This script takes data (based on subjectID) from
% ../output_vep_kids/output_vep_kids_uncleaned, and saves information
% to the struct which gets saved in ../output_vep_kids/.
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
dataDir = [cd filesep 'output_vep_kids_uncleaned'];

%% go thru subjects

sID = {'NS_DXK_L_12'}
% sID = {...
% 'NS_APB_G_13',...
% 'NS_BXK_L_09',...
% 'NS_DXK_L_12',...
% 'NS_DXN_G_10',...
% 'NS_ERQ_G_10',...
% 'NS_HKN_L_10',...
% 'NS_IKJ_L_12',...
% 'NS_MKB_G_13',...
% 'NS_MQZ_G_11',...
% 'NS_MYB_G_10',...
% 'NS_QIU_G_09',...
% 'NS_YJP_L_15'};


for i = 1:length(sID)

    disp(['Working on ' sID{i} ' (' num2str(i) ' of ' num2str(length(sID)) ')'])

    %% load
    runfile = dir([dataDir filesep sID{i} '*EEG*.mat']);    % contains data saved on linux computer
    vepfile = dir([dataDir filesep sID{i} '*hapvepv.mat']); % contains eeg data

    if any([isempty(runfile) isempty(vepfile)])
    disp('File load error, there should be a *hapvepv.mat & *_EEG_*.mat file for participant')
    return;
    end

    disp(' .. processing run')
    rundata = load([runfile.folder filesep runfile.name]);
    vepdata = load([vepfile.folder filesep vepfile.name]);

    % initialize structure
    dataVep = struct();

    % fill fields that some of the other analysis scripts expect
    dataVep.binocular = rundata.binocular;
    dataVep.calculatedFrameRate = rundata.calculatedFrameRate;
    dataVep.carrierFreq = rundata.carrierFreq;
    dataVep.config = rundata.config;
    dataVep.display = rundata.display;
    dataVep.experiment = rundata.experiment; % contains presented contrast
    dataVep.hz = rundata.hz;
    dataVep.ifi = rundata.ifi;
    dataVep.offsetL = rundata.offsetL;
    dataVep.offsetR = rundata.offsetR;
    dataVep.openedScreenRes = rundata.openedScreenRes;
    dataVep.sID = rundata.sID;
    dataVep.stim = rundata.stim;
    dataVep.t = rundata.t;

    % initialize response vectors for VEP data
    dataVep.experiment.binoResponseStart = [];
    dataVep.experiment.response = [];
    dataVep.experiment.binoResponseEnd = [];

    % rescale everything 0-1 (trial-by-trial)
    A = [vepdata.experiment.binoResponseStart vepdata.experiment.response vepdata.experiment.binoResponseEnd];
    rowmin = min(A,[],2);
    rowmax = max(A,[],2);
    A_scaled = rescale(A,'InputMin',rowmin,'InputMax',rowmax);

        tempStart = A_scaled(:,1:size(vepdata.experiment.binoResponseStart,2));
        tempResponse = A_scaled(:, size(tempStart,2)+1 : (size(tempStart,2)) + size(vepdata.experiment.response,2));
        tempEnd = A_scaled(:, size(tempStart, 2)+size(tempResponse, 2)+1 : end);


    dataVep.experiment.response = tempResponse;
    dataVep.experiment.binoResponseStart = tempStart;
    dataVep.experiment.binoResponseEnd = tempEnd;
    dataVep.conditionInfo.notes = {'experiment.response contains VEP data'};


    save([cd filesep sID{i} '-vep-data.mat'], 'dataVep');


end

clear; close all;