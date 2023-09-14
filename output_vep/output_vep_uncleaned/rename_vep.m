% rename vep files
clear
% Some notes on renamed files, was > now:
% - AM_LE_AM_16 > BD_LE_AM_16
% - AM_RE_BA_15 > AM_LE_BA_15
% - AM_LE_HU_63 > BD_LE_HU_63
% - AM_RE_L_IA_22 > had an orthogonal run, renamed AM_RE_IA_24 and merged
%   into the folder with that name for the original run
% - BD_YT_63 > NS_YT_63

%old = {'AM_LE_AM_16'}
%new = {'BD_LE_AM_16'}

%old = {'AM_RE_BA_15'}
%new = {'AM_LE_BA_15'}

%old = {'AM_LE_HU_63'}
%new = {'BD_LE_HU_63'}

%old = {'AM_RE_L_IA_22'}
%new = {'AM_RE_IA_24'}

% old = {'BD_YT_63'}
% new = {'NS_YT_63'}

%old = {'AM_RE_CK_23'}
%new = {'AM_RE_AK_25'}

old = {'AM_LE_HH_35', 'AM_RE_AK_25', 'AM_RE_XV_19'}

new = {'BD_LE_HH_35', 'BD_RE_AK_25', 'NS_XV_19'}

% being very specific so no one can accidentally run this 
fileDir = '/Users/bear/Documents/UW/allthegit/dynamic_contrast/output_vep/output_vep_uncleaned/';

for i = 1:size(old,2)
    fileMove = dir([fileDir old{i} '*.mat']);

    disp(['Found ' num2str(size(fileMove,1))  ' files, renaming ' old{i} ' to ' new{i}])

    for f = 1:size(fileMove,1)
        tmp = strrep([fileMove(f).folder filesep fileMove(f).name], old{i}, new{i});
        movefile([fileMove(f).folder filesep fileMove(f).name], tmp);

    end
end