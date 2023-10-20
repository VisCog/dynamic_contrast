% rename vep files
clear
% Some notes on renamed files, was > now:
% - AM_LE_AM_16 > BD_LE_AM_16 (strab with equal vision, not amblyopia)
% - AM_RE_BA_15 > AM_LE_BA_15 (reported incorrect eye on screening)
% - AM_LE_HU_63 > BD_LE_HU_63 (strab with equal vision, not amblyopia)
% - AM_RE_L_IA_22 > orthogonal run, renamed AM_RE_IA_24 and merged
% - BD_YT_63 > NS_YT_63 (screening said possible childhood disorder, KTH
%                       confirmed can be a control)
% - AM_RE_CK_23 > BD_RE_AK_25 (initials and age typo'd at first entry, plus
%                             has strab with equal vision not amblyopia)
% - AM_RE_XV_19 > NS_XV_19 (had teenage vision problems that don't turn out
%                          to be amblyopia, likely just needed glasses)
% - AM_LE_HH_35 > BD_LE_HH_35 (strab with equal vision, not amblyopia)



old = {'AM_RE_BA_15', 'AM_RE_L_IA_22', 'AM_LE_HH_35', 'AM_RE_CK_23', 'AM_RE_XV_19', 'AM_LE_AM_16', 'AM_LE_HU_63', 'BD_YT_63'}
new = {'AM_LE_BA_15', 'AM_RE_IA_24',   'BD_LE_HH_35', 'BD_RE_AK_25', 'NS_XV_19',    'BD_LE_AM_16', 'BD_LE_HU_63', 'NS_YT_63'}





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





% rename hapvepv2.mat > hapvepv.mat
clear
f = dir(['/Users/bear/Documents/UW/allthegit/dynamic_contrast/output_vep/output_vep_uncleaned/*hapvepv2.mat']);
for i = 1:length(f)
    rnThis = [f(i).folder filesep f(i).name];
    rnTo = strrep(rnThis, 'hapvepv2.mat', 'hapvepv.mat');
    movefile(rnThis, rnTo);
end