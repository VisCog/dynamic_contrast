clear
close all

sCondition = 'congruent';

%%% load files
vepFolder = [cd filesep 'fitdata_vep' filesep 'model_fits'];     % vep
psyFolder = [cd filesep 'fitdata_vep_psychophysics' filesep 'model_fits']; %psychophysics

vepFiles = dir([vepFolder filesep '*' sCondition '.mat']);
psyFiles = dir([psyFolder filesep '*' sCondition '.mat']);

% grouping variable
idx_psy_ns = find(cellfun(@(x) strcmpi(x(1:2),'NS'), {psyFiles.name}));
idx_psy_am = find(cellfun(@(x) strcmpi(x(1:2),'AM'), {psyFiles.name}));
idx_psy_bd = find(cellfun(@(x) strcmpi(x(1:2),'BD'), {psyFiles.name}));
% grp_psy = nan(size(psyFiles));
% grp_psy(idx_psy_am) = 1;
% grp_psy(idx_psy_bd) = 2;
% grp_psy(idx_psy_ns) = 3;

idx_vep_ns = find(cellfun(@(x) strcmpi(x(1:2),'NS'), {vepFiles.name}));
idx_vep_am = find(cellfun(@(x) strcmpi(x(1:2),'AM'), {vepFiles.name}));
idx_vep_bd = find(cellfun(@(x) strcmpi(x(1:2),'BD'), {vepFiles.name}));
% grp_vep = nan(size(vepFiles));
% grp_vep(idx_psy_am) = 1;
% grp_vep(idx_psy_bd) = 2;
% grp_vep(idx_psy_ns) = 3;

% load fits and compile
all_psy = [];
for i = 1:length(psyFiles)
    load([psyFiles(i).folder filesep psyFiles(i).name]);
    %rm the extra junk
    tmp = rmfield(p, {'abs','clean_range', 'n_good', 'startT', 'junk', 'penalizeDelay', 'costflag'...
        'p', 'tau', 'm', 'offset', 'meanResponseScaled'...
        ...'predModel_meanModel', ...
        ... 'meanModelErr', 'predModel_maxModel', 'maxModelErr', 'predModel_mnmxwghtModel', ...
        ...   'minkModelErr', 'mink_wghtd_p', 'mink_wghtd_w', 'predModel_minkwghtModel', 'minkwghtModelErr' ...
        ...'minkModelErr'...
        ...'mnmxwghtModelErr', 'predModel_minkModel','predModel_softmax'...
        'predModel_meanModel', 'predModel_meanModel_wghtd'...
        'predModel_maxModel_wghtd', 'predModel_maxModel'...
        'predModel_softmax'...
        });
    if tmp.k(1) == 1 % le fe
        tmp.kAE = tmp.k(2);
        tmp.uAE = tmp.U(3);
        tmp.uFE = tmp.U(2);
        tmp.AE = 'RE';
    elseif tmp.k(2) == 1 % re fe
        tmp.kAE = tmp.k(1);
        tmp.uAE = tmp.U(2);
        tmp.uFE = tmp.U(3);
        tmp.AE = 'LE';
    end
    if tmp.k(1) == tmp.k(2)
        disp(['-------flag ' tmp.sID])
    end
    tmp.grp = tmp.sID(1:2);
    if tmp.k_arvo(1) == 1 % le fe
        tmp.kAE_arvo = tmp.k_arvo(2);
    elseif tmp.k_arvo(2) == 1 % le fe
        tmp.kAE_arvo = tmp.k_arvo(1);
    end
    all_psy = [all_psy tmp];
    clear p tmp
end



all_vep = [];
for i = 1:length(vepFiles)
    load([vepFiles(i).folder filesep vepFiles(i).name]);
    %rm the extra junk
    tmp = rmfield(p, {'abs','clean_range', 'n_good', 'startT', 'junk', 'penalizeDelay', 'costflag'...
        'p', 'tau', 'm', 'offset', 'meanResponseScaled'...
        ...'predModel_meanModel', ...
        ... 'meanModelErr', 'predModel_maxModel', 'maxModelErr', 'predModel_mnmxwghtModel', ...
        ...   'minkModelErr', 'mink_wghtd_p', 'mink_wghtd_w', 'predModel_minkwghtModel', 'minkwghtModelErr' ...
        ...'minkModelErr'...
        ...'mnmxwghtModelErr', 'predModel_minkModel','predModel_softmax'...
        'predModel_meanModel', 'predModel_meanModel_wghtd'...
        'predModel_maxModel_wghtd', 'predModel_maxModel'...
        'predModel_softmax'...
        });
    if tmp.k(1) == 1 % le fe
        tmp.kAE = tmp.k(2);
        tmp.uAE = tmp.U(3);
        tmp.uFE = tmp.U(2);
        tmp.AE = 'RE';
    elseif tmp.k(2) == 1 % re fe
        tmp.kAE = tmp.k(1);
        tmp.uAE = tmp.U(2);
        tmp.uFE = tmp.U(3);
        tmp.AE = 'LE';
    end
    if tmp.k(1) == tmp.k(2)
        disp(['-------flag ' tmp.sID])
    end
    tmp.grp = tmp.sID(1:2);
    if tmp.k_arvo(1) == 1 % le fe
        tmp.kAE_arvo = tmp.k_arvo(2);
    elseif tmp.k_arvo(2) == 1 % le fe
        tmp.kAE_arvo = tmp.k_arvo(1);
    end
    all_vep = [all_vep tmp];
    clear p tmp
end

%originally was going to do here but easier reproduce in R .. save out
all_psy = struct2table(all_psy);
all_vep = struct2table(all_vep);

writetable(all_psy, [psyFolder filesep 'scirepfits_' sCondition '_psy.csv'])
writetable(all_vep, [vepFolder filesep 'scirepfits_' sCondition '_vep.csv'])