sID = {...
'AM_LE_BA_15', ...
'AM_LE_LF_16', ...
'AM_LE_QK_24', ...
'AM_LE_QQ_37', ...
'AM_RE_CI_20', ...
'AM_RE_IA_24', ...
'AM_RE_ID_16', ...
'AM_RE_KC_24', ...
'AM_RE_LK_25', ...
'AM_RE_WH_18', ...
'AM_RE_YA_26', ...
'BD_LE_AM_16', ...
'BD_LE_HH_35', ...
'BD_LE_HU_63', ...
'BD_RE_AK_25', ...
'NS_AA_17', ...
'NS_AC_48', ...
'NS_CW_25', ...
'NS_EE_21', ...
'NS_HJ_24', ...
'NS_JK_17', ...
'NS_JX_19', ...
'NS_KQ_17', ...
'NS_PB_25', ...
'NS_PP_17', ...
'NS_RA_37', ...
'NS_TN_23', ...
'NS_XV_19', ...
'NS_YJ_23', ...
'NS_YJ_28', ...
'NS_YN_17', ...
'NS_YT_63', ...
};

% %%
% % % %% to update sID list based on available data
% clear
% idlist = dir([cd filesep 'output_vep' filesep 'output_vep_uncleaned' filesep '*hapvepv.mat']);
% idx = contains({idlist.name}, 'AM_') ...  rm the example toy datasets
%     | contains({idlist.name}, 'BD_') ...
%     | contains({idlist.name}, 'NS_'); 
% idlist = idlist(idx);
% idlist = {idlist.name};
% idlist = strrep(idlist, '_orthog_hapvepv.mat', '');
% idlist = strrep(idlist, '_hapvepv.mat', '');
% idlist = unique(idlist);
% %print a thing to copy-and-paste to replace the above
% for i = 1:length(idlist)
%     disp(['''' idlist{i} ''', ...'])
% end
