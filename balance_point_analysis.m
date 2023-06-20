% balance_point_analysis
%
% compares the estimate of balance point from our model to balance point
% estimates made previously in the literature
%
% written IF KM 07/03/2023

close all; clear all
MIDRANGEFLAG = 0; % if using model fits based on small discrepancies between left and right eye
% load participants we're fitting:
subjectList_psychophysics; % puts variable called sID in workspace

group(1).bp = []; group(3).bp = []; group(2).bp = [];
ct = [0 0 0 ];
for i = 1:length(sID)
    clear tmp p contrast;
    if ~MIDRANGEFLAG
        load(['fitdata_psychophysics/', sID{i} ]);
        tmp.p = p; clear p
    else
        tmp =   load(['fitdata_psychophysics/midrange_fits/', sID{i}, '-midrange.mat' ], 'p');
    end
    p = tmp.p;

    if p.k(2)<p.k(1) % swop the eyes so amblyopic left
        p.k = p.k([2, 1]);
        u = p.U;
        p.U(2) = u(3);
        p.U(3) = u(2);
    end

    if contains(p.sID(1:4) ,'AM');          clr = 'r'; g = 1; offset = -.2; group(g).name = 'amblyope';
    elseif contains(p.sID(1:4) ,'BD');  clr = 'g'; g = 2;     offset  = 0 ;group(g).name  = 'strab';
    elseif contains(p.sID(1:4) ,'NS');      clr = 'b'; offset = 0.2; g = 3; group(g).name = 'control';
    end
    ct(g) = ct(g) + 1;
    % find contrast for which perceived contrast is equal across the two
    % eyes
    contrast(:, 1) = 0:.001:1;
    contrast(:, 2) = 1-contrast(:, 1);
    out = fit_model(p, contrast);
    tmp = abs([out(:, 1)-out(:, 2)]);
    ind = find(tmp==min(tmp));
    group(g).bp(ct(g), :) =contrast(ind,:);

    % find the point where the phase is intermediate, as a function of
    % contrast in the amblyopic eye
    clear contrast
    contrastAE = [ .075 1.5 3 6 12 24 48 96 ]./100; % contrast
    for e  = 1:length(contrastAE)
        contrast(:, 1) = contrastAE(e).* ones(1000, 1);
        contrast(:, 2) = linspace(0, contrastAE(e), 1000); % fellow eye is going to be same or lesser contrast

        out = fit_model(p, contrast);
        tmp = abs([out(:, 1)-out(:, 2)]);
        ind = find(tmp==min(tmp));
        group(g).ding(e).contrastAE(ct(g))  = contrastAE(e);
        group(g).ding(e).contrastFE(ct(g))   = contrast(ind, 2);
    end

    for expflag = 1:2
        if  expflag == 1
            % Beylerian et al. 2020. Interocular suppressive interactions in amblyopia depend on spatial frequency
            % Mid spatial frequency 1.31 cpd
            % mask set at 5x it's detection threshold (supplementary table 1)
            group(1).thresh(expflag).c_mon_thr_AE  = 1/70;   % monocular thresholds
            group(1).thresh(expflag).c_masked_AE = 1/17; %  experimental values, contrast needed to see masked grating
            group(1).thresh(expflag).c_masker_FE = 0.05527*5; % masking contrasts
            group(3).thresh( expflag).c_mon_thr_AE =  1/130;
            group(3).thresh(expflag).c_masked_AE = 1/30;
            group(3).thresh(expflag).c_masker_FE =  0.03329*5;

        else
            % Zhou et al 2018. Amblyopic Suppression: Passive Attenuation, Enhanced Dichoptic Masking by the Fellow Eye or Reduced Dichoptic Masking by the Amblyopic Eye?
            % mask set at 5x it's detection threshold
            group(1).thresh(expflag).c_mon_thr_AE  = 1/60;  % monocular thresholds
            group(1).thresh(expflag).c_masked_AE = 1/22; %  experimental values, contrast needed to see masked grating
            group(1).thresh(expflag).c_masker_FE = (1/20)*5; % masking contrasts
            group(3).thresh(expflag).c_mon_thr_AE =  1/130;
            group(3).thresh(expflag).c_masked_AE = 1/35;
            group(3).thresh(expflag).c_masker_FE =  (1/40)*5;
        end
        if g~=2 % no data for strabs
            out = fit_model(p, [group(g).thresh(expflag).c_mon_thr_AE 0 ] );
            group(g).thresh(expflag).pcAEthr(ct(g))  = max(out); % find the perceived contrast at threshold

            % contrast needed to reach AE threshold when unmasked
            c_nomask= find_threshold(p, group(g).thresh(expflag).pcAEthr(ct(g)),  0, 1); % initialize value
            % contrast needed to reach AE threshold when masked
            c = find_threshold(p, group(g).thresh(expflag).pcAEthr(ct(g)),  group(g).thresh(expflag).c_masker_FE, 1); % initialize value

            group(g).thresh(expflag).c_masked_AE_sim(ct(g)) = c;
            group(g).thresh(expflag).c_nomask_AE_sim(ct(g)) = c_nomask;
            group(g).thresh(expflag).Sim_dB_El(ct(g))  = 20*log10( c/c_nomask);
        end
    end
end

% get the numbers
for g = 1:3
    group(g).bp_mean = mean([cat(1,group(g).bp)]);
    group(g).bp_ste = std([cat(1,group(g).bp)])/sqrt(ct(g));

    for e = 1:length(group(g).ding)
        group(g).ding(e).bp_mean = mean([cat(1,group(g).ding(e).contrastFE)]);
        group(g).ding(e).bp_ste = std([cat(1,group(g).ding(e).contrastFE)])/sqrt(ct(g));
    end
    if g~=2
        for expflag = 1:2
            group(g).thresh(expflag).c_masked_AE_sim_mean  =  mean([cat(1,group(g).thresh(expflag).c_masked_AE_sim)]);
            group(g).thresh(expflag).c_masked_AE_sim_ste =  std([cat(1,group(g).thresh(expflag).c_masked_AE_sim)])/sqrt(ct(g));

            group(g).thresh(expflag).Exp_dB_El = 20*log10( group(g).thresh(expflag).c_masked_AE/group(g).thresh(expflag).c_mon_thr_AE );
            group(g).thresh(expflag).Exp_dB_El_mean = mean([cat(1,  group(g).thresh(expflag).Sim_dB_El)]);
            group(g).thresh(expflag).Exp_dB_El_ste = std([cat(1,  group(g).thresh(expflag).Sim_dB_El)])./sqrt(ct(g));
        end
    end
end


function out = fit_model(p, in)
% returns the perceived contrast in left and right eye as a function of the
% contrasts in each eye
[out, p] = b_s.linear_attenuation(p, in);
[out, p] = b_s.normalization(p, out);
end

function c = find_threshold(p, per_contrast, contrast_other_eye, eye)
contrast = ones(5000, 2);
contrast(:, eye) = linspace(0, 1, 5000);
contrast(:, mod(eye, 2)+1) = contrast_other_eye.*ones(1, 5000);
out = fit_model(p, contrast);
tmp = abs([out(:, eye)-per_contrast]);
ind = find(tmp==min(tmp));
c = contrast( ind, eye);
end

