%% plot_vep_topo_correlation
%
% Correlates the psychophysical response with the VEP response from each
% electrode.
%
% Makes a topo plot of correlation strength
%
% Can look at binocular trial portion only, or full trial portions
%
% 
% IMPORTANT: no model calibration is done for either modality.
%            also, none of the EEG data have been rescaled into the
%            interval [0,1] - it's the raw from Mark.
%
% Note this does trial-by-trial correlation (i.e. only the VEP trials where
% they were using the joystick are used) -- the very large n means tiny
% correlations are significant -- can use simulation to determine
% meaningful cutoff

%%

clear

binoTrialOnly = 1;      % 1 = binocular sections only, 0 = full trial
plotIndividuals = 0;    % 1 = one plot for each participant, 0 = no plots
plotGroupmean = 1;      % 1 = make plot for group mean, 0 = no plot

%%

subjectList_vep;

group = 'NS';
idx = cellfun(@(x) strcmpi(x(1:2),group), sID); % select group
sID = sID(idx);

%%

all_r = []; % store all participant correlations

for i=1:length(sID)

    vepdata = load([cd filesep 'output_vep' filesep 'output_vep_uncleaned' filesep sID{i} '_hapvepv2ar.mat']);
    psydata = load([cd filesep 'output_vep_psychophysics' filesep sID{i} '-motor-congruent.mat']);

    if binoTrialOnly == 0
        titlestring = 'all dichoptic/bino/mono trial portions';

        %reshape psychophysics into vector
        psyvector = psydata.congruentMotor.experiment;
        jsu = psydata.congruentMotor.experiment.joystickActuallyUsed;
        psyvector = [psyvector.binoResponseStart(jsu,:) psyvector.response(jsu,:) psyvector.binoResponseEnd(jsu,:)];
        psyvector = reshape(psyvector', [], 1);

        %reshape vep into vector
        chandata = [vepdata.experiment.binoResponseStart(jsu,:,:) vepdata.experiment.response(jsu,:,:) vepdata.experiment.binoResponseEnd(jsu,:,:)];
        [ntr, tpts, nchan] = size(chandata);
        chandata = permute(chandata, [3, 2, 1]);
        chandata = reshape(chandata, [nchan, tpts * ntr])';

    elseif binoTrialOnly == 1 % same but for binocular sections only
                        titlestring = 'only binocular trial portions';
        psyvector = psydata.congruentMotor.experiment;
        jsu = psydata.congruentMotor.experiment.joystickActuallyUsed;
        psyvector = [psyvector.binoResponseStart(jsu,:) psyvector.binoResponseEnd(jsu,:)];
        psyvector = reshape(psyvector', [], 1);
        chandata = [vepdata.experiment.binoResponseStart(jsu,:,:) vepdata.experiment.binoResponseEnd(jsu,:,:)];
        [ntr, tpts, nchan] = size(chandata);
        chandata = permute(chandata, [3, 2, 1]);
        chandata = reshape(chandata, [nchan, tpts * ntr])';
    end

    % correlate vectors
    [r,p] = corr(psyvector, chandata, 'rows','pairwise');

    if plotIndividuals == 1

        %%%%% figure
        fig1=figure(1); clf; hold on;
        if i==1
            set(fig1, 'Position', [100 100 600 800]);
        end

        % histogram of correlation values
        subplot(4,1,1); hold on; set(gca, 'FontSize', 14);
        nbins = 20;
        h=histogram(r,nbins);
        h.FaceColor = [0.3 0.6 1];
        h.EdgeColor = 'k';
        h.LineWidth = 1.5;
        edges = h.BinEdges;
        for binIdx = 1:nbins % label electrodes
            binIndices = find(r >= edges(binIdx) & r < edges(binIdx + 1));
            binLabels = {vepdata.chanlocs(binIndices).labels};
            labelX = mean([edges(binIdx), edges(binIdx + 1)]);
            text(labelX, 0.02, strjoin(binLabels, '\n'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10);
        end
        ylim([0 max(h.BinCounts+1)]);
        ylabel('frequency')
        xlabel('correlation with psychophysical response')
        set(gca, 'YGrid', 'on');

        title([sID{i} ', ' titlestring], 'Interpreter','none')

        % electrode correlation heatmap
        subplot(4, 1, 2:4); hold on; set(gca, 'FontSize', 14);

        usecolormap = parula(64);
        usecolormap = rescale(usecolormap, 0.5, 1);
        colormap(usecolormap);

        color_index = ceil(abs(r) * (size(usecolormap, 1)));
        col = usecolormap(color_index, :);

        scatter([vepdata.chanlocs.X], [vepdata.chanlocs.Y], 500, col, 'o', 'filled', 'MarkerEdgeColor', [.9 .9 .9], 'LineWidth', 2);
        for c = 1:size(vepdata.chanlocs,2)
            text(vepdata.chanlocs(c).X, vepdata.chanlocs(c).Y, strjoin({vepdata.chanlocs(c).labels, num2str(r(c), 2)}, ': '),...
                'FontSize', 11, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontWeight','bold');
        end
        ylabel('right-to-left');
        xlabel('posterior-to-anterior');
        cb = colorbar;
        cb.Label.String = 'abs(correlation) with psychophysical response';
        saveas(fig1, [cd filesep 'plots_vep_psychophysics' filesep sID{i} '-vep_psychophysics-correlations.fig'])
        if i~=max(length(sID))
            input('Press enter to see the next plot\n')
        end

    end % end plotIndividuals

    all_r(i,:) = r;

end

if plotGroupmean == 1
    mean_r = mean(all_r);

    %%%%% figure
    fig3=figure(3); clf; hold on;
    if plotIndividuals==0
        set(fig3, 'Position', [100 100 600 800]);
    end
    % histogram of correlation values
    subplot(4,1,1); hold on; set(gca, 'FontSize', 14);
    nbins = 20;
    h=histogram(mean_r,nbins);
    h.FaceColor = [0.3 0.6 1];
    h.EdgeColor = 'k';
    h.LineWidth = 1.5;
    edges = h.BinEdges;
    for binIdx = 1:nbins % label electrodes
        binIndices = find(mean_r >= edges(binIdx) & mean_r < edges(binIdx + 1));
        binLabels = {vepdata.chanlocs(binIndices).labels};
        labelX = mean([edges(binIdx), edges(binIdx + 1)]);
        text(labelX, 0.02, strjoin(binLabels, '\n'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10);
    end
    ylim([0 max(h.BinCounts+1)]);
    ylabel('frequency')
    xlabel('mean correlation with psychophysical response')
    set(gca, 'YGrid', 'on');

    title([group ' group mean, ' titlestring], 'Interpreter','none')

    % electrode correlation heatmap
    subplot(4, 1, 2:4); hold on; set(gca, 'FontSize', 14);

    usecolormap = parula(64);
    usecolormap = rescale(usecolormap, 0.5, 1);
    colormap(usecolormap);

    color_index = ceil(abs(mean_r) * (size(usecolormap, 1)));
    col = usecolormap(color_index, :);

    scatter([vepdata.chanlocs.X], [vepdata.chanlocs.Y], 500, col, 'o', 'filled', 'MarkerEdgeColor', [.9 .9 .9], 'LineWidth', 2);
    for c = 1:size(vepdata.chanlocs,2)
        text(vepdata.chanlocs(c).X, vepdata.chanlocs(c).Y, strjoin({vepdata.chanlocs(c).labels, num2str(mean_r(c), 2)}, ': '),...
            'FontSize', 11, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontWeight','bold');
    end
    ylabel('right-to-left');
    xlabel('posterior-to-anterior');
    cb = colorbar;
    cb.Label.String = 'mean abs(correlation) with psychophysical response';
    saveas(fig3, [cd filesep 'plots_vep_psychophysics' filesep 'groupmean-vep_psychophysics-correlations.fig'])

end



