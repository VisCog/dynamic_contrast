clear;

% plot settings
nplotcols = 2;
ntrialspercol = 14;

files = dir('*-vep-data.mat');

for f = 1:length(files)

    try
        load(files(f).name);

        bstart = dataVep.experiment.binoResponseStart;
        middle = dataVep.experiment.response;
        bend = dataVep.experiment.binoResponseEnd;
        response = [bstart middle bend];
        clear bstart middle bend

        bstart = dataVep.experiment.BEcontrastStart;
        lcontrast = dataVep.experiment.LEcontrast;
        rcontrast = dataVep.experiment.REcontrast;
        bend = dataVep.experiment.BEcontrastEnd;
        contrastLE = [bstart lcontrast bend];
        contrastRE = [bstart rcontrast bend];
        clear bstart lcontrast rcontrast bend

        t = size(response, 2);

        fig = figure(1);
        set(fig, 'Position', [0 0 900 800]);
        clf;
        topAxes = axes('Position', [0.1, 0.4, 0.8, 0.55]);
        bottomAxes = axes('Position', [0.1, 0.1, 0.8, 0.25]);

        % Plot in the top tile - each trial, by column
        axes(topAxes);
        hold on;
        for g = 1:nplotcols
            for i = 1:ntrialspercol
                idx = (g - 1) * ntrialspercol + i;
                rowspace = (ntrialspercol * (g - 1)) - (i - 1) * 0.2;
                plot([1:t] + (t * 1.4 * (g - 1)), contrastLE(idx, :) - idx + rowspace, 'Color', [.5 .5 .7], 'LineWidth', 2);
                plot([1:t] + (t * 1.4 * (g - 1)), contrastRE(idx, :) - idx + rowspace, 'Color', [0.8 .8 .5], 'LineWidth', 2);
                plot([1:t] + (t * 1.4 * (g - 1)), response(idx, :) - idx + rowspace, 'k', 'LineWidth', 2);
                text(-100 + (t * 1.4 * (g - 1)), 0.8 - idx + rowspace, num2str(idx));
            end
        end

        set(topAxes, 'XTickLabel', '');
        set(topAxes, 'YTickLabel', '');
        xlim(topAxes, [-300 ceil(t * g * 1.2 / 100) * 100]);
        ylim(topAxes, [round(-ntrialspercol * 1.2 - 1) 2]);
        text(100, 1, 'individual trials:', 'color', 'k');
        text(600, 1, 'response in black', 'color', 'k', 'FontWeight', 'bold');
        text(1200, 1, 'left contrast in blue', 'color', [0, 0, 0.5], 'FontWeight', 'bold');
        text(1800, 1, 'right contrast in gold', 'color', [0.5, 0.4, 0], 'FontWeight', 'bold');
        title(topAxes, files(f).name, 'Interpreter', 'none');


        % Plot in the bottom tile - averages for left vs right fast
        axes(bottomAxes); hold on;

        leftEyeFast = nanmean(response(dataVep.config.fastEye == 0, :));
        rightEyeFast = nanmean(response(dataVep.config.fastEye == 1, :));

        % denote bino sections
        xline(size(dataVep.experiment.BEcontrastStart, 2), 'k--', 'Color', [0.5, 0.5, 0.5], 'Parent', bottomAxes);
        xline(t - size(dataVep.experiment.BEcontrastStart, 2), 'k--', 'Color', [0.5, 0.5, 0.5], 'Parent', bottomAxes);


        fastSin = ((sin(2*pi*([1:length(dataVep.experiment.response)]/30)/6)+1)/2);
        slowSin = ((sin(2*pi*([1:length(dataVep.experiment.response)]/30)/8)+1)/2);

        fastSin = [dataVep.experiment.BEcontrastStart(1,:) fastSin dataVep.experiment.BEcontrastEnd(1,:)];
        slowSin = [dataVep.experiment.BEcontrastStart(1,:) slowSin dataVep.experiment.BEcontrastEnd(1,:)];

        c1=plot(bottomAxes,1:t, fastSin, 'Color', [.7 .7 .7], 'LineWidth',2,'DisplayName','fast eye');
        c2=plot(bottomAxes,1:t, slowSin, 'Color', [.7 .7 .7], 'LineWidth',2, 'LineStyle','--', 'DisplayName','slow eye');

        p1=plot(bottomAxes, 1:t, leftEyeFast, 'Color', [0, 0, 0.5], 'LineWidth', 2, 'DisplayName', 'left fast mean');
        p2=plot(bottomAxes, 1:t, rightEyeFast, 'Color', [0.5, 0.4, 0], 'LineWidth', 2, 'DisplayName', 'right fast mean');

        xlim(bottomAxes, [1, t]);
        ylim(bottomAxes, [0, 1]);
        ylabel(bottomAxes, 'response');
        xlabel(bottomAxes, 'sample t');
        legend([c1 c2 p1 p2], 'Location', 'EastOutside');

        % Save the figure
        saveas(fig, [cd filesep 'plots' filesep strrep(files(f).name, '.mat', '.fig')]);
        if ~(f==length(files))
            input('Press enter to look at next participant');
        end
    catch ME
        ME.message
        disp(['Error with ' files(f).name ', skipping'])
    end

end

