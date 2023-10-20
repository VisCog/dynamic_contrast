clear;

% plot settings
nplotcols = 4;
ntrialspercol = 14;
sCond = 'congruent';

files = dir(['*-motor-' sCond '*']);



for f = 1:length(files)
    load(files(f).name);

    switch lower(sCond)
        case 'congruent'
            data = congruentMotor;
        case 'orthogonal'
            data = orthogonalMotor;
    end

    bstart = data.experiment.binoResponseStart;
    middle = data.experiment.response;
    bend = data.experiment.binoResponseEnd;
    response = [bstart middle bend];
    clear bstart middle bend

    bstart = data.experiment.BEcontrastStart;
    lcontrast = data.experiment.LEcontrast;
    rcontrast = data.experiment.REcontrast;
    bend = data.experiment.BEcontrastEnd;
    contrastLE = [bstart lcontrast bend];
    contrastRE = [bstart rcontrast bend];
    clear bstart lcontrast rcontrast bend

    t = size(response, 2);

    fig = figure(1);
    set(fig, 'Position', [0 0 900 800]);
    clf;
    topAxes = axes('Position', [0.1, 0.4, 0.8, 0.55]);
    bottomAxes = axes('Position', [0.1, 0.1, 0.8, 0.25]);


    % Plot in the top tile
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
    xlim(topAxes, [-300 ceil(t * g * 1.2 / 1000) * 1000]);
    ylim(topAxes, [round(-ntrialspercol * 1.2 - 1) 2]);
    text(100, 1, 'individual trials:', 'color', 'k');
    text(1000, 1, 'response in black', 'color', 'k', 'FontWeight', 'bold');
    text(2000, 1, 'left contrast in blue', 'color', [0, 0, 0.5], 'FontWeight', 'bold');
    text(3000, 1, 'right contrast in gold', 'color', [0.5, 0.4, 0], 'FontWeight', 'bold');
    title(topAxes, files(f).name, 'Interpreter', 'none');


    % Plot in the bottom tile
     axes(bottomAxes); hold on;

    leftEyeFast = nanmean(response(data.config.fastEye == 0 & data.experiment.joystickActuallyUsed == 1, :));
    rightEyeFast = nanmean(response(data.config.fastEye == 1 & data.experiment.joystickActuallyUsed == 1, :));
    
    % denote bino sections
    xline(size(data.experiment.BEcontrastStart, 2), 'k--', 'Color', [0.5, 0.5, 0.5], 'Parent', bottomAxes);
    xline(t - size(data.experiment.BEcontrastStart, 2), 'k--', 'Color', [0.5, 0.5, 0.5], 'Parent', bottomAxes);


    fastSin = ((sin(2*pi*([1:length(data.experiment.response)]/30)/6)+1)/2);
    slowSin = ((sin(2*pi*([1:length(data.experiment.response)]/30)/8)+1)/2);

    fastSin = [data.experiment.BEcontrastStart(1,:) fastSin data.experiment.BEcontrastEnd(1,:)];
    slowSin = [data.experiment.BEcontrastStart(1,:) slowSin data.experiment.BEcontrastEnd(1,:)];

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
    input('Press enter to look at next participant');
end
