function p = plotShadyError(x, means, error, color, yscale)

    if strcmp(yscale, 'log')
        curve1 = max(means + error, 1);
        curve2 = max(means - error, 1);
        xaxis = [x; flipud(x)];
        
        inBetween = [curve1; flipud(curve2)];

        fill(xaxis, inBetween, color{1});
        hold on;

        means(means < 1) = 1;
        p = plot(x, means, color{2}, 'LineWidth', 1);
        set(gca, 'YScale', 'log' )

    else
        curve1 = means + error;
        curve2 = means - error;
        xaxis = [x; flipud(x)];
        inBetween = [curve1; flipud(curve2)];

        fill(xaxis, inBetween, color{1}, 'FaceAlpha', 0.5, 'LineStyle', 'none');
        hold on;

        p = plot(x, means, 'Color', color{2}, 'LineWidth', 1);
        set(gca, 'YScale', 'default')
    end