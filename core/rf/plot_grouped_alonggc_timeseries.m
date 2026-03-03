function h = plot_grouped_alonggc_timeseries(G, varargin)
% PLOT_GROUPED_ALONGGC_TIMESERIES
% Plot one time-course subplot per along_GC group using a smooth colormap.

p = inputParser;
p.addParameter('plotMode', 'subplots', @(x) ischar(x) || isstring(x)); % 'subplots'|'overlay'
p.addParameter('timeRef', 'center', @(x) ischar(x) || isstring(x)); % 'start'|'center'|'end'
p.addParameter('onsetMs', 0, @(x) isnumeric(x) && isscalar(x));
p.addParameter('cmapName', 'parula', @(x) ischar(x) || isstring(x));
p.addParameter('lineWidth', 1.8, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('useAbs', false, @(x) islogical(x) && isscalar(x));
p.addParameter('layoutCols', 2, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.parse(varargin{:});
opt = p.Results;

assert(isstruct(G) && isfield(G, 'groupMeanSigned') && isfield(G, 'timeWindows'), ...
    'G must be grouped output from build_grouped_alonggc_polygons_allbins.');

Y = double(G.groupMeanSigned);
if opt.useAbs
    Y = abs(Y);
end
[nGroups, nBins] = size(Y);

tw = double(G.timeWindows);
assert(size(tw,1) == nBins && size(tw,2) == 2, 'G.timeWindows must be [nBins x 2].');

ref = lower(string(opt.timeRef));
switch ref
    case "start"
        t = tw(:,1);
    case "end"
        t = tw(:,2);
    otherwise
        t = mean(tw,2);
end

cmapFun = str2func(char(opt.cmapName));
try
    C = cmapFun(nGroups);
catch
    C = parula(nGroups);
end

v = Y(isfinite(Y));
if isempty(v)
    yL = [-1 1];
else
    lo = min(v);
    hi = max(v);
    if lo == hi
        pad = max(1e-3, 0.1*abs(lo) + 1e-3);
        yL = [lo-pad, hi+pad];
    else
        pad = 0.08*(hi-lo);
        yL = [lo-pad, hi+pad];
    end
end

nCols = max(1, round(opt.layoutCols));
nRows = ceil(nGroups / nCols);

fig = figure('Color', 'w');
mode = lower(string(opt.plotMode));
if mode == "overlay"
    ax = axes('Parent', fig); hold(ax, 'on');
    for g = 1:nGroups
        lbl = sprintf('G%d', g);
        if isfield(G, 'groupSummary') && istable(G.groupSummary) && height(G.groupSummary) >= g
            a0 = double(G.groupSummary.alongMin(g));
            a1 = double(G.groupSummary.alongMax(g));
            lbl = sprintf('G%d [%.2f, %.2f]', g, a0, a1);
        end
        plot(ax, t, Y(g,:), 'Color', C(g,:), 'LineWidth', opt.lineWidth, 'DisplayName', lbl);
    end
    xline(ax, opt.onsetMs, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.0, 'DisplayName', 'onset');
    if ~opt.useAbs
        yline(ax, 0, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.0, 'DisplayName', 'zero');
    end
    xlim(ax, [min(t) max(t)]);
    ylim(ax, yL);
    grid(ax, 'on');
    set(ax, 'XColor', [0.2 0.2 0.2], 'YColor', [0.2 0.2 0.2], 'Box', 'off');
    xlabel(ax, 'Time (ms)');
    if opt.useAbs
        ylabel(ax, 'Mean |T-D|');
        title(ax, 'Grouped mean activity time-courses (all groups)');
    else
        ylabel(ax, 'Mean signed \Delta (T-D)');
        title(ax, 'Grouped mean signed \Delta (T-D) time-courses (all groups)');
    end
    legend(ax, 'Location', 'eastoutside');

    h = struct();
    h.fig = fig;
    h.tiled = [];
    h.ax = ax;
else
    tl = tiledlayout(nRows, nCols, 'Padding', 'compact', 'TileSpacing', 'compact');
    ax = gobjects(nGroups,1);
    for g = 1:nGroups
        ax(g) = nexttile(tl, g);
        hold(ax(g), 'on');
        plot(ax(g), t, Y(g,:), 'Color', C(g,:), 'LineWidth', opt.lineWidth);
        xline(ax(g), opt.onsetMs, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.0);
        if ~opt.useAbs
            yline(ax(g), 0, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.0);
        end
        xlim(ax(g), [min(t) max(t)]);
        ylim(ax(g), yL);
        grid(ax(g), 'on');
        set(ax(g), 'XColor', [0.2 0.2 0.2], 'YColor', [0.2 0.2 0.2], 'Box', 'off');

        ttl = sprintf('Group %d', g);
        if isfield(G, 'groupSummary') && istable(G.groupSummary) && height(G.groupSummary) >= g
            nComb = double(G.groupSummary.nComb(g));
            a0 = double(G.groupSummary.alongMin(g));
            a1 = double(G.groupSummary.alongMax(g));
            ttl = sprintf('Group %d | n=%d | along %.2f-%.2f', g, nComb, a0, a1);
        end
        title(ax(g), ttl, 'Color', C(g,:), 'FontWeight', 'bold', 'FontSize', 9);
    end

    xlabel(tl, 'Time (ms)');
    if opt.useAbs
        ylabel(tl, 'Mean |T-D|');
        title(tl, 'Grouped mean activity time-courses');
    else
        ylabel(tl, 'Mean signed \Delta (T-D)');
        title(tl, 'Grouped mean signed \Delta (T-D) time-courses');
    end

    h = struct();
    h.fig = fig;
    h.tiled = tl;
    h.ax = ax;
end

h.t = t;
h.Y = Y;
h.colors = C;
h.yL = yL;

end
