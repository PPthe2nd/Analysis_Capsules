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
p.addParameter('smoothW', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('halfMaxFrac', 0.5, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
p.addParameter('halfMaxSearchStartMs', 0, @(x) isnumeric(x) && isscalar(x));
p.addParameter('showHalfMaxMarkers', true, @(x) islogical(x) && isscalar(x));
p.addParameter('halfMaxMarkerSize', 6.5, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.parse(varargin{:});
opt = p.Results;

assert(isstruct(G) && isfield(G, 'groupMeanSigned') && isfield(G, 'timeWindows'), ...
    'G must be grouped output from build_grouped_alonggc_polygons_allbins.');

Y = double(G.groupMeanSigned);
if opt.useAbs
    Y = abs(Y);
end
[nGroups, nBins] = size(Y);
Yplot = Y;
smoothW = max(1, round(opt.smoothW));
if smoothW > 1
    Yplot = smooth_rows_movmean_omitnan(Yplot, smoothW);
end

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

% First crossing time to half-max (on smoothed curves, search window constrained by time).
tHalf = nan(nGroups,1);
yHalf = nan(nGroups,1);
yMax = nan(nGroups,1);
iHalf = nan(nGroups,1);
for g = 1:nGroups
    yg = Yplot(g,:).';
    mSearch = isfinite(t) & isfinite(yg) & (t >= opt.halfMaxSearchStartMs);
    if ~any(mSearch)
        continue;
    end
    tS = t(mSearch);
    yS = yg(mSearch);
    [mx, iMx] = max(yS);
    if ~isfinite(mx) || mx <= 0 || isempty(iMx)
        continue;
    end
    yMax(g) = mx;
    yH = opt.halfMaxFrac * mx;
    yHalf(g) = yH;
    i0 = find(yS >= yH, 1, 'first');
    if isempty(i0)
        continue;
    end
    if i0 == 1
        tHalf(g) = tS(1);
        iHalf(g) = find(mSearch, 1, 'first');
    else
        t1 = tS(i0-1); t2 = tS(i0);
        y1 = yS(i0-1); y2 = yS(i0);
        if isfinite(y2-y1) && abs(y2-y1) > 0
            tHalf(g) = t1 + (yH - y1) * (t2 - t1) / (y2 - y1);
        else
            tHalf(g) = t2;
        end
        idxGlobal = find(mSearch);
        iHalf(g) = idxGlobal(i0);
    end
end

v = Yplot(isfinite(Yplot));
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
        plot(ax, t, Yplot(g,:), 'Color', C(g,:), 'LineWidth', opt.lineWidth, 'DisplayName', lbl);
        if opt.showHalfMaxMarkers && isfinite(tHalf(g)) && isfinite(yHalf(g))
            plot(ax, tHalf(g), yHalf(g), 'o', ...
                'MarkerSize', opt.halfMaxMarkerSize, ...
                'MarkerFaceColor', C(g,:), ...
                'MarkerEdgeColor', [0.1 0.1 0.1], ...
                'LineWidth', 0.8, ...
                'HandleVisibility', 'off');
        end
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
        title(ax, sprintf('Grouped mean activity time-courses (all groups, smoothW=%d)', smoothW));
    else
        ylabel(ax, 'Mean signed \Delta (T-D)');
        title(ax, sprintf('Grouped mean signed \Delta (T-D) time-courses (all groups, smoothW=%d)', smoothW));
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
        plot(ax(g), t, Yplot(g,:), 'Color', C(g,:), 'LineWidth', opt.lineWidth);
        if opt.showHalfMaxMarkers && isfinite(tHalf(g)) && isfinite(yHalf(g))
            plot(ax(g), tHalf(g), yHalf(g), 'o', ...
                'MarkerSize', opt.halfMaxMarkerSize, ...
                'MarkerFaceColor', C(g,:), ...
                'MarkerEdgeColor', [0.1 0.1 0.1], ...
                'LineWidth', 0.8);
        end
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
        title(tl, sprintf('Grouped mean activity time-courses (smoothW=%d)', smoothW));
    else
        ylabel(tl, 'Mean signed \Delta (T-D)');
        title(tl, sprintf('Grouped mean signed \Delta (T-D) time-courses (smoothW=%d)', smoothW));
    end

    h = struct();
    h.fig = fig;
    h.tiled = tl;
    h.ax = ax;
end

h.t = t;
h.Y = Yplot;
h.YRaw = Y;
h.smoothW = smoothW;
h.halfMaxFrac = opt.halfMaxFrac;
h.halfMaxSearchStartMs = opt.halfMaxSearchStartMs;
h.tHalf = tHalf;
h.yHalf = yHalf;
h.yMax = yMax;
h.iHalf = iHalf;
h.summary = table((1:nGroups)', tHalf, yHalf, yMax, ...
    'VariableNames', {'groupIdx','tHalf','yHalf','yMaxSmooth'});
h.colors = C;
h.yL = yL;

end

function Ysm = smooth_rows_movmean_omitnan(Y, w)
Ysm = nan(size(Y));
w = max(1, round(w));
if w <= 1
    Ysm = Y;
    return;
end
[nRows, nCols] = size(Y);
halfLo = floor((w-1)/2);
halfHi = ceil((w-1)/2);
for r = 1:nRows
    yr = Y(r,:);
    for c = 1:nCols
        i1 = max(1, c-halfLo);
        i2 = min(nCols, c+halfHi);
        v = yr(i1:i2);
        v = v(isfinite(v));
        if ~isempty(v)
            Ysm(r,c) = mean(v);
        end
    end
end
end
