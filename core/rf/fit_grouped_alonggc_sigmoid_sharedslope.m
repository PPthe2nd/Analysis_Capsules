function F = fit_grouped_alonggc_sigmoid_sharedslope(G, varargin)
% FIT_GROUPED_ALONGGC_SIGMOID_SHAREDSLOPE
% Jointly fit one sigmoid per group with:
%   - shared slope parameter tau across groups
%   - group-specific amplitude A_g and midpoint t50_g
%   - fixed baseline at 0
%
% Model:
%   y_g(t) = A_g ./ (1 + exp(-(t - t50_g)/tau))
%
% Inputs:
%   G : grouped output struct from build_grouped_alonggc_polygons_allbins
%
% Name/value options:
%   'timeRef'            : 'start'|'center'|'end' (default 'center')
%   'fitStartMs'         : default 0
%   'fitEndMs'           : default 500
%   'smoothW'            : moving-average width in bins (default 3)
%   'useAbs'             : fit abs(trace) (default false)
%   'minPointsPerGroup'  : minimum finite bins in fit range (default 6)
%   'minTau'             : lower bound for shared tau (default 5)
%   'maxTau'             : upper bound for shared tau (default 400)
%   'initTau'            : initial shared tau (default 40)
%   't50PadMs'           : extra bound padding for t50 (default 100)
%   'maxAmpFactor'       : upper bound factor for amplitudes (default 3)
%   'plotMode'           : 'overlay'|'none' (default 'overlay')
%   'cmapName'           : default 'parula'
%   'lineWidth'          : raw/smoothed line width (default 1.8)
%   'fitLineWidthFactor' : fit width multiplier (default 2.2)
%   'rawFadeToWhite'     : 0..1 (default 0.65)
%   'showT50Markers'     : default true
%   't50MarkerSize'      : default 7
%   'showT50Regression'  : default true
%   'regressionX'        : 'alongMid'|'groupIdx'|'nComb' (default 'alongMid')
%   'excludeHighestXFromRegression' : default false
%   'onsetMs'            : default 0
%   'verbose'            : default true

p = inputParser;
p.addParameter('timeRef', 'center', @(x) ischar(x) || isstring(x));
p.addParameter('fitStartMs', 0, @(x) isnumeric(x) && isscalar(x));
p.addParameter('fitEndMs', 500, @(x) isnumeric(x) && isscalar(x));
p.addParameter('smoothW', 3, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('useAbs', false, @(x) islogical(x) && isscalar(x));
p.addParameter('minPointsPerGroup', 6, @(x) isnumeric(x) && isscalar(x) && x >= 3);
p.addParameter('minTau', 5, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('maxTau', 400, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('initTau', 40, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('t50PadMs', 100, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('maxAmpFactor', 3, @(x) isnumeric(x) && isscalar(x) && x > 1);
p.addParameter('plotMode', 'overlay', @(x) ischar(x) || isstring(x));
p.addParameter('cmapName', 'parula', @(x) ischar(x) || isstring(x));
p.addParameter('lineWidth', 1.8, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('fitLineWidthFactor', 2.2, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('rawFadeToWhite', 0.65, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
p.addParameter('showT50Markers', true, @(x) islogical(x) && isscalar(x));
p.addParameter('t50MarkerSize', 7, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('showT50Regression', true, @(x) islogical(x) && isscalar(x));
p.addParameter('regressionX', 'alongMid', @(x) ischar(x) || isstring(x));
p.addParameter('excludeHighestXFromRegression', false, @(x) islogical(x) && isscalar(x));
p.addParameter('onsetMs', 0, @(x) isnumeric(x) && isscalar(x));
p.addParameter('verbose', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;

assert(isstruct(G) && isfield(G, 'groupMeanSigned') && isfield(G, 'timeWindows'), ...
    'G must be grouped output from build_grouped_alonggc_polygons_allbins.');
assert(exist('lsqcurvefit', 'file') == 2, 'lsqcurvefit not found (Optimization Toolbox required).');

Y = double(G.groupMeanSigned);
if opt.useAbs
    Y = abs(Y);
end
[nGroups, nBins] = size(Y);
Yraw = Y;

smoothW = max(1, round(opt.smoothW));
if smoothW > 1
    Y = smooth_rows_movmean_omitnan(Y, smoothW);
end

tw = double(G.timeWindows);
assert(size(tw,1) == nBins && size(tw,2) == 2, 'G.timeWindows must be [nBins x 2].');
switch lower(string(opt.timeRef))
    case "start"
        t = tw(:,1);
    case "end"
        t = tw(:,2);
    otherwise
        t = mean(tw,2);
end
t = t(:);

fitMask = isfinite(t) & (t >= opt.fitStartMs) & (t <= opt.fitEndMs);
assert(any(fitMask), 'No bins in fit range [%.1f, %.1f] ms.', opt.fitStartMs, opt.fitEndMs);

% Build joint dataset over groups with enough points.
grpValid = false(nGroups,1);
grpLocal = nan(nGroups,1);
tData = [];
yData = [];
gData = [];
A0 = [];
t500 = [];

for g = 1:nGroups
    yg = Y(g,:).';
    ok = fitMask & isfinite(yg);
    if nnz(ok) < round(opt.minPointsPerGroup)
        continue;
    end

    tFit = t(ok);
    yFit = yg(ok);
    m = max(yFit, [], 'omitnan');
    if ~isfinite(m)
        continue;
    end
    m = max(m, 1e-4);

    yHalf = 0.5 * m;
    iHalf = find(yFit >= yHalf, 1, 'first');
    if isempty(iHalf)
        t50init = median(tFit, 'omitnan');
    else
        t50init = tFit(iHalf);
    end

    grpValid(g) = true;
    grpLocal(g) = numel(A0) + 1;
    A0(end+1,1) = m; %#ok<AGROW>
    t500(end+1,1) = t50init; %#ok<AGROW>

    tData = [tData; tFit]; %#ok<AGROW>
    yData = [yData; yFit]; %#ok<AGROW>
    gData = [gData; repmat(grpLocal(g), numel(tFit), 1)]; %#ok<AGROW>
end

nFitGroups = nnz(grpValid);
assert(nFitGroups >= 1, 'No groups had enough finite points for sigmoid fit.');

tau0 = min(max(opt.initTau, opt.minTau), opt.maxTau);
p0 = [A0; t500; tau0];

ampGlobal = max(yData, [], 'omitnan');
if ~isfinite(ampGlobal) || ampGlobal <= 0
    ampGlobal = max(A0, [], 'omitnan');
end
if ~isfinite(ampGlobal) || ampGlobal <= 0
    ampGlobal = 1;
end
Aub = max(1e-3, opt.maxAmpFactor * ampGlobal);

tMin = min(tData, [], 'omitnan');
tMax = max(tData, [], 'omitnan');
lb = [zeros(nFitGroups,1); repmat(tMin - opt.t50PadMs, nFitGroups,1); opt.minTau];
ub = [repmat(Aub, nFitGroups,1); repmat(tMax + opt.t50PadMs, nFitGroups,1); opt.maxTau];

X = [tData, gData];
optsLSQ = optimoptions('lsqcurvefit', 'Display', 'off');
modelFun = @(pp,xx) sigmoid_shared_model(pp, xx, nFitGroups);
[pHat, ~, residual, exitflag, output] = lsqcurvefit(modelFun, p0, X, yData, lb, ub, optsLSQ); %#ok<ASGLU>

Ahat = pHat(1:nFitGroups);
t50hatLocal = pHat(nFitGroups+1 : 2*nFitGroups);
tauHat = pHat(end);

% Expand fitted params back to all groups.
AhatAll = nan(nGroups,1);
t50hatAll = nan(nGroups,1);
rmse = nan(nGroups,1);
nFitPts = zeros(nGroups,1);
yFitAll = nan(nGroups, nBins);
status = strings(nGroups,1);
status(:) = "not_fitted";
regX = nan(nGroups,1);

for g = 1:nGroups
    if ~grpValid(g)
        continue;
    end
    k = grpLocal(g);
    AhatAll(g) = Ahat(k);
    t50hatAll(g) = t50hatLocal(k);
    yFitAll(g,:) = Ahat(k) ./ (1 + exp(-(t(:).' - t50hatLocal(k)) ./ tauHat));

    yg = Y(g,:).';
    ok = fitMask & isfinite(yg);
    nFitPts(g) = nnz(ok);
    if nFitPts(g) > 0
        e = yFitAll(g,ok) - yg(ok).';
        rmse(g) = sqrt(mean(e.^2, 'omitnan'));
    end
    status(g) = "ok";
end

% Regression x-variable (default: along_GC midpoint per group).
xMode = lower(string(opt.regressionX));
if isfield(G, 'groupSummary') && istable(G.groupSummary) && height(G.groupSummary) >= nGroups
    switch xMode
        case "ncomb"
            regX = double(G.groupSummary.nComb(1:nGroups));
        case "groupidx"
            regX = (1:nGroups).';
        otherwise
            if ismember('alongMin', G.groupSummary.Properties.VariableNames) && ...
                    ismember('alongMax', G.groupSummary.Properties.VariableNames)
                regX = 0.5 * (double(G.groupSummary.alongMin(1:nGroups)) + ...
                              double(G.groupSummary.alongMax(1:nGroups)));
            else
                regX = (1:nGroups).';
            end
    end
else
    regX = (1:nGroups).';
end

% Regression mask with optional exclusion of the highest-x point.
regMask = isfinite(regX) & isfinite(t50hatAll);
excludedGroupIdx = NaN;
if opt.excludeHighestXFromRegression && nnz(regMask) >= 3
    idxGood = find(regMask);
    xGood = regX(idxGood);
    [~, iMax] = max(xGood);
    excludedGroupIdx = idxGood(iMax);
    regMask(excludedGroupIdx) = false;
end

regN = nnz(regMask);
regSlope = NaN;
regIntercept = NaN;
regR2 = NaN;
if regN >= 2
    pLin = polyfit(regX(regMask), t50hatAll(regMask), 1);
    regSlope = pLin(1);
    regIntercept = pLin(2);
    yHat = polyval(pLin, regX(regMask));
    ssRes = sum((t50hatAll(regMask) - yHat).^2);
    ssTot = sum((t50hatAll(regMask) - mean(t50hatAll(regMask))).^2);
    if ssTot > 0
        regR2 = 1 - ssRes/ssTot;
    end
else
    pLin = [NaN NaN];
end

% Optional plot
fig = [];
ax = [];
figReg = [];
axReg = [];
mode = lower(string(opt.plotMode));
if mode == "overlay"
    cmapFun = str2func(char(opt.cmapName));
    try
        C = cmapFun(nGroups);
    catch
        C = parula(nGroups);
    end

    fig = figure('Color', 'w');
    ax = axes('Parent', fig); hold(ax, 'on');
    fitLW = max(opt.lineWidth * opt.fitLineWidthFactor, opt.lineWidth + 0.8);
    rawW = max(0.8, 0.9 * opt.lineWidth);

    for g = 1:nGroups
        lbl = sprintf('G%d', g);
        if isfield(G, 'groupSummary') && istable(G.groupSummary) && height(G.groupSummary) >= g
            a0 = double(G.groupSummary.alongMin(g));
            a1 = double(G.groupSummary.alongMax(g));
            lbl = sprintf('G%d [%.2f, %.2f]', g, a0, a1);
        end

        rawCol = (1 - opt.rawFadeToWhite) * C(g,:) + opt.rawFadeToWhite * [1 1 1];
        plot(ax, t, Y(g,:), '-', 'Color', rawCol, 'LineWidth', rawW, 'DisplayName', lbl);

        if all(isfinite(yFitAll(g,:)))
            plot(ax, t, yFitAll(g,:), '-', 'Color', C(g,:), 'LineWidth', fitLW, ...
                'HandleVisibility', 'off');
            if opt.showT50Markers && isfinite(t50hatAll(g)) && isfinite(AhatAll(g))
                plot(ax, t50hatAll(g), 0, 'o', ...
                    'MarkerSize', opt.t50MarkerSize, ...
                    'MarkerFaceColor', C(g,:), ...
                    'MarkerEdgeColor', [0.1 0.1 0.1], ...
                    'LineWidth', 0.8, ...
                    'HandleVisibility', 'off');
            end
        end
    end

    xline(ax, opt.onsetMs, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.0, 'DisplayName', 'onset');
    yline(ax, 0, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.0, 'DisplayName', 'zero');
    grid(ax, 'on');
    box(ax, 'off');
    xlabel(ax, 'Time (ms)');
    if opt.useAbs
        ylabel(ax, 'Mean |\Delta(T-D)|');
    else
        ylabel(ax, 'Mean signed \Delta (T-D)');
    end
    title(ax, sprintf('Shared-slope sigmoid fits (tau=%.2f ms, smoothW=%d)', tauHat, smoothW));
    legend(ax, 'Location', 'eastoutside');

    if opt.showT50Regression
        if regN >= 2
            xr = regX(regMask);
            figReg = figure('Color', 'w');
            axReg = axes('Parent', figReg); hold(axReg, 'on');
            for g = 1:nGroups
                if isfinite(regX(g)) && isfinite(t50hatAll(g))
                    if regMask(g)
                        mkFace = C(g,:);
                        mkEdge = [0.1 0.1 0.1];
                    else
                        mkFace = [1 1 1];
                        mkEdge = [0.35 0.35 0.35];
                    end
                    plot(axReg, regX(g), t50hatAll(g), 'o', ...
                        'MarkerSize', opt.t50MarkerSize, ...
                        'MarkerFaceColor', mkFace, ...
                        'MarkerEdgeColor', mkEdge, ...
                        'LineWidth', 0.8);
                end
            end
            xs = linspace(min(xr), max(xr), 200);
            ys = polyval(pLin, xs);
            plot(axReg, xs, ys, '-', 'Color', [0.1 0.1 0.1], 'LineWidth', 2.0);
            grid(axReg, 'on');
            box(axReg, 'off');

            switch xMode
                case "ncomb"
                    xlabel(axReg, 'Number of growth cones (nComb)');
                case "groupidx"
                    xlabel(axReg, 'Group index');
                otherwise
                    xlabel(axReg, 'Number of growth cones (along\_GC midpoint)');
            end
            ylabel(axReg, 't50 (ms)');
            title(axReg, sprintf('t50 regression: slope=%.3g ms/unit, R^2=%.3f', regSlope, regR2));
        end
    end
end

summary = table((1:nGroups)', AhatAll, t50hatAll, repmat(tauHat, nGroups,1), ...
    0.5*AhatAll, rmse, nFitPts, status, ...
    'VariableNames', {'groupIdx','A','t50','tauShared','yHalfModel','rmse','nFitPoints','status'});

if isfield(G, 'groupSummary') && istable(G.groupSummary) && height(G.groupSummary) >= nGroups
    summary.alongMin = double(G.groupSummary.alongMin(1:nGroups));
    summary.alongMax = double(G.groupSummary.alongMax(1:nGroups));
end

F = struct();
F.summary = summary;
F.sharedTau = tauHat;
F.paramsLocal = pHat;
F.groupValid = grpValid;
F.groupLocalIndex = grpLocal;
F.t = t;
F.y = Y;
F.yRaw = Yraw;
F.yFit = yFitAll;
F.fitMask = fitMask;
F.options = opt;
F.exitflag = exitflag;
F.output = output;
F.residual = residual;
F.fig = fig;
F.ax = ax;
F.regressionX = regX;
F.figRegression = figReg;
F.axRegression = axReg;
F.regressionMask = regMask;
F.regression = struct( ...
    'mode', char(xMode), ...
    'slope', regSlope, ...
    'intercept', regIntercept, ...
    'r2', regR2, ...
    'n', regN, ...
    'excludeHighestXFromRegression', opt.excludeHighestXFromRegression, ...
    'excludedGroupIdx', excludedGroupIdx);

if opt.verbose
    fprintf('Shared-slope sigmoid fit: %d/%d groups fitted | tau=%.4g ms\n', ...
        nFitGroups, nGroups, tauHat);
    disp(summary(:, {'groupIdx','A','t50','tauShared','rmse','nFitPoints','status'}));
    if isstruct(F.regression)
        fprintf('t50 regression (%s): slope=%.6g, intercept=%.6g, R^2=%.4f, n=%d\n', ...
            F.regression.mode, F.regression.slope, F.regression.intercept, F.regression.r2, F.regression.n);
        if isfinite(F.regression.excludedGroupIdx)
            fprintf('  excluded highest-x group from regression: G%d\n', F.regression.excludedGroupIdx);
        end
    end
end

end

function y = sigmoid_shared_model(pp, X, nGroupsFit)
t = X(:,1);
g = round(X(:,2));
A = pp(1:nGroupsFit);
t50 = pp(nGroupsFit+1 : 2*nGroupsFit);
tau = pp(end);

g(g < 1) = 1;
g(g > nGroupsFit) = nGroupsFit;
y = A(g) ./ (1 + exp(-(t - t50(g)) ./ tau));
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
