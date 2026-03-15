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
%   'useT50SEWeightedRegression' : default true
%   'showT50SEErrorBars' : default true
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
p.addParameter('useT50SEWeightedRegression', true, @(x) islogical(x) && isscalar(x));
p.addParameter('showT50SEErrorBars', true, @(x) islogical(x) && isscalar(x));
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
t50ResidSD = nan(nGroups,1);
t50SlopeAtHalf = nan(nGroups,1);
t50SE = nan(nGroups,1);
regWeightAll = nan(nGroups,1);
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
        dofLocal = nFitPts(g) - 2;
        if dofLocal > 0
            rssLocal = sum(e.^2, 'omitnan');
            if isfinite(rssLocal) && rssLocal >= 0
                t50ResidSD(g) = sqrt(rssLocal / dofLocal);
            end
        end
    end
    if isfinite(AhatAll(g)) && isfinite(tauHat) && tauHat > 0
        t50SlopeAtHalf(g) = abs(AhatAll(g)) / (4 * tauHat);
    end
    if isfinite(t50ResidSD(g)) && isfinite(t50SlopeAtHalf(g)) && t50SlopeAtHalf(g) > 0
        t50SE(g) = t50ResidSD(g) / t50SlopeAtHalf(g);
    end
    status(g) = "ok";
end

regWeightAll = 1 ./ (t50SE.^2);
regWeightAll(~isfinite(regWeightAll) | regWeightAll <= 0) = NaN;

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
regMaskBase = isfinite(regX) & isfinite(t50hatAll);
regMaskWeighted = regMaskBase & isfinite(regWeightAll);
regUseWeighted = opt.useT50SEWeightedRegression && (nnz(regMaskWeighted) >= 2);
if regUseWeighted
    regMask = regMaskWeighted;
else
    regMask = regMaskBase;
end
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
regP = NaN;
regT = NaN;
regDf = NaN;
regSlopeSE = NaN;
regWeightMode = "unweighted";
if regN >= 2
    xr = regX(regMask);
    yr = t50hatAll(regMask);
    if regUseWeighted
        wr = regWeightAll(regMask);
        [regIntercept, regSlope, regSlopeSE, regP, regT, regDf, regR2] = ...
            weighted_linear_regression(xr, yr, wr);
        pLin = [regSlope regIntercept];
        regWeightMode = "inverse_t50se2";
    else
        pLin = polyfit(xr, yr, 1);
        regSlope = pLin(1);
        regIntercept = pLin(2);
        yHat = polyval(pLin, xr);
        ssRes = sum((yr - yHat).^2);
        ssTot = sum((yr - mean(yr)).^2);
        if ssTot > 0
            regR2 = 1 - ssRes/ssTot;
        end

        if regN >= 3
            xCtr = xr - mean(xr);
            sxx = sum(xCtr.^2);
            regDf = regN - 2;
            if isfinite(sxx) && sxx > 0 && isfinite(ssRes) && regDf > 0
                mse = ssRes / regDf;
                if isfinite(mse) && mse >= 0
                    regSlopeSE = sqrt(mse / sxx);
                    if isfinite(regSlopeSE) && regSlopeSE > 0
                        regT = regSlope / regSlopeSE;
                        xBeta = regDf / (regDf + regT.^2);
                        regP = betainc(xBeta, regDf/2, 0.5);
                    elseif isfinite(regSlope) && regSlope == 0
                        regT = 0;
                        regP = 1;
                    end
                end
            end
        end
        regWeightMode = "unweighted";
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
                    if opt.showT50SEErrorBars && isfinite(t50SE(g)) && t50SE(g) > 0
                        errorbar(axReg, regX(g), t50hatAll(g), t50SE(g), ...
                            'Color', [0.65 0.65 0.65], 'LineWidth', 1.0, 'CapSize', 0, ...
                            'HandleVisibility', 'off');
                    end
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
            if regUseWeighted
                title(axReg, sprintf('Weighted t50 regression: slope=%.3g ms/unit, p=%s, R^2=%.3f', ...
                    regSlope, format_p_value(regP), regR2));
            else
                title(axReg, sprintf('t50 regression: slope=%.3g ms/unit, p=%s, R^2=%.3f', ...
                    regSlope, format_p_value(regP), regR2));
            end
        end
    end
end

summary = table((1:nGroups)', AhatAll, t50hatAll, repmat(tauHat, nGroups,1), ...
    0.5*AhatAll, rmse, nFitPts, t50ResidSD, t50SlopeAtHalf, t50SE, regWeightAll, status, ...
    'VariableNames', {'groupIdx','A','t50','tauShared','yHalfModel','rmse','nFitPoints', ...
    'residSD','slopeAtT50','t50SE','regWeight','status'});

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
F.t50ResidSD = t50ResidSD;
F.t50SlopeAtHalf = t50SlopeAtHalf;
F.t50SE = t50SE;
F.regressionWeights = regWeightAll;
F.regression = struct( ...
    'mode', char(xMode), ...
    'slope', regSlope, ...
    'intercept', regIntercept, ...
    'r2', regR2, ...
    'pValueSlope', regP, ...
    'tStatSlope', regT, ...
    'df', regDf, ...
    'stderrSlope', regSlopeSE, ...
    'n', regN, ...
    'weightMode', char(regWeightMode), ...
    'excludeHighestXFromRegression', opt.excludeHighestXFromRegression, ...
    'excludedGroupIdx', excludedGroupIdx);

if opt.verbose
    fprintf('Shared-slope sigmoid fit: %d/%d groups fitted | tau=%.4g ms\n', ...
        nFitGroups, nGroups, tauHat);
    disp(summary(:, {'groupIdx','A','t50','t50SE','tauShared','rmse','nFitPoints','status'}));
    if isstruct(F.regression)
        fprintf(['t50 regression (%s, %s): slope=%.6g, intercept=%.6g, p=%s, ' ...
                 'R^2=%.4f, n=%d\n'], ...
            F.regression.mode, F.regression.weightMode, ...
            F.regression.slope, F.regression.intercept, ...
            format_p_value(F.regression.pValueSlope), F.regression.r2, F.regression.n);
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

function [b0, b1, seB1, pB1, tB1, df, r2] = weighted_linear_regression(x, y, w)
x = double(x(:));
y = double(y(:));
w = double(w(:));

ok = isfinite(x) & isfinite(y) & isfinite(w) & (w > 0);
x = x(ok);
y = y(ok);
w = w(ok);

n = numel(x);
b0 = NaN;
b1 = NaN;
seB1 = NaN;
pB1 = NaN;
tB1 = NaN;
df = NaN;
r2 = NaN;
if n < 2
    return;
end

X = [ones(n,1), x];
XtW = X' .* w';
XtWX = XtW * X;
XtWy = XtW * y;
if rcond(XtWX) <= eps
    return;
end
beta = XtWX \ XtWy;
b0 = beta(1);
b1 = beta(2);

yHat = X * beta;
res = y - yHat;
wrss = sum(w .* (res.^2));
yBarW = sum(w .* y) / sum(w);
wtss = sum(w .* ((y - yBarW).^2));
if isfinite(wtss) && wtss > 0
    r2 = 1 - wrss / wtss;
end

df = n - 2;
if df <= 0
    return;
end
mse = wrss / df;
covBeta = mse * inv(XtWX); %#ok<MINV>
seB1 = sqrt(max(covBeta(2,2), 0));
if isfinite(seB1) && seB1 > 0
    tB1 = b1 / seB1;
    xBeta = df / (df + tB1.^2);
    pB1 = betainc(xBeta, df/2, 0.5);
elseif isfinite(b1) && b1 == 0
    tB1 = 0;
    pB1 = 1;
end
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

function s = format_p_value(p)
if ~isfinite(p)
    s = 'NaN';
elseif p < 1e-3
    s = '<1e-3';
else
    s = sprintf('%.3f', p);
end
end
