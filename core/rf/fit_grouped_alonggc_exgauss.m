function F = fit_grouped_alonggc_exgauss(G, varargin)
% FIT_GROUPED_ALONGGC_EXGAUSS
% Fit exGauss_mod to each grouped along_GC trace.
%
% Uses model: exGauss_mod(p,t), p = [mu sigma alpha c d]
% Fits each group independently via lsqcurvefit.

p = inputParser;
p.addParameter('timeRef', 'center', @(x) ischar(x) || isstring(x)); % 'start'|'center'|'end'
p.addParameter('fitStartMs', -200, @(x) isnumeric(x) && isscalar(x));
p.addParameter('fitEndMs', 500, @(x) isnumeric(x) && isscalar(x));
p.addParameter('smoothW', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('minSigma', 5, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('maxSigma', 400, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('plotMode', 'overlay', @(x) ischar(x) || isstring(x)); % 'overlay'|'none'
p.addParameter('cmapName', 'parula', @(x) ischar(x) || isstring(x));
p.addParameter('lineWidth', 1.8, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('fitLineWidthFactor', 2.2, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('rawFadeToWhite', 0.65, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
p.addParameter('showT33Markers', true, @(x) islogical(x) && isscalar(x));
p.addParameter('t33MarkerSize', 7.5, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('sigThreshold', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x > 0));
p.addParameter('sigMinConsecutiveBins', 3, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('sigSearchStartMs', 0, @(x) isnumeric(x) && isscalar(x) && isfinite(x));
p.addParameter('sigUseAbs', true, @(x) islogical(x) && isscalar(x));
p.addParameter('verbose', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;

assert(isstruct(G) && isfield(G, 'groupMeanSigned') && isfield(G, 'timeWindows'), ...
    'G must be grouped output from build_grouped_alonggc_polygons_allbins.');
assert(exist('exGauss_mod', 'file') == 2, 'exGauss_mod not found on MATLAB path.');
assert(exist('lsqcurvefit', 'file') == 2, 'lsqcurvefit not found (Optimization Toolbox required).');

Y = double(G.groupMeanSigned);
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
t = t(:);

fitMask = isfinite(t) & (t >= opt.fitStartMs) & (t <= opt.fitEndMs);
assert(any(fitMask), 'No bins in fit range [%.1f, %.1f] ms.', opt.fitStartMs, opt.fitEndMs);

% Colors
cmapFun = str2func(char(opt.cmapName));
try
    C = cmapFun(nGroups);
catch
    C = parula(nGroups);
end

% Fit outputs
P = nan(nGroups, 5);
rmse = nan(nGroups,1);
t33Abs = nan(nGroups,1);
tPeakAbs = nan(nGroups,1);
peakAbs = nan(nGroups,1);
exitflag = nan(nGroups,1);
yFitAll = nan(nGroups, nBins);
status = strings(nGroups,1);
tSigOnset = nan(nGroups,1);
idxSigOnset = nan(nGroups,1);

% "Significant" threshold for onset detection.
sigThr = NaN;
if ~isempty(opt.sigThreshold)
    sigThr = double(opt.sigThreshold);
elseif isfield(G, 'calibration') && isstruct(G.calibration) && ...
        isfield(G.calibration, 'thresholdPreQ') && isfinite(G.calibration.thresholdPreQ) && G.calibration.thresholdPreQ > 0
    sigThr = double(G.calibration.thresholdPreQ);
end
minRun = round(opt.sigMinConsecutiveBins);

optsLSQ = optimoptions('lsqcurvefit', 'Display', 'off');
modfun = @(pp,tt) exGauss_mod(pp,tt);

for g = 1:nGroups
    y = Y(g,:).';
    ok = fitMask & isfinite(y);
    tfit = t(ok);
    yfit = y(ok);
    if numel(tfit) < 6
        status(g) = "too_few_points";
        continue;
    end

    if opt.smoothW > 1
        yfitUse = smoothdata(yfit, 'movmean', round(opt.smoothW));
    else
        yfitUse = yfit;
    end

    [~, iPeak] = max(abs(yfitUse));
    mu0 = tfit(iPeak);
    sigma0 = min(max(20, opt.minSigma), opt.maxSigma);
    alpha0 = 0.01;

    A = max(1e-3, 5 * max(abs(yfitUse)));
    sgn = sign(yfitUse(iPeak));
    if sgn == 0, sgn = 1; end
    c0 = 0.2 * A * sgn;
    d0 = 0;

    p0 = [mu0, sigma0, alpha0, c0, d0];
    lb = [min(tfit), opt.minSigma, 0, -A, -A];
    ub = [max(tfit), opt.maxSigma, 10,  A,  A];

    try
        [pHat, ~, ~, ef] = lsqcurvefit(modfun, p0, tfit, yfitUse, lb, ub, optsLSQ);
        yHatFit = modfun(pHat, tfit);
        yHatAll = modfun(pHat, t);
        P(g,:) = pHat;
        yFitAll(g,:) = yHatAll;
        rmse(g) = sqrt(mean((yHatFit - yfitUse).^2, 'omitnan'));
        exitflag(g) = ef;

        ya = abs(yHatFit);
        mx = max(ya, [], 'omitnan');
        peakAbs(g) = mx;
        if isfinite(mx) && mx > 0
            thr = 0.33 * mx;
            idx = find(ya >= thr, 1, 'first');
            if isempty(idx)
                t33Abs(g) = NaN;
            elseif idx == 1
                t33Abs(g) = tfit(1);
            else
                t33Abs(g) = interp1(ya(idx-1:idx), tfit(idx-1:idx), thr, 'linear');
            end
            [~, ip] = max(ya);
            tPeakAbs(g) = tfit(ip);
        end
        status(g) = "ok";
    catch
        status(g) = "fit_failed";
    end

    % First "significant onset" bin: first of >=minRun consecutive bins above threshold.
    if isfinite(sigThr) && sigThr > 0
        if opt.sigUseAbs
            ySig = abs(y);
        else
            ySig = y;
        end
        sigMask = isfinite(ySig) & (ySig > sigThr) & (t >= opt.sigSearchStartMs);
        idx0 = first_consecutive_true(sigMask, minRun);
        if isfinite(idx0)
            idxSigOnset(g) = idx0;
            tSigOnset(g) = t(idx0);
        end
    end
end

if opt.verbose
    fprintf('exGauss grouped fit: %d/%d groups fitted.\n', nnz(status=="ok"), nGroups);
    if isfinite(sigThr)
        fprintf(['Significance onset rule: first of >=%d consecutive bins with %s > %.6g ' ...
                 '(search from %.1f ms)\n'], ...
            minRun, ternary(opt.sigUseAbs, '|value|', 'value'), sigThr, opt.sigSearchStartMs);
    end
end

fig = [];
ax = [];
mode = lower(string(opt.plotMode));
if mode == "overlay"
    fig = figure('Color', 'w');
    ax = axes('Parent', fig); hold(ax, 'on');

    fitLW = max(opt.lineWidth * opt.fitLineWidthFactor, opt.lineWidth + 0.8);
    rawW = max(0.8, 0.9 * opt.lineWidth);

    for g = 1:nGroups
        lbl = sprintf('G%d data', g);
        if isfield(G, 'groupSummary') && istable(G.groupSummary) && height(G.groupSummary) >= g
            a0 = double(G.groupSummary.alongMin(g));
            a1 = double(G.groupSummary.alongMax(g));
            lbl = sprintf('G%d [%.2f, %.2f]', g, a0, a1);
        end
        rawCol = (1 - opt.rawFadeToWhite) * C(g,:) + opt.rawFadeToWhite * [1 1 1];
        plot(ax, t, Y(g,:), '-', 'Color', rawCol, 'LineWidth', rawW, 'DisplayName', lbl);
        if all(isfinite(yFitAll(g,:)))
            % Fit is emphasized: thicker and continuous.
            plot(ax, t, yFitAll(g,:), '-', 'Color', C(g,:), 'LineWidth', fitLW, ...
                'HandleVisibility', 'off');
        end
    end
    xline(ax, 0, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.0, 'DisplayName', 'onset');
    yline(ax, 0, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.0, 'DisplayName', 'zero');

    % Ensure y-limits include 0 so t33 markers can sit on the x-axis line.
    allVals = [Y(:); yFitAll(:); 0];
    allVals = allVals(isfinite(allVals));
    if ~isempty(allVals)
        yLo = min(allVals);
        yHi = max(allVals);
        if yLo == yHi
            pad = max(1e-3, 0.2 * abs(yLo) + 1e-3);
            yLo = yLo - pad;
            yHi = yHi + pad;
        else
            pad = 0.08 * (yHi - yLo);
            yLo = yLo - pad;
            yHi = yHi + pad;
        end
        yLo = min(yLo, 0);
        yHi = max(yHi, 0);
        ylim(ax, [yLo yHi]);
    end

    % Mark each group's t33 on the x-axis (y=0) with matching color.
    if opt.showT33Markers
        for g = 1:nGroups
            if isfinite(t33Abs(g))
                plot(ax, t33Abs(g), 0, 'o', ...
                    'MarkerSize', opt.t33MarkerSize, ...
                    'MarkerFaceColor', C(g,:), ...
                    'MarkerEdgeColor', [0.1 0.1 0.1], ...
                    'LineWidth', 0.8, ...
                    'HandleVisibility', 'off');
            end
        end
    end

    grid(ax, 'on'); box(ax, 'off');
    xlabel(ax, 'Time (ms)');
    ylabel(ax, 'Mean signed \Delta (T-D)');
    title(ax, 'Grouped traces with exGauss fits');
    legend(ax, 'Location', 'eastoutside');
end

summary = table((1:nGroups)', ...
    P(:,1), P(:,2), P(:,3), P(:,4), P(:,5), ...
    rmse, t33Abs, tPeakAbs, peakAbs, ...
    repmat(sigThr, nGroups, 1), repmat(minRun, nGroups, 1), idxSigOnset, tSigOnset, ...
    exitflag, status, ...
    'VariableNames', {'groupIdx','mu','sigma','alpha','c','d','rmse','t33Abs','tPeakAbs','peakAbs', ...
                      'sigThreshold','sigMinRunBins','idxSigOnset','tSigOnset','exitflag','status'});

F = struct();
F.summary = summary;
F.params = P;
F.t = t;
F.y = Y;
F.yFit = yFitAll;
F.colors = C;
F.fig = fig;
F.ax = ax;
F.sigThreshold = sigThr;
F.sigMinRunBins = minRun;
F.idxSigOnset = idxSigOnset;
F.tSigOnset = tSigOnset;

end

function idx0 = first_consecutive_true(mask, minRun)
idx0 = NaN;
mask = logical(mask(:));
n = numel(mask);
if n < minRun
    return;
end
for i = 1:(n - minRun + 1)
    if all(mask(i:i+minRun-1))
        idx0 = i;
        return;
    end
end
end

function y = ternary(cond, a, b)
if cond
    y = a;
else
    y = b;
end
end
