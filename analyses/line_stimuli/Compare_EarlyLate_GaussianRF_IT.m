function CMP = Compare_EarlyLate_GaussianRF_IT(Monkey, Opts)
% Compare early and late IT Gaussian RF estimates site-by-site.

if nargin < 1 || isempty(Monkey)
    Monkey = 1; % 1 = Nilson, 2 = Figaro
end
if nargin < 2 || isempty(Opts)
    Opts = struct();
end

Opts = normalize_compare_opts_local(Opts);

alpha = 0.05;
figNum = 105;

cfg = config();

if Monkey == 1
    monkeySuffix = "N";
elseif Monkey == 2
    monkeySuffix = "F";
else
    error('Compare_EarlyLate_GaussianRF_IT:InvalidMonkey', ...
        'Monkey must be 1 (Nilson) or 2 (Figaro).');
end

gaussPath = fullfile(cfg.matDir, sprintf('GaussianOccupancy_Tuning_IT_%s.mat', char(monkeySuffix)));
assert(exist(gaussPath, 'file') == 2, ...
    'Missing %s. Run GaussianOccupancy_Tuning_IT.m first.', gaussPath);

Sg = load(gaussPath, 'OUT');
assert(isfield(Sg, 'OUT') && isstruct(Sg.OUT), '%s must contain struct OUT.', gaussPath);
OUT = Sg.OUT;
assert(isfield(OUT, 'RFrange') && isfield(OUT, 'FitSpatialEarly') && isfield(OUT, 'FitSpatialLate'), ...
    'Gaussian OUT is missing required fields.');

FitEarly = OUT.FitSpatialEarly;
FitLate = OUT.FitSpatialLate;
RFrange = OUT.RFrange(:);
nSites = numel(RFrange);

gxEarly = [FitEarly.centerX].';
gyEarly = [FitEarly.centerY].';
gxLate = [FitLate.centerX].';
gyLate = [FitLate.centerY].';
sigmaEarly = [FitEarly.sigmaPx].';
sigmaLate = [FitLate.sigmaPx].';
pEarly = [FitEarly.pValueApprox].';
pLate = [FitLate.pValueApprox].';
veEarly = [FitEarly.r2TrainPct].';
veLate = [FitLate.r2TrainPct].';

validEarly = isfinite(gxEarly) & isfinite(gyEarly) & isfinite(sigmaEarly) & ...
    (sigmaEarly > 0) & (sigmaEarly <= Opts.MaxGaussianSigmaPx) & ...
    isfinite(veEarly) & (veEarly >= Opts.MinSpatialVE);
validLate = isfinite(gxLate) & isfinite(gyLate) & isfinite(sigmaLate) & ...
    (sigmaLate > 0) & (sigmaLate <= Opts.MaxGaussianSigmaPx) & ...
    isfinite(veLate) & (veLate >= Opts.MinSpatialVE);
validBoth = validEarly & validLate;
validUnion = validEarly | validLate;

sigEarly = validEarly & isfinite(pEarly) & (pEarly < alpha);
sigLate = validLate & isfinite(pLate) & (pLate < alpha);
sigBoth = sigEarly & sigLate;
sigEarlyOnly = sigEarly & ~sigLate;
sigLateOnly = sigLate & ~sigEarly;
sigNeither = validUnion & ~sigEarly & ~sigLate;

if Opts.RequireBothSpatialSig
    primaryMask = sigBoth;
    primaryLabel = "both spatially significant";
else
    primaryMask = validBoth;
    primaryLabel = "all valid paired fits";
end

dxAll = gxLate - gxEarly;
dyAll = gyLate - gyEarly;
distAll = hypot(dxAll, dyAll);
meanSigmaAll = 0.5 * (sigmaEarly + sigmaLate);
normDistAll = distAll ./ meanSigmaAll;
log2RatioAll = log2(sigmaLate ./ sigmaEarly);
[radialShiftAll, tangentialShiftAll] = decompose_shift_local(gxEarly, gyEarly, gxLate, gyLate);

sigmaListPx = extract_sigma_list_local(OUT, sigmaEarly, sigmaLate);
sigmaIdxEarly = map_sigma_to_idx_local(sigmaEarly, sigmaListPx);
sigmaIdxLate = map_sigma_to_idx_local(sigmaLate, sigmaListPx);
sigmaStepDiffAll = sigmaIdxLate - sigmaIdxEarly;

centerStableMask = primaryMask & isfinite(normDistAll) & (normDistAll <= Opts.ShiftNormThresh);
centerShiftedMask = primaryMask & isfinite(normDistAll) & (normDistAll > Opts.ShiftNormThresh);
sizeStableMask = primaryMask & isfinite(sigmaStepDiffAll) & ...
    (abs(sigmaStepDiffAll) < Opts.MinSigmaStepForSizeChange);
sizeExpandedMask = primaryMask & isfinite(sigmaStepDiffAll) & ...
    (sigmaStepDiffAll >= Opts.MinSigmaStepForSizeChange);
sizeContractedMask = primaryMask & isfinite(sigmaStepDiffAll) & ...
    (sigmaStepDiffAll <= -Opts.MinSigmaStepForSizeChange);
sizeSameSigmaMask = primaryMask & (sigmaStepDiffAll == 0);
sizeOneStepMask = primaryMask & (abs(sigmaStepDiffAll) == 1);
sizeMultiStepMask = primaryMask & (abs(sigmaStepDiffAll) >= 2);

allValidSummary = summarize_paired_window(gxEarly, gyEarly, sigmaEarly, gxLate, gyLate, sigmaLate, ...
    dxAll, dyAll, distAll, normDistAll, radialShiftAll, tangentialShiftAll, ...
    log2RatioAll, sigmaStepDiffAll, validBoth, Opts.nShuffle, Opts.nEccBins);
primarySummary = summarize_paired_window(gxEarly, gyEarly, sigmaEarly, gxLate, gyLate, sigmaLate, ...
    dxAll, dyAll, distAll, normDistAll, radialShiftAll, tangentialShiftAll, ...
    log2RatioAll, sigmaStepDiffAll, primaryMask, Opts.nShuffle, Opts.nEccBins);

CMP = struct();
CMP.monkeySuffix = monkeySuffix;
CMP.alpha = alpha;
CMP.gaussPath = gaussPath;
CMP.filterOptions = Opts;
CMP.primaryLabel = primaryLabel;
CMP.globalSiteInR = RFrange;
CMP.validEarly = validEarly;
CMP.validLate = validLate;
CMP.validBoth = validBoth;
CMP.validUnion = validUnion;
CMP.sigEarly = sigEarly;
CMP.sigLate = sigLate;
CMP.sigBoth = sigBoth;
CMP.sigEarlyOnly = sigEarlyOnly;
CMP.sigLateOnly = sigLateOnly;
CMP.sigNeither = sigNeither;
CMP.primaryMask = primaryMask;
CMP.centerStableMask = centerStableMask;
CMP.centerShiftedMask = centerShiftedMask;
CMP.sizeStableMask = sizeStableMask;
CMP.sizeExpandedMask = sizeExpandedMask;
CMP.sizeContractedMask = sizeContractedMask;
CMP.sizeSameSigmaMask = sizeSameSigmaMask;
CMP.sizeOneStepMask = sizeOneStepMask;
CMP.sizeMultiStepMask = sizeMultiStepMask;
CMP.sigmaListPx = sigmaListPx;
CMP.dx = dxAll;
CMP.dy = dyAll;
CMP.dist = distAll;
CMP.normDist = normDistAll;
CMP.radialShift = radialShiftAll;
CMP.tangentialShift = tangentialShiftAll;
CMP.log2Ratio = log2RatioAll;
CMP.sigmaStepDiff = sigmaStepDiffAll;
CMP.allValid = allValidSummary;
CMP.primary = primarySummary;

fprintf('IT Gaussian RF early-vs-late comparison (%s)\n', char(monkeySuffix));
fprintf('Filters: %s\n', char(describe_filter_opts_local(Opts)));
fprintf('Finite Gaussian fits | early=%d | late=%d | both=%d | union=%d | excluded=%d\n', ...
    nnz(validEarly), nnz(validLate), nnz(validBoth), nnz(validUnion), nSites - nnz(validUnion));
fprintf('Spatially significant at p < %.2f | early=%d | late=%d | both=%d | early-only=%d | late-only=%d | neither=%d\n', ...
    alpha, nnz(sigEarly), nnz(sigLate), nnz(sigBoth), nnz(sigEarlyOnly), nnz(sigLateOnly), nnz(sigNeither));
print_paired_summary('All valid paired fits', allValidSummary);
print_paired_summary(sprintf('Primary cohort (%s)', char(primaryLabel)), primarySummary);
fprintf('Primary cohort center classification | stable=%d | shifted (norm dist > %.2f sigma)= %d\n', ...
    nnz(centerStableMask), Opts.ShiftNormThresh, nnz(centerShiftedMask));
fprintf('Primary cohort size classification | stable=%d | expanded=%d | contracted=%d | same sigma=%d | one step=%d | multi-step=%d\n', ...
    nnz(sizeStableMask), nnz(sizeExpandedMask), nnz(sizeContractedMask), ...
    nnz(sizeSameSigmaMask), nnz(sizeOneStepMask), nnz(sizeMultiStepMask));

figTitle = sprintf('IT Gaussian RF early vs late (%s)', char(monkeySuffix));
figure(figNum); clf;
set(gcf, 'Color', 'w', 'Name', figTitle, 'NumberTitle', 'off', ...
    'Tag', 'IT_Gaussian_RF_early_vs_late');

useTiled = exist('tiledlayout', 'file') == 2;
if useTiled
    tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
end

displayMask = primaryMask;
displayLabel = primaryLabel;
displaySummary = primarySummary;
displayShiftMask = centerShiftedMask;
if ~any(displayMask) && any(validBoth)
    displayMask = validBoth;
    displayLabel = "all valid paired fits";
    displaySummary = allValidSummary;
    displayShiftMask = false(size(displayMask));
end

[xLim, yLim] = choose_center_limits_local(gxEarly(displayMask), gyEarly(displayMask), ...
    gxLate(displayMask), gyLate(displayMask));

if useTiled, nexttile; else, subplot(2, 3, 1); end
draw_center_arrow_panel(gxEarly, gyEarly, gxLate, gyLate, displayMask, displayShiftMask, ...
    xLim, yLim, displaySummary, displayLabel);

if useTiled, nexttile; else, subplot(2, 3, 2); end
draw_norm_shift_hist_panel(normDistAll(validBoth), normDistAll(primaryMask), Opts.ShiftNormThresh, primarySummary);

if useTiled, nexttile; else, subplot(2, 3, 3); end
draw_shuffle_panel(primarySummary, displaySummary, primaryLabel, displayLabel);

if useTiled, nexttile; else, subplot(2, 3, 4); end
draw_sigma_scatter_panel(sigmaEarly, sigmaLate, validBoth, primaryMask, primarySummary);

if useTiled, nexttile; else, subplot(2, 3, 5); end
draw_log_ratio_panel(log2RatioAll(validBoth), log2RatioAll(primaryMask), sigmaStepDiffAll(primaryMask), primarySummary);

if useTiled, nexttile; else, subplot(2, 3, 6); end
draw_summary_bar_panel(sigEarlyOnly, sigBoth, sigLateOnly, sigNeither, ...
    centerStableMask, centerShiftedMask, sizeStableMask, sizeExpandedMask, sizeContractedMask, ...
    primaryLabel);

sgText = sprintf('%s | %s | %s', figTitle, char(primaryLabel), char(describe_filter_opts_local(Opts)));
if exist('sgtitle', 'file') == 2
    sgtitle(sgText);
else
    annotation('textbox', [0.08 0.955 0.84 0.04], ...
        'String', sgText, 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end
end

function W = summarize_paired_window(gxEarly, gyEarly, sigmaEarly, gxLate, gyLate, sigmaLate, ...
    dxAll, dyAll, distAll, normDistAll, radialShiftAll, tangentialShiftAll, ...
    log2RatioAll, sigmaStepDiffAll, mask, nShuffle, nEccBins)
W = struct();
W.mask = mask(:);
W.gxEarly = gxEarly(mask);
W.gyEarly = gyEarly(mask);
W.gxLate = gxLate(mask);
W.gyLate = gyLate(mask);
W.sigmaEarly = sigmaEarly(mask);
W.sigmaLate = sigmaLate(mask);
W.dx = dxAll(mask);
W.dy = dyAll(mask);
W.dist = distAll(mask);
W.normDist = normDistAll(mask);
W.radialShift = radialShiftAll(mask);
W.tangentialShift = tangentialShiftAll(mask);
W.log2Ratio = log2RatioAll(mask);
W.sigmaStepDiff = sigmaStepDiffAll(mask);
W.nPaired = nnz(mask);

if W.nPaired >= 2
    W.rX = corr(W.gxEarly, W.gxLate, 'Rows', 'complete');
    W.rY = corr(W.gyEarly, W.gyLate, 'Rows', 'complete');
    W.rEcc = corr(hypot(W.gxEarly, W.gyEarly), hypot(W.gxLate, W.gyLate), 'Rows', 'complete');
    W.rSigmaLog = corr(log10(W.sigmaEarly), log10(W.sigmaLate), 'Rows', 'complete');
    W.rhoSigma = corr(W.sigmaEarly, W.sigmaLate, 'Type', 'Spearman', 'Rows', 'complete');
else
    W.rX = NaN;
    W.rY = NaN;
    W.rEcc = NaN;
    W.rSigmaLog = NaN;
    W.rhoSigma = NaN;
end

W.medianDist = median(W.dist, 'omitnan');
W.medianNormDist = median(W.normDist, 'omitnan');
W.meanNormDist = mean(W.normDist, 'omitnan');
W.meanDx = mean(W.dx, 'omitnan');
W.meanDy = mean(W.dy, 'omitnan');
W.medianRadialShift = median(W.radialShift, 'omitnan');
W.meanRadialShift = mean(W.radialShift, 'omitnan');
W.medianTangentialShift = median(W.tangentialShift, 'omitnan');
W.meanTangentialShift = mean(W.tangentialShift, 'omitnan');
W.medianLog2Ratio = median(W.log2Ratio, 'omitnan');
W.meanLog2Ratio = mean(W.log2Ratio, 'omitnan');
W.medianStepDiff = median(W.sigmaStepDiff, 'omitnan');

if W.nPaired >= 2 && nShuffle > 0
    eccMid = hypot(0.5 * (W.gxEarly + W.gxLate), 0.5 * (W.gyEarly + W.gyLate));
    [W.nullMedianDist, W.binIdx] = build_distance_shuffle_null( ...
        W.gxEarly, W.gyEarly, eccMid, W.gxLate, W.gyLate, nShuffle, nEccBins);
    W.shuffleP = (1 + nnz(W.nullMedianDist <= W.medianDist)) / (numel(W.nullMedianDist) + 1);
else
    W.nullMedianDist = NaN;
    W.binIdx = ones(W.nPaired, 1);
    W.shuffleP = NaN;
end
end

function print_paired_summary(label, W)
fprintf(['%s | paired=%d | median dist=%.2f px | median norm dist=%.2f sigma | ' ...
    'median log2(late/early)=%.2f | shuffle p=%.4f | rX=%.3f | rY=%.3f | rSigma(log)=%.3f\n'], ...
    label, W.nPaired, W.medianDist, W.medianNormDist, W.medianLog2Ratio, ...
    W.shuffleP, W.rX, W.rY, W.rSigmaLog);
fprintf('%s | mean dx=%.2f | mean dy=%.2f | median radial shift=%.2f px | median tangential shift=%.2f px\n', ...
    label, W.meanDx, W.meanDy, W.medianRadialShift, W.medianTangentialShift);
end

function [radialShift, tangentialShift] = decompose_shift_local(x1, y1, x2, y2)
dx = x2 - x1;
dy = y2 - y1;
refX = 0.5 * (x1 + x2);
refY = 0.5 * (y1 + y2);
refNorm = hypot(refX, refY);

radialShift = nan(size(dx));
tangentialShift = nan(size(dx));
good = isfinite(dx) & isfinite(dy) & isfinite(refNorm) & (refNorm > 0);
if any(good)
    ux = refX(good) ./ refNorm(good);
    uy = refY(good) ./ refNorm(good);
    radialShift(good) = dx(good) .* ux + dy(good) .* uy;
    tangentialShift(good) = -dx(good) .* uy + dy(good) .* ux;
end
end

function sigmaListPx = extract_sigma_list_local(OUT, sigmaEarly, sigmaLate)
sigmaListPx = [];
if isfield(OUT, 'Library') && isstruct(OUT.Library)
    if isfield(OUT.Library, 'sigmaListPx') && ~isempty(OUT.Library.sigmaListPx)
        sigmaListPx = double(OUT.Library.sigmaListPx(:));
    elseif isfield(OUT.Library, 'sigmaPx') && ~isempty(OUT.Library.sigmaPx)
        sigmaListPx = unique(double(OUT.Library.sigmaPx(:)));
    end
end
if isempty(sigmaListPx)
    sigmaListPx = unique([sigmaEarly(isfinite(sigmaEarly)); sigmaLate(isfinite(sigmaLate))]);
end
sigmaListPx = sort(sigmaListPx(:));
end

function idx = map_sigma_to_idx_local(vals, sigmaListPx)
idx = nan(size(vals));
if isempty(sigmaListPx)
    return;
end
for i = 1:numel(vals)
    if ~isfinite(vals(i))
        continue;
    end
    [d, k] = min(abs(sigmaListPx - vals(i))); %#ok<ASGLU>
    if isfinite(d) && d <= max(1e-6, abs(vals(i)) * 1e-6)
        idx(i) = k;
    end
end
end

function [nullMedianDist, binIdx] = build_distance_shuffle_null(xEarly, yEarly, eccRef, xLate, yLate, nShuffle, nEccBins)
n = numel(xEarly);
nullMedianDist = nan(nShuffle, 1);
binIdx = assign_ecc_bins_local(eccRef, nEccBins);

if n < 2
    return;
end

for iShuf = 1:nShuffle
    perm = zeros(n, 1);
    for b = unique(binIdx(:)).'
        idx = find(binIdx == b);
        if numel(idx) == 1
            perm(idx) = idx;
        else
            perm(idx) = idx(randperm(numel(idx)));
        end
    end
    d = hypot(xLate(perm) - xEarly, yLate(perm) - yEarly);
    nullMedianDist(iShuf) = median(d, 'omitnan');
end
end

function binIdx = assign_ecc_bins_local(ecc, nEccBins)
ecc = ecc(:);
binIdx = ones(size(ecc));
if numel(ecc) < max(4, nEccBins)
    return;
end

edges = prctile(ecc, linspace(0, 100, nEccBins + 1));
edges(1) = -inf;
edges(end) = inf;
edges = make_strictly_increasing_local(edges);

tmp = discretize(ecc, edges);
if any(isnan(tmp))
    tmp(isnan(tmp)) = 1;
end
binIdx = tmp(:);
end

function edges = make_strictly_increasing_local(edges)
for i = 2:numel(edges)
    if ~(edges(i) > edges(i-1))
        edges(i) = edges(i-1) + max(1e-6, abs(edges(i-1)) * 1e-6);
    end
end
end

function [xLim, yLim] = choose_center_limits_local(x1, y1, x2, y2)
x = [x1(:); x2(:)];
y = [y1(:); y2(:)];
x = x(isfinite(x));
y = y(isfinite(y));

if isempty(x) || isempty(y)
    xLim = [-150 150];
    yLim = [-150 150];
    return;
end

pad = 15;
xHalf = prctile(abs(x), 98);
yHalf = prctile(abs(y), 98);
if ~isfinite(xHalf) || xHalf <= 0
    xHalf = max(abs(x));
end
if ~isfinite(yHalf) || yHalf <= 0
    yHalf = max(abs(y));
end

xHalf = max(100, ceil((xHalf + pad) / 25) * 25);
yHalf = max(100, ceil((yHalf + pad) / 25) * 25);

xLim = [-xHalf xHalf];
yLim = [-yHalf yHalf];
end

function draw_center_arrow_panel(gxEarly, gyEarly, gxLate, gyLate, mask, shiftedMask, xLim, yLim, W, label)
if ~any(mask)
    text(0.5, 0.5, 'No paired sites for center comparison', ...
        'HorizontalAlignment', 'center');
    axis off;
    return;
end

hold on;
idx = find(mask);
isShifted = shiftedMask(idx);

for i = 1:numel(idx)
    k = idx(i);
    if isShifted(i)
        lineColor = [0.85 0.10 0.10];
        lineWidth = 1.2;
    else
        lineColor = [0.70 0.70 0.70];
        lineWidth = 0.8;
    end
    plot([gxEarly(k) gxLate(k)], [gyEarly(k) gyLate(k)], '-', ...
        'Color', lineColor, 'LineWidth', lineWidth);
end

hEarly = plot(gxEarly(mask), gyEarly(mask), 'o', ...
    'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', ...
    'MarkerSize', 5, 'LineStyle', 'none');
hLateStable = plot(gxLate(mask & ~shiftedMask), gyLate(mask & ~shiftedMask), 'o', ...
    'MarkerFaceColor', [0.55 0.55 0.55], 'MarkerEdgeColor', 'none', ...
    'MarkerSize', 5, 'LineStyle', 'none');
hLateShift = plot(gxLate(mask & shiftedMask), gyLate(mask & shiftedMask), 'o', ...
    'MarkerFaceColor', [0.85 0.10 0.10], 'MarkerEdgeColor', 'k', ...
    'MarkerSize', 6, 'LineStyle', 'none');

plot([0 0], yLim, 'k:');
plot(xLim, [0 0], 'k:');
xlim(xLim);
ylim(yLim);
axis equal;
grid on;
xlabel('x center (px)');
ylabel('y center (px)');
title(sprintf('Center shifts | %s | paired=%d', char(label), W.nPaired));
legend([hEarly, hLateStable, hLateShift], ...
    {'Early center', 'Late center', 'Late center (shifted)'}, ...
    'Location', 'southeast');

txt = sprintf('median dist = %.1f px | median norm = %.2f sigma', ...
    W.medianDist, W.medianNormDist);
text(xLim(1) + 0.03 * diff(xLim), yLim(2) - 0.06 * diff(yLim), txt, 'FontSize', 9);
end

function draw_norm_shift_hist_panel(normDistValid, normDistPrimary, shiftThresh, W)
if isempty(normDistValid)
    text(0.5, 0.5, 'No paired sites for normalized shift histogram', ...
        'HorizontalAlignment', 'center');
    axis off;
    return;
end

valsAll = normDistValid(isfinite(normDistValid) & (normDistValid >= 0));
if isempty(valsAll)
    valsAll = 0;
end
hi = max([4; shiftThresh * 1.2; max(valsAll); prctile(valsAll, 98) * 1.1]);
hi = max(1, ceil(hi / 0.25) * 0.25);
edges = linspace(0, hi, 26);

hAll = histogram(normDistValid, edges, 'FaceColor', [0.82 0.82 0.82], 'EdgeColor', 'none');
hold on;
if ~isempty(normDistPrimary)
    hPrimary = histogram(normDistPrimary, edges, 'FaceColor', [0.85 0.10 0.10], 'EdgeColor', 'none', ...
        'FaceAlpha', 0.60);
else
    hPrimary = plot(nan, nan, 's', 'MarkerFaceColor', [0.85 0.10 0.10], ...
        'MarkerEdgeColor', 'k', 'MarkerSize', 6, 'LineStyle', 'none');
end
hThresh = xline(shiftThresh, 'k--', 'LineWidth', 1.2);
grid on;
xlabel('Center shift / mean sigma');
ylabel('N sites');
title('Normalized center shift');
txt = sprintf('primary median = %.2f sigma | mean = %.2f sigma', ...
    W.medianNormDist, W.meanNormDist);
text(min(xlim) + 0.03 * diff(xlim), max(ylim) - 0.08 * diff(ylim), txt, 'FontSize', 9);
legend([hAll, hPrimary, hThresh], {'All valid paired fits', 'Primary cohort', 'Shift threshold'}, ...
    'Location', 'northeast');
end

function draw_shuffle_panel(Wprimary, Wdisplay, primaryLabel, displayLabel)
if isfinite(Wprimary.shuffleP) && numel(Wprimary.nullMedianDist) > 1
    vals = Wprimary.nullMedianDist(isfinite(Wprimary.nullMedianDist));
    histogram(vals, 35, 'FaceColor', [0.82 0.82 0.82], 'EdgeColor', 'none');
    hold on;
    yl = ylim;
    plot([Wprimary.medianDist Wprimary.medianDist], yl, 'r-', 'LineWidth', 1.6);
    ylim(yl);
    grid on;
    xlabel('Median paired center distance (px)');
    ylabel('N shuffles');
    title(sprintf('Shuffle null | %s', char(primaryLabel)));
    txt = sprintf('obs = %.1f px | shuffle p = %.3g', Wprimary.medianDist, Wprimary.shuffleP);
    text(min(xlim) + 0.03 * diff(xlim), yl(2) - 0.08 * diff(yl), txt, 'FontSize', 9);
    return;
end

if isfinite(Wdisplay.shuffleP) && numel(Wdisplay.nullMedianDist) > 1
    vals = Wdisplay.nullMedianDist(isfinite(Wdisplay.nullMedianDist));
    histogram(vals, 35, 'FaceColor', [0.82 0.82 0.82], 'EdgeColor', 'none');
    hold on;
    yl = ylim;
    plot([Wdisplay.medianDist Wdisplay.medianDist], yl, 'r-', 'LineWidth', 1.6);
    ylim(yl);
    grid on;
    xlabel('Median paired center distance (px)');
    ylabel('N shuffles');
    title(sprintf('Shuffle null | %s', char(displayLabel)));
    txt = sprintf('obs = %.1f px | shuffle p = %.3g', Wdisplay.medianDist, Wdisplay.shuffleP);
    text(min(xlim) + 0.03 * diff(xlim), yl(2) - 0.08 * diff(yl), txt, 'FontSize', 9);
else
    text(0.5, 0.5, 'Not enough paired sites for shuffle test', ...
        'HorizontalAlignment', 'center');
    axis off;
end
end

function draw_sigma_scatter_panel(sigmaEarly, sigmaLate, validBoth, primaryMask, W)
valsX = sigmaEarly(validBoth);
valsY = sigmaLate(validBoth);
if isempty(valsX)
    text(0.5, 0.5, 'No paired sites for sigma comparison', ...
        'HorizontalAlignment', 'center');
    axis off;
    return;
end

lim = choose_positive_limits_local([valsX; valsY], 0.99);
plot(valsX, valsY, 'o', 'MarkerFaceColor', [0.78 0.78 0.78], ...
    'MarkerEdgeColor', 'none', 'MarkerSize', 5, 'LineStyle', 'none');
hold on;
plot(sigmaEarly(primaryMask), sigmaLate(primaryMask), 'o', ...
    'MarkerFaceColor', [0.85 0.10 0.10], 'MarkerEdgeColor', 'k', ...
    'MarkerSize', 6, 'LineStyle', 'none');
plot([lim(1) lim(2)], [lim(1) lim(2)], 'k--', 'LineWidth', 1.2);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlim(lim);
ylim(lim);
axis square;
grid on;
xlabel('Early sigma (px)');
ylabel('Late sigma (px)');
title(sprintf('Sigma comparison | paired=%d', nnz(validBoth)));
txt = sprintf('\\rho_s = %.2f | r_{log} = %.2f', W.rhoSigma, W.rSigmaLog);
text(exp(log(lim(1)) + 0.04 * diff(log(lim))), ...
    exp(log(lim(2)) - 0.08 * diff(log(lim))), ...
    txt, 'FontSize', 9);
legend({'All valid paired fits', 'Primary cohort', 'y = x'}, 'Location', 'southeast');
end

function draw_log_ratio_panel(log2RatioValid, log2RatioPrimary, sigmaStepDiffPrimary, W)
if isempty(log2RatioValid)
    text(0.5, 0.5, 'No paired sites for size-change histogram', ...
        'HorizontalAlignment', 'center');
    axis off;
    return;
end

allVals = log2RatioValid(isfinite(log2RatioValid));
if isempty(allVals)
    allVals = 0;
end
lim = max(1, ceil(max(abs([allVals(:); prctile(allVals, [2 98]).']))));
edges = linspace(-lim, lim, 25);

hAll = histogram(log2RatioValid, edges, 'FaceColor', [0.82 0.82 0.82], 'EdgeColor', 'none');
hold on;
if ~isempty(log2RatioPrimary)
    hPrimary = histogram(log2RatioPrimary, edges, 'FaceColor', [0.85 0.10 0.10], 'EdgeColor', 'none', ...
        'FaceAlpha', 0.60);
else
    hPrimary = plot(nan, nan, 's', 'MarkerFaceColor', [0.85 0.10 0.10], ...
        'MarkerEdgeColor', 'k', 'MarkerSize', 6, 'LineStyle', 'none');
end
hZero = xline(0, 'k-');
grid on;
xlabel('log_2(late sigma / early sigma)');
ylabel('N sites');
title('Size change');
sameSigma = nnz(sigmaStepDiffPrimary == 0);
oneStep = nnz(abs(sigmaStepDiffPrimary) == 1);
multiStep = nnz(abs(sigmaStepDiffPrimary) >= 2);
txt = sprintf('median = %.2f | same=%d | 1-step=%d | multi=%d', ...
    W.medianLog2Ratio, sameSigma, oneStep, multiStep);
text(min(xlim) + 0.03 * diff(xlim), max(ylim) - 0.08 * diff(ylim), txt, 'FontSize', 9);
legend([hAll, hPrimary, hZero], {'All valid paired fits', 'Primary cohort', 'No size change'}, ...
    'Location', 'northwest');
end

function draw_summary_bar_panel(sigEarlyOnly, sigBoth, sigLateOnly, sigNeither, ...
    centerStableMask, centerShiftedMask, sizeStableMask, sizeExpandedMask, sizeContractedMask, primaryLabel)
vals = [nnz(sigEarlyOnly), nnz(sigBoth), nnz(sigLateOnly), nnz(sigNeither), ...
    nnz(centerStableMask), nnz(centerShiftedMask), ...
    nnz(sizeStableMask), nnz(sizeExpandedMask), nnz(sizeContractedMask)];
pos = [1 2 3 4 6 7 9 10 11];
cols = [ ...
    0.25 0.45 0.85; ...
    0.65 0.25 0.75; ...
    0.90 0.45 0.20; ...
    0.72 0.72 0.72; ...
    0.55 0.55 0.55; ...
    0.85 0.10 0.10; ...
    0.55 0.55 0.55; ...
    0.20 0.60 0.25; ...
    0.85 0.35 0.20];

hold on;
for i = 1:numel(vals)
    bar(pos(i), vals(i), 0.75, 'FaceColor', cols(i,:), 'EdgeColor', 'none');
end
set(gca, 'XTick', pos, ...
    'XTickLabel', {'E only', 'Both sig', 'L only', 'Neither', ...
                   'Center stable', 'Center shift', ...
                   'Size stable', 'Expand', 'Contract'});
xtickangle(30);
grid on;
ylabel('N sites');
title(sprintf('State summary | primary: %s', char(primaryLabel)));
end

function lim = choose_positive_limits_local(vals, pctKeep)
vals = vals(isfinite(vals) & vals > 0);
if isempty(vals)
    lim = [1 10];
    return;
end
lo = prctile(vals, 1);
hi = prctile(vals, pctKeep * 100);
lo = max(lo, min(vals(vals > 0)));
hi = max(hi, lo * 1.5);
lim = [10^(floor(log10(lo))), 10^(ceil(log10(hi)))];
if lim(1) <= 0 || ~all(isfinite(lim)) || lim(1) >= lim(2)
    lim = [min(vals) max(vals)];
end
end

function Opts = normalize_compare_opts_local(Opts)
defaults = struct();
defaults.MaxGaussianSigmaPx = inf;
defaults.MinSpatialVE = -inf;
defaults.RequireBothSpatialSig = true;
defaults.ShiftNormThresh = 1.0;
defaults.MinSigmaStepForSizeChange = 1;
defaults.nShuffle = 5000;
defaults.nEccBins = 4;

fn = fieldnames(defaults);
for i = 1:numel(fn)
    if ~isfield(Opts, fn{i}) || isempty(Opts.(fn{i}))
        Opts.(fn{i}) = defaults.(fn{i});
    end
end
end

function txt = describe_filter_opts_local(Opts)
parts = {};
if isfinite(Opts.MaxGaussianSigmaPx)
    parts{end+1} = sprintf('Gaussian sigma <= %.1f', Opts.MaxGaussianSigmaPx); %#ok<AGROW>
end
if isfinite(Opts.MinSpatialVE)
    parts{end+1} = sprintf('VE >= %.1f%%', Opts.MinSpatialVE); %#ok<AGROW>
end
if Opts.RequireBothSpatialSig
    parts{end+1} = 'primary cohort = p < 0.05 in both windows'; %#ok<AGROW>
else
    parts{end+1} = 'primary cohort = all valid paired fits'; %#ok<AGROW>
end
parts{end+1} = sprintf('shift threshold = %.2f sigma', Opts.ShiftNormThresh); %#ok<AGROW>
parts{end+1} = sprintf('size threshold = %d sigma-step', Opts.MinSigmaStepForSizeChange); %#ok<AGROW>
txt = strjoin(parts, ', ');
end
