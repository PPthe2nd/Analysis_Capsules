function CMP = Compare_RFcenter_IT_Gaussian_vs_RFm(Monkey, UseResponsiveSubset, Opts)
% Compare RF.m center estimates to the best Gaussian center in IT.

if nargin < 1 || isempty(Monkey)
    Monkey = 1; % 1 = Nilson, 2 = Figaro
end
if nargin < 2 || isempty(UseResponsiveSubset)
    UseResponsiveSubset = true;
end
if nargin < 3 || isempty(Opts)
    Opts = struct();
end

Opts = normalize_compare_opts_local(Opts);

alpha = 0.05;
nShuffle = 5000;
nEccBins = 4;
figNum = 104;

cfg = config();

if Monkey == 1
    monkeySuffix = "N";
    rfFile = resolve_rf_file_local(cfg, 'THINGS_RF1s_N.mat', 'Mr Nilson');
elseif Monkey == 2
    monkeySuffix = "F";
    rfFile = resolve_rf_file_local(cfg, 'THINGS_RF1s_F.mat', 'Figaro');
else
    error('Compare_RFcenter_IT_Gaussian_vs_RFm:InvalidMonkey', ...
        'Monkey must be 1 (Nilson) or 2 (Figaro).');
end

gaussPath = fullfile(cfg.matDir, sprintf('GaussianOccupancy_Tuning_IT_%s.mat', char(monkeySuffix)));
colorTunePath = fullfile(cfg.matDir, sprintf('ColorTune_balanced_IT_%s.mat', char(monkeySuffix)));

assert(exist(rfFile, 'file') == 2, 'Missing RF file: %s', rfFile);
assert(exist(gaussPath, 'file') == 2, ...
    'Missing %s. Run GaussianOccupancy_Tuning_IT.m first.', gaussPath);

Srf = load(rfFile, 'all_centrex', 'all_centrey', 'all_szx', 'all_szy', ...
    'all_test_corr', 'all_test_corr_noisecorr');
Sg = load(gaussPath, 'OUT');
assert(isfield(Sg, 'OUT') && isstruct(Sg.OUT), '%s must contain struct OUT.', gaussPath);
OUT = Sg.OUT;
assert(isfield(OUT, 'RFrange') && isfield(OUT, 'FitSpatialEarly') && isfield(OUT, 'FitSpatialLate'), ...
    'Gaussian OUT is missing required fields.');

rfXAll = double(Srf.all_centrex(:));
rfYAll = double(Srf.all_centrey(:));
rfSizeEqAll = sqrt(double(Srf.all_szx(:)) .* double(Srf.all_szy(:)));
rfGlobal = OUT.RFrange(:);

assert(max(rfGlobal) <= numel(rfXAll), ...
    'RF indices exceed the RF file length.');

rfX = rfXAll(rfGlobal);
rfY = rfYAll(rfGlobal);
rfSizeEq = rfSizeEqAll(rfGlobal);
rfEcc = hypot(rfX, rfY);
rfValid = isfinite(rfX) & isfinite(rfY) & abs(rfX) <= 200 & abs(rfY) <= 200;

[rfTestCorr, corrLabel] = extract_rf_confidence_local(Srf, rfGlobal, Opts.RFcorrField);
sizeMask = isfinite(rfSizeEq) & (rfSizeEq > 0) & (rfSizeEq <= Opts.MaxRFsizeEq);
hasFiniteRFcorr = ~isempty(rfTestCorr) && any(isfinite(rfTestCorr));
if isempty(rfTestCorr) || ~isfinite(Opts.MinRFtestCorr) || ~hasFiniteRFcorr
    confMask = true(size(rfGlobal));
else
    confMask = isfinite(rfTestCorr) & (rfTestCorr >= Opts.MinRFtestCorr);
end

responsiveMask = true(size(rfGlobal));
subsetLabel = "all valid RF centers";
if UseResponsiveSubset
    if exist(colorTunePath, 'file') == 2
        Sct = load(colorTunePath, 'ColorTune');
        if isfield(Sct, 'ColorTune') && isstruct(Sct.ColorTune) && ...
                isfield(Sct.ColorTune, 'keepSites') && isfield(Sct.ColorTune, 'RFrange')
            globalKeep = Sct.ColorTune.RFrange(Sct.ColorTune.keepSites(:));
            responsiveMask = ismember(rfGlobal, globalKeep);
            if isfield(Sct.ColorTune, 'thr') && isfinite(Sct.ColorTune.thr)
                subsetLabel = sprintf('responsive IT sites (bestSNR > %.2f)', Sct.ColorTune.thr);
            else
                subsetLabel = 'responsive IT sites';
            end
            fprintf('Using responsive subset from %s: %d / %d sites\n', ...
                colorTunePath, nnz(responsiveMask), numel(responsiveMask));
        else
            fprintf('ColorTune file exists but is missing keepSites/RFrange. Using all valid RF centers.\n');
        end
    else
        fprintf('ColorTune file not found (%s). Using all valid RF centers.\n', colorTunePath);
    end
end

baseMask = rfValid & responsiveMask & sizeMask & confMask;

gxEarly = [OUT.FitSpatialEarly.centerX].';
gyEarly = [OUT.FitSpatialEarly.centerY].';
gxLate = [OUT.FitSpatialLate.centerX].';
gyLate = [OUT.FitSpatialLate.centerY].';
pSigmaEarly = [OUT.FitSpatialEarly.sigmaPx].';
pSigmaLate = [OUT.FitSpatialLate.sigmaPx].';
pEarly = [OUT.FitSpatialEarly.pValueApprox].';
pLate = [OUT.FitSpatialLate.pValueApprox].';

gaussSizeMaskEarly = isfinite(pSigmaEarly) & (pSigmaEarly <= Opts.MaxGaussianSigmaPx);
gaussSizeMaskLate = isfinite(pSigmaLate) & (pSigmaLate <= Opts.MaxGaussianSigmaPx);

validEarly = baseMask & isfinite(gxEarly) & isfinite(gyEarly) & gaussSizeMaskEarly;
validLate = baseMask & isfinite(gxLate) & isfinite(gyLate) & gaussSizeMaskLate;
sigEarly = validEarly & isfinite(pEarly) & (pEarly < alpha);
sigLate = validLate & isfinite(pLate) & (pLate < alpha);

CMP = struct();
CMP.monkeySuffix = monkeySuffix;
CMP.alpha = alpha;
CMP.nShuffle = nShuffle;
CMP.nEccBins = nEccBins;
CMP.useResponsiveSubset = UseResponsiveSubset;
CMP.subsetLabel = subsetLabel;
CMP.rfFile = rfFile;
CMP.gaussPath = gaussPath;
CMP.colorTunePath = colorTunePath;
CMP.filterOptions = Opts;
CMP.globalSiteInR = rfGlobal;
CMP.rfSizeEq = rfSizeEq;
CMP.rfTestCorr = rfTestCorr;
CMP.rfConfidenceLabel = corrLabel;
CMP.baseMask = baseMask;
CMP.early = summarize_center_window(rfX, rfY, rfEcc, gxEarly, gyEarly, validEarly, sigEarly, nShuffle, nEccBins);
CMP.late = summarize_center_window(rfX, rfY, rfEcc, gxLate, gyLate, validLate, sigLate, nShuffle, nEccBins);

fprintf('IT RF center comparison: RF.m vs Gaussian (%s)\n', char(monkeySuffix));
fprintf('Subset: %s\n', char(subsetLabel));
fprintf('Filters: %s\n', char(describe_filter_opts_local(Opts, corrLabel)));
if isfinite(Opts.MinRFtestCorr) && ~hasFiniteRFcorr
    fprintf('Requested RF confidence filter ignored: %s has no finite values for these IT rows.\n', char(corrLabel));
end
fprintf('Base-mask retained %d / %d sites before early/late Gaussian filters\n', ...
    nnz(baseMask), numel(baseMask));
print_center_summary('Early', CMP.early);
print_center_summary('Late', CMP.late);

figTitle = sprintf('IT RF centers: RF.m vs Gaussian RF (%s)', char(monkeySuffix));
figure(figNum); clf;
set(gcf, 'Color', 'w', 'Name', figTitle, 'NumberTitle', 'off', ...
    'Tag', 'IT_RFcenter_compare_RFm_vs_Gaussian');
fprintf('Opened figure %d: %s\n', figNum, figTitle);

useTiled = exist('tiledlayout', 'file') == 2;
if useTiled
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
end

[xLim, yLim] = choose_center_limits(CMP.early, CMP.late);

if useTiled, nexttile; else, subplot(2, 2, 1); end
draw_center_panel(CMP.early, xLim, yLim, alpha, 'Early');

if useTiled, nexttile; else, subplot(2, 2, 2); end
draw_center_panel(CMP.late, xLim, yLim, alpha, 'Late');

if useTiled, nexttile; else, subplot(2, 2, 3); end
draw_shuffle_panel(CMP.early, alpha, 'Early');

if useTiled, nexttile; else, subplot(2, 2, 4); end
draw_shuffle_panel(CMP.late, alpha, 'Late');

if exist('sgtitle', 'file') == 2
    sgtitle(sprintf('%s | %s | %s', figTitle, char(subsetLabel), ...
        char(describe_filter_opts_local(Opts, corrLabel))));
else
    annotation('textbox', [0.08 0.955 0.84 0.04], ...
        'String', sprintf('%s | %s | %s', figTitle, char(subsetLabel), ...
        char(describe_filter_opts_local(Opts, corrLabel))), ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end
end

function W = summarize_center_window(rfX, rfY, rfEcc, gaussX, gaussY, validMask, sigMask, nShuffle, nEccBins)
W = struct();
W.rfX = rfX(:);
W.rfY = rfY(:);
W.gaussX = gaussX(:);
W.gaussY = gaussY(:);
W.rfEcc = rfEcc(:);
W.valid = validMask(:);
W.sig = sigMask(:);

W.rfXValid = rfX(validMask);
W.rfYValid = rfY(validMask);
W.gaussXValid = gaussX(validMask);
W.gaussYValid = gaussY(validMask);
W.dxValid = W.gaussXValid - W.rfXValid;
W.dyValid = W.gaussYValid - W.rfYValid;
W.distValid = hypot(W.dxValid, W.dyValid);
W.nValid = nnz(validMask);
W.nSig = nnz(sigMask);

if W.nValid >= 2
    W.rX = corr(W.rfXValid, W.gaussXValid, 'Rows', 'complete');
    W.rY = corr(W.rfYValid, W.gaussYValid, 'Rows', 'complete');
    W.rEcc = corr(hypot(W.rfXValid, W.rfYValid), hypot(W.gaussXValid, W.gaussYValid), ...
        'Rows', 'complete');
    W.meanDist = mean(W.distValid, 'omitnan');
    W.medianDist = median(W.distValid, 'omitnan');
    W.meanDx = mean(W.dxValid, 'omitnan');
    W.meanDy = mean(W.dyValid, 'omitnan');
    [W.nullMedianDist, W.binIdx] = build_distance_shuffle_null( ...
        W.rfXValid, W.rfYValid, W.rfEcc(validMask), W.gaussXValid, W.gaussYValid, nShuffle, nEccBins);
    W.shuffleP = (1 + nnz(W.nullMedianDist <= W.medianDist)) / (numel(W.nullMedianDist) + 1);
else
    W.rX = NaN;
    W.rY = NaN;
    W.rEcc = NaN;
    W.meanDist = NaN;
    W.medianDist = NaN;
    W.meanDx = NaN;
    W.meanDy = NaN;
    W.nullMedianDist = NaN;
    W.binIdx = ones(W.nValid, 1);
    W.shuffleP = NaN;
end
end

function [nullMedianDist, binIdx] = build_distance_shuffle_null(rfX, rfY, rfEcc, gaussX, gaussY, nShuffle, nEccBins)
n = numel(rfX);
nullMedianDist = NaN(nShuffle, 1);
binIdx = assign_ecc_bins(rfEcc, nEccBins);

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
    d = hypot(gaussX(perm) - rfX, gaussY(perm) - rfY);
    nullMedianDist(iShuf) = median(d, 'omitnan');
end
end

function binIdx = assign_ecc_bins(ecc, nEccBins)
ecc = ecc(:);
binIdx = ones(size(ecc));
if numel(ecc) < max(4, nEccBins)
    return;
end

edges = prctile(ecc, linspace(0, 100, nEccBins + 1));
edges(1) = -inf;
edges(end) = inf;
edges = make_strictly_increasing(edges);

tmp = discretize(ecc, edges);
if any(isnan(tmp))
    tmp(isnan(tmp)) = 1;
end
binIdx = tmp(:);
end

function edges = make_strictly_increasing(edges)
for i = 2:numel(edges)
    if ~(edges(i) > edges(i-1))
        edges(i) = edges(i-1) + max(1e-6, abs(edges(i-1)) * 1e-6);
    end
end
end

function print_center_summary(label, W)
fprintf(['%s | valid=%d | spatially significant=%d | median dist=%.2f px | ' ...
    'shuffle p=%.4f | rX=%.3f | rY=%.3f | rEcc=%.3f | mean dx=%.2f | mean dy=%.2f\n'], ...
    label, W.nValid, W.nSig, W.medianDist, W.shuffleP, W.rX, W.rY, W.rEcc, W.meanDx, W.meanDy);
end

function [xLim, yLim] = choose_center_limits(W1, W2)
rfX = [W1.rfXValid; W2.rfXValid];
rfY = [W1.rfYValid; W2.rfYValid];
rfX = rfX(isfinite(rfX));
rfY = rfY(isfinite(rfY));

if isempty(rfX) || isempty(rfY)
    xLim = [-150 150];
    yLim = [-150 150];
    return;
end

pad = 15;
xHalf = prctile(abs(rfX), 98);
yHalf = prctile(abs(rfY), 98);

if ~isfinite(xHalf) || xHalf <= 0
    xHalf = max(abs(rfX));
end
if ~isfinite(yHalf) || yHalf <= 0
    yHalf = max(abs(rfY));
end

xHalf = max(100, ceil((xHalf + pad) / 25) * 25);
yHalf = max(100, ceil((yHalf + pad) / 25) * 25);

xLim = [-xHalf xHalf];
yLim = [-yHalf yHalf];
end

function draw_center_panel(W, xLim, yLim, alpha, panelTitle)
hold on;

rfX = W.rfXValid;
rfY = W.rfYValid;
gx = W.gaussXValid;
gy = W.gaussYValid;
sig = W.sig(W.valid);
[gxPlot, gyPlot, gClip] = clip_points_to_box(gx, gy, xLim, yLim);

for i = 1:numel(rfX)
    if sig(i)
        lineColor = [0.85 0.10 0.10];
    else
        lineColor = [0.75 0.75 0.75];
    end
    plot([rfX(i) gxPlot(i)], [rfY(i) gyPlot(i)], '-', 'Color', lineColor, 'LineWidth', 0.8);
end

hRF = plot(rfX, rfY, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', ...
    'MarkerSize', 5, 'LineStyle', 'none');
hGAll = plot(gxPlot(~sig & ~gClip), gyPlot(~sig & ~gClip), 'o', 'MarkerFaceColor', [0.72 0.72 0.72], ...
    'MarkerEdgeColor', 'none', 'MarkerSize', 5, 'LineStyle', 'none');
hGSig = plot(gxPlot(sig & ~gClip), gyPlot(sig & ~gClip), 'o', 'MarkerFaceColor', [0.85 0.10 0.10], ...
    'MarkerEdgeColor', 'k', 'MarkerSize', 6, 'LineStyle', 'none');
plot(gxPlot(~sig & gClip), gyPlot(~sig & gClip), 's', 'MarkerFaceColor', [0.72 0.72 0.72], ...
    'MarkerEdgeColor', [0.35 0.35 0.35], 'MarkerSize', 6, 'LineStyle', 'none');
plot(gxPlot(sig & gClip), gyPlot(sig & gClip), 's', 'MarkerFaceColor', [0.85 0.10 0.10], ...
    'MarkerEdgeColor', 'k', 'MarkerSize', 7, 'LineStyle', 'none');

plot([0 0], yLim, 'k:');
plot(xLim, [0 0], 'k:');

xlim(xLim);
ylim(yLim);
axis equal;
grid on;
xlabel('x center (px)');
ylabel('y center (px)');
title(sprintf('%s | valid=%d | p<%.2f', panelTitle, W.nValid, alpha));
legend([hRF, hGAll, hGSig], {'RF.m center', 'Gaussian center', 'Gaussian center (spatial sig)'}, ...
    'Location', 'southeast');

txt = sprintf('median dist = %.1f px | shuffle p = %.3g', W.medianDist, W.shuffleP);
text(xLim(1) + 0.03 * diff(xLim), yLim(2) - 0.06 * diff(yLim), txt, 'FontSize', 9);
if any(gClip)
    text(xLim(1) + 0.03 * diff(xLim), yLim(1) + 0.05 * diff(yLim), ...
        'square Gaussian markers = clipped to display', 'FontSize', 8);
end
end

function draw_shuffle_panel(W, alpha, panelTitle)
vals = W.nullMedianDist(isfinite(W.nullMedianDist));
if isempty(vals)
    text(0.5, 0.5, 'Not enough valid sites for shuffle test', ...
        'HorizontalAlignment', 'center');
    axis off;
    return;
end

histogram(vals, 35, 'FaceColor', [0.82 0.82 0.82], 'EdgeColor', 'none');
hold on;
yl = ylim;
plot([W.medianDist W.medianDist], yl, 'r-', 'LineWidth', 1.6);
ylim(yl);
grid on;
xlabel('Median paired center distance (px)');
ylabel('N shuffles');
title(sprintf('%s shuffle null | valid=%d | p<%.2f', panelTitle, W.nValid, alpha));
txt = sprintf('obs = %.1f px | shuffle p = %.3g', W.medianDist, W.shuffleP);
text(min(xlim) + 0.03 * diff(xlim), yl(2) - 0.08 * diff(yl), txt, 'FontSize', 9);
end

function rfPath = resolve_rf_file_local(cfg, fileName, monkeySubdir)
candidates = { ...
    fullfile(cfg.matDir, fileName), ...
    fullfile(cfg.dataRoot, monkeySubdir, fileName) ...
};
rfPath = '';
for i = 1:numel(candidates)
    if exist(candidates{i}, 'file') == 2
        rfPath = candidates{i};
        return;
    end
end
error('Compare_RFcenter_IT_Gaussian_vs_RFm:MissingRFFile', ...
    ['Could not find %s. Checked:\n - %s\n - %s'], ...
    fileName, candidates{1}, candidates{2});
end

function [xPlot, yPlot, isClipped] = clip_points_to_box(x, y, xLim, yLim)
xPlot = x;
yPlot = y;
xPlot(xPlot < xLim(1)) = xLim(1);
xPlot(xPlot > xLim(2)) = xLim(2);
yPlot(yPlot < yLim(1)) = yLim(1);
yPlot(yPlot > yLim(2)) = yLim(2);
isClipped = isfinite(x) & isfinite(y) & ...
    ((x < xLim(1)) | (x > xLim(2)) | (y < yLim(1)) | (y > yLim(2)));
end

function Opts = normalize_compare_opts_local(Opts)
defaults = struct();
defaults.MaxRFsizeEq = inf;
defaults.MinRFtestCorr = -inf;
defaults.RFcorrField = 'all_test_corr';
defaults.MaxGaussianSigmaPx = inf;

fn = fieldnames(defaults);
for i = 1:numel(fn)
    if ~isfield(Opts, fn{i}) || isempty(Opts.(fn{i}))
        Opts.(fn{i}) = defaults.(fn{i});
    end
end
end

function [rfTestCorr, corrLabel] = extract_rf_confidence_local(Srf, rfGlobal, corrField)
rfTestCorr = [];
corrLabel = '';
if isfield(Srf, corrField)
    rfTestCorr = double(Srf.(corrField)(:));
    rfTestCorr = rfTestCorr(rfGlobal);
    corrLabel = corrField;
elseif isfield(Srf, 'all_test_corr')
    rfTestCorr = double(Srf.all_test_corr(:));
    rfTestCorr = rfTestCorr(rfGlobal);
    corrLabel = 'all_test_corr';
end
end

function txt = describe_filter_opts_local(Opts, corrLabel)
parts = {};
if isfinite(Opts.MaxRFsizeEq)
    parts{end+1} = sprintf('RF.m eq size <= %.1f', Opts.MaxRFsizeEq); %#ok<AGROW>
end
if isfinite(Opts.MinRFtestCorr)
    if isempty(corrLabel)
        parts{end+1} = sprintf('RF corr >= %.2f requested (field missing)', Opts.MinRFtestCorr); %#ok<AGROW>
    else
        parts{end+1} = sprintf('%s >= %.2f', corrLabel, Opts.MinRFtestCorr); %#ok<AGROW>
    end
end
if isfinite(Opts.MaxGaussianSigmaPx)
    parts{end+1} = sprintf('Gaussian sigma <= %.1f', Opts.MaxGaussianSigmaPx); %#ok<AGROW>
end
if isempty(parts)
    txt = 'no extra size/confidence filter';
else
    txt = strjoin(parts, ', ');
end
end
