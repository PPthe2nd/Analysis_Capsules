function CMP = Compare_ColorPreference_IT_Gaussian_vs_RFm(Monkey)
% Compare IT color preference from RF.m-center tuning and Gaussian RF tuning.

if nargin < 1 || isempty(Monkey)
    Monkey = 1; % 1 = Nilson, 2 = Figaro
end

alpha = 0.05;
plotRange = [-0.6 0.6];
figNum = 102;

cfg = config();

if Monkey == 1
    monkeySuffix = "N";
elseif Monkey == 2
    monkeySuffix = "F";
else
    error('Compare_ColorPreference_IT_Gaussian_vs_RFm:InvalidMonkey', ...
        'Monkey must be 1 (Nilson) or 2 (Figaro).');
end

rfPath = fullfile(cfg.matDir, sprintf('ColorTune_balanced_IT_%s.mat', char(monkeySuffix)));
gaussPath = fullfile(cfg.matDir, sprintf('GaussianOccupancy_Tuning_IT_%s.mat', char(monkeySuffix)));

assert(exist(rfPath, 'file') == 2, ...
    'Missing %s. Run ColorTuning_Capsules_IT.m first.', rfPath);
assert(exist(gaussPath, 'file') == 2, ...
    'Missing %s. Run GaussianOccupancy_Tuning_IT.m first.', gaussPath);

Srf = load(rfPath, 'ColorTune');
Sg = load(gaussPath, 'OUT');
assert(isfield(Srf, 'ColorTune') && isstruct(Srf.ColorTune), ...
    '%s must contain struct ColorTune.', rfPath);
assert(isfield(Sg, 'OUT') && isstruct(Sg.OUT), ...
    '%s must contain struct OUT.', gaussPath);

ColorTune = Srf.ColorTune;
OUT = Sg.OUT;

assert(isfield(ColorTune, 'RFrange') && ~isempty(ColorTune.RFrange), ...
    'ColorTune must contain RFrange.');
assert(isfield(OUT, 'RFrange') && ~isempty(OUT.RFrange), ...
    'OUT must contain RFrange.');
assert(isfield(OUT, 'PairColorEarly') && isfield(OUT, 'PairColorLate'), ...
    'OUT must contain PairColorEarly and PairColorLate.');
assert(isfield(OUT.PairColorEarly, 'colorIndex') && isfield(OUT.PairColorEarly, 'pValuePooled') && ...
       isfield(OUT.PairColorLate, 'colorIndex') && isfield(OUT.PairColorLate, 'pValuePooled'), ...
    ['Gaussian IT output is missing pooled RF color fields. ' ...
     'Rerun GaussianOccupancy_Tuning_IT.m.']);

rfRFrange = ColorTune.RFrange(:);
gaussRFrange = OUT.RFrange(:);
[globalSiteInR, idxRF, idxGauss] = intersect(rfRFrange, gaussRFrange, 'stable');
assert(~isempty(globalSiteInR), 'No overlapping IT sites found between RF-center and Gaussian analyses.');

CMP = struct();
CMP.monkeySuffix = monkeySuffix;
CMP.alpha = alpha;
CMP.globalSiteInR = globalSiteInR;
CMP.idxRF = idxRF;
CMP.idxGauss = idxGauss;
CMP.early = build_compare_window(ColorTune.early, OUT.PairColorEarly, idxRF, idxGauss, alpha);
CMP.late = build_compare_window(ColorTune.late, OUT.PairColorLate, idxRF, idxGauss, alpha);

fprintf('IT color-preference comparison (%s)\n', char(monkeySuffix));
fprintf('Matched IT sites: %d\n', numel(globalSiteInR));
print_window_summary('Early', CMP.early);
print_window_summary('Late', CMP.late);

figTitle = sprintf('IT color preference: RF.m center vs Gaussian RF (%s)', char(monkeySuffix));
figure(figNum); clf;
set(gcf, 'Color', 'w', 'Name', figTitle, 'NumberTitle', 'off', ...
    'Tag', 'IT_color_pref_compare_RFm_vs_Gaussian');

useTiled = exist('tiledlayout', 'file') == 2;
if useTiled
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
end

if useTiled, nexttile; else, subplot(1, 2, 1); end
draw_compare_panel(CMP.early, plotRange, alpha, 'Early');

if useTiled, nexttile; else, subplot(1, 2, 2); end
draw_compare_panel(CMP.late, plotRange, alpha, 'Late');

if exist('sgtitle', 'file') == 2
    sgtitle(figTitle);
else
    annotation('textbox', [0.12 0.955 0.76 0.04], 'String', figTitle, ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end
end

function W = build_compare_window(ColorWin, PairColor, idxRF, idxGauss, alpha)
rfCI = ColorWin.colorIndex(idxRF);
rfP = ColorWin.p(idxRF);
gaussCI = [PairColor(idxGauss).colorIndex].';
gaussP = [PairColor(idxGauss).pValuePooled].';

valid = isfinite(rfCI) & isfinite(gaussCI);
W = struct();
W.rfCI = rfCI(:);
W.gaussCI = gaussCI(:);
W.rfP = rfP(:);
W.gaussP = gaussP(:);
W.valid = valid(:);
W.rfSig = valid(:) & isfinite(rfP(:)) & (rfP(:) < alpha);
W.gaussSig = valid(:) & isfinite(gaussP(:)) & (gaussP(:) < alpha);
W.noneSig = W.valid & ~W.rfSig & ~W.gaussSig;
W.rfOnly = W.valid & W.rfSig & ~W.gaussSig;
W.gaussOnly = W.valid & ~W.rfSig & W.gaussSig;
W.bothSig = W.valid & W.rfSig & W.gaussSig;
W.nValid = nnz(W.valid);

if nnz(W.valid) >= 2
    W.rPearson = corr(W.rfCI(W.valid), W.gaussCI(W.valid), 'Rows', 'complete');
else
    W.rPearson = NaN;
end

sameSign = sign(W.rfCI) == sign(W.gaussCI);
nonZero = (W.rfCI ~= 0) & (W.gaussCI ~= 0);
W.nSameSignBothSig = nnz(W.bothSig & sameSign & nonZero);
W.nBothSigNonZero = nnz(W.bothSig & nonZero);
end

function print_window_summary(label, W)
fprintf('%s | valid=%d | RF only sig=%d | Gaussian only sig=%d | both sig=%d | r=%.3f\n', ...
    label, W.nValid, nnz(W.rfOnly), nnz(W.gaussOnly), nnz(W.bothSig), W.rPearson);
if W.nBothSigNonZero > 0
    fprintf('%s | same sign among both-significant nonzero sites: %d / %d\n', ...
        label, W.nSameSignBothSig, W.nBothSigNonZero);
end
end

function draw_compare_panel(W, plotRange, alpha, panelTitle)
lo = plotRange(1);
hi = plotRange(2);

x = W.rfCI(W.valid);
y = W.gaussCI(W.valid);
catNone = W.noneSig(W.valid);
catRF = W.rfOnly(W.valid);
catGauss = W.gaussOnly(W.valid);
catBoth = W.bothSig(W.valid);

[xPlot, xClip] = clip_to_range(x, lo, hi);
[yPlot, yClip] = clip_to_range(y, lo, hi);
isClipped = xClip | yClip;

hold on;
hNone = plot_category(xPlot, yPlot, catNone & ~isClipped, 'o', [0.75 0.75 0.75], 'none', 5);
hRF = plot_category(xPlot, yPlot, catRF & ~isClipped, 'o', [0.20 0.45 0.85], 'k', 6);
hGauss = plot_category(xPlot, yPlot, catGauss & ~isClipped, 'o', [0.90 0.55 0.15], 'k', 6);
hBoth = plot_category(xPlot, yPlot, catBoth & ~isClipped, 'o', [0.85 0.10 0.10], 'k', 6);
hClip = plot_category(xPlot, yPlot, isClipped, 's', 'none', [0.35 0.35 0.35], 6);
hDiag = plot([lo hi], [lo hi], 'k--', 'LineWidth', 1.2);
xline(0, 'k:');
yline(0, 'k:');
xlim([lo hi]);
ylim([lo hi]);
axis square;
grid on;
xlabel('RF.m center color index (yellow - purple)');
ylabel('Gaussian RF color index (yellow - purple)');
title(sprintf('%s | valid=%d | p<%.2f', panelTitle, W.nValid, alpha));
legend([hRF, hGauss, hBoth], ...
    {'RF.m center only', 'Gaussian only', 'Both significant'}, ...
    'Location', 'southeast');

txt = sprintf('r = %.2f', W.rPearson);
if W.nBothSigNonZero > 0
    txt = sprintf('%s | same sign (both sig): %d/%d', txt, W.nSameSignBothSig, W.nBothSigNonZero);
end
text(lo + 0.03 * (hi - lo), hi - 0.06 * (hi - lo), txt, 'FontSize', 9);
end

function h = plot_category(x, y, mask, marker, faceColor, edgeColor, markerSize)
if ischar(faceColor) || (isstring(faceColor) && strcmp(faceColor, "none"))
    h = plot(x(mask), y(mask), marker, 'MarkerFaceColor', 'none', ...
        'MarkerEdgeColor', edgeColor, 'MarkerSize', markerSize, 'LineStyle', 'none');
else
    h = plot(x(mask), y(mask), marker, 'MarkerFaceColor', faceColor, ...
        'MarkerEdgeColor', edgeColor, 'MarkerSize', markerSize, 'LineStyle', 'none');
end
if ~any(mask)
    if ischar(faceColor) || (isstring(faceColor) && strcmp(faceColor, "none"))
        h = plot(nan, nan, marker, 'MarkerFaceColor', 'none', ...
            'MarkerEdgeColor', edgeColor, 'MarkerSize', markerSize, 'LineStyle', 'none');
    else
        h = plot(nan, nan, marker, 'MarkerFaceColor', faceColor, ...
            'MarkerEdgeColor', edgeColor, 'MarkerSize', markerSize, 'LineStyle', 'none');
    end
end
end

function [valsPlot, isClipped] = clip_to_range(vals, lo, hi)
valsPlot = vals;
valsPlot(vals < lo) = lo;
valsPlot(vals > hi) = hi;
isClipped = isfinite(vals) & ((vals < lo) | (vals > hi));
end
