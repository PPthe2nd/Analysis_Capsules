function CMP = Compare_RFsize_IT_Gaussian_vs_RFm(Monkey)
% Compare RF.m size estimates to the best Gaussian sigma in IT.

if nargin < 1 || isempty(Monkey)
    Monkey = 1; % 1 = Nilson, 2 = Figaro
end

alpha = 0.05;
figNum = 103;

cfg = config();

if Monkey == 1
    monkeySuffix = "N";
    rfFile = resolve_rf_file_local(cfg, 'THINGS_RF1s_N.mat', 'Mr Nilson');
elseif Monkey == 2
    monkeySuffix = "F";
    rfFile = resolve_rf_file_local(cfg, 'THINGS_RF1s_F.mat', 'Figaro');
else
    error('Compare_RFsize_IT_Gaussian_vs_RFm:InvalidMonkey', ...
        'Monkey must be 1 (Nilson) or 2 (Figaro).');
end

gaussPath = fullfile(cfg.matDir, sprintf('GaussianOccupancy_Tuning_IT_%s.mat', char(monkeySuffix)));
assert(exist(rfFile, 'file') == 2, 'Missing RF file: %s', rfFile);
assert(exist(gaussPath, 'file') == 2, ...
    'Missing %s. Run GaussianOccupancy_Tuning_IT.m first.', gaussPath);

Srf = load(rfFile, 'all_centrex', 'all_centrey', 'all_szx', 'all_szy');
Sg = load(gaussPath, 'OUT');
assert(isfield(Sg, 'OUT') && isstruct(Sg.OUT), '%s must contain struct OUT.', gaussPath);
OUT = Sg.OUT;
assert(isfield(OUT, 'RFrange') && isfield(OUT, 'FitSpatialEarly') && isfield(OUT, 'FitSpatialLate'), ...
    'Gaussian OUT is missing required fields.');

x = double(Srf.all_centrex(:));
y = double(Srf.all_centrey(:));
sx = double(Srf.all_szx(:));
sy = double(Srf.all_szy(:));
rfGlobal = OUT.RFrange(:);

assert(max(rfGlobal) <= numel(x), ...
    'RF indices exceed the RF file length.');

rfValid = isfinite(x(rfGlobal)) & isfinite(y(rfGlobal)) & isfinite(sx(rfGlobal)) & isfinite(sy(rfGlobal)) & ...
    (sx(rfGlobal) > 0) & (sy(rfGlobal) > 0) & abs(x(rfGlobal)) <= 200 & abs(y(rfGlobal)) <= 200;

% Use the geometric mean of the x/y RF sigmas as an isotropic RF.m size proxy.
rfSizeEq = sqrt(sx(rfGlobal) .* sy(rfGlobal));

sigmaEarly = [OUT.FitSpatialEarly.sigmaPx].';
sigmaLate = [OUT.FitSpatialLate.sigmaPx].';
pEarly = [OUT.FitSpatialEarly.pValueApprox].';
pLate = [OUT.FitSpatialLate.pValueApprox].';

validEarly = rfValid & isfinite(rfSizeEq) & isfinite(sigmaEarly) & (rfSizeEq > 0) & (sigmaEarly > 0);
validLate = rfValid & isfinite(rfSizeEq) & isfinite(sigmaLate) & (rfSizeEq > 0) & (sigmaLate > 0);
sigEarly = validEarly & isfinite(pEarly) & (pEarly < alpha);
sigLate = validLate & isfinite(pLate) & (pLate < alpha);

CMP = struct();
CMP.monkeySuffix = monkeySuffix;
CMP.alpha = alpha;
CMP.rfFile = rfFile;
CMP.gaussPath = gaussPath;
CMP.globalSiteInR = rfGlobal;
CMP.rfValid = rfValid;
CMP.rfSizeEq = rfSizeEq;
CMP.early = summarize_window(rfSizeEq, sigmaEarly, validEarly, sigEarly);
CMP.late = summarize_window(rfSizeEq, sigmaLate, validLate, sigLate);

fprintf('IT RF size comparison: RF.m vs Gaussian (%s)\n', char(monkeySuffix));
print_size_summary('Early', CMP.early);
print_size_summary('Late', CMP.late);

figTitle = sprintf('IT RF size: RF.m vs Gaussian RF (%s)', char(monkeySuffix));
figure(figNum); clf;
set(gcf, 'Color', 'w', 'Name', figTitle, 'NumberTitle', 'off', ...
    'Tag', 'IT_RFsize_compare_RFm_vs_Gaussian');
fprintf('Opened figure %d: %s\n', figNum, figTitle);

useTiled = exist('tiledlayout', 'file') == 2;
if useTiled
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
end

allX = [CMP.early.xValid; CMP.late.xValid];
allY = [CMP.early.yValid; CMP.late.yValid];
xLim = choose_positive_limits(allX, 0.98);
yLim = choose_positive_limits(allY, 1.00);

if useTiled, nexttile; else, subplot(1, 2, 1); end
draw_size_panel(CMP.early, xLim, yLim, alpha, 'Early');

if useTiled, nexttile; else, subplot(1, 2, 2); end
draw_size_panel(CMP.late, xLim, yLim, alpha, 'Late');

if exist('sgtitle', 'file') == 2
    sgtitle(figTitle);
else
    annotation('textbox', [0.12 0.955 0.76 0.04], 'String', figTitle, ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end
end

function W = summarize_window(rfSizeEq, gaussSigma, validMask, sigMask)
W = struct();
W.xAll = rfSizeEq(:);
W.yAll = gaussSigma(:);
W.valid = validMask(:);
W.sig = sigMask(:);
W.xValid = rfSizeEq(validMask);
W.yValid = gaussSigma(validMask);
W.xSig = rfSizeEq(sigMask);
W.ySig = gaussSigma(sigMask);
W.nValid = nnz(validMask);
W.nSig = nnz(sigMask);
if W.nValid >= 2
    W.rhoSpearman = corr(W.xValid, W.yValid, 'Type', 'Spearman', 'Rows', 'complete');
    W.rPearsonLog = corr(log10(W.xValid), log10(W.yValid), 'Rows', 'complete');
    W.medianRatio = median(W.yValid ./ W.xValid, 'omitnan');
else
    W.rhoSpearman = NaN;
    W.rPearsonLog = NaN;
    W.medianRatio = NaN;
end
end

function print_size_summary(label, W)
fprintf('%s | valid=%d | spatially significant=%d | Spearman rho=%.3f | log-Pearson r=%.3f | median ratio=%.3f\n', ...
    label, W.nValid, W.nSig, W.rhoSpearman, W.rPearsonLog, W.medianRatio);
end

function draw_size_panel(W, xLim, yLim, alpha, panelTitle)
[xPlot, xClip] = clip_to_positive_range(W.xValid, xLim);
[yPlot, yClip] = clip_to_positive_range(W.yValid, yLim);
sig = W.sig(W.valid);
isClipped = xClip | yClip;

hold on;
hAll = plot(xPlot(~sig & ~isClipped), yPlot(~sig & ~isClipped), 'o', ...
    'MarkerFaceColor', [0.78 0.78 0.78], 'MarkerEdgeColor', 'none', ...
    'MarkerSize', 5, 'LineStyle', 'none');
hSig = plot(xPlot(sig & ~isClipped), yPlot(sig & ~isClipped), 'o', ...
    'MarkerFaceColor', [0.85 0.10 0.10], 'MarkerEdgeColor', 'k', ...
    'MarkerSize', 6, 'LineStyle', 'none');
hClip = plot(xPlot(isClipped), yPlot(isClipped), 's', ...
    'MarkerFaceColor', 'none', 'MarkerEdgeColor', [0.35 0.35 0.35], ...
    'MarkerSize', 6, 'LineStyle', 'none');
hDiag = plot([max(xLim(1), yLim(1)) min(xLim(2), yLim(2))], [max(xLim(1), yLim(1)) min(xLim(2), yLim(2))], ...
    'k--', 'LineWidth', 1.2);

set(gca, 'XScale', 'log', 'YScale', 'log');
xlim(xLim);
ylim(yLim);
axis square;
grid on;
xlabel('RF.m equivalent size  sqrt(sz_x * sz_y)');
ylabel('Best Gaussian sigma (px)');
title(sprintf('%s | valid=%d | p<%.2f', panelTitle, W.nValid, alpha));
legend([hAll, hSig, hClip, hDiag], ...
    {'All valid sites', 'Spatially significant Gaussian', 'Clipped outlier', 'y = x'}, ...
    'Location', 'southeast');

txt = sprintf('\\rho_s = %.2f | r_{log} = %.2f | median y/x = %.2f', ...
    W.rhoSpearman, W.rPearsonLog, W.medianRatio);
text(exp(log(xLim(1)) + 0.04 * (log(xLim(2)) - log(xLim(1)))), ...
    exp(log(yLim(2)) - 0.08 * (log(yLim(2)) - log(yLim(1)))), ...
    txt, 'FontSize', 9);
end

function lim = choose_positive_limits(vals, pctKeep)
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

function [valsPlot, isClipped] = clip_to_positive_range(vals, lim)
valsPlot = vals;
valsPlot(vals < lim(1)) = lim(1);
valsPlot(vals > lim(2)) = lim(2);
isClipped = isfinite(vals) & ((vals < lim(1)) | (vals > lim(2)));
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
error('Compare_RFsize_IT_Gaussian_vs_RFm:MissingRFFile', ...
    ['Could not find %s. Checked:\n - %s\n - %s'], ...
    fileName, candidates{1}, candidates{2});
end
