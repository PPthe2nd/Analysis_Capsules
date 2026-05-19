function OUT = compare_area_session_responses_3bin(areaName, Monkey, Puser)
% COMPARE_AREA_SESSION_RESPONSES_3BIN
% Session-1 vs session-2 comparison for V1/V4 line-task responses in the
% standard three coarse windows [-200 0], [40 240], [300 500] ms.

if nargin < 2 || isempty(Monkey)
    Monkey = 1;
end
if nargin < 3 || isempty(Puser)
    Puser = struct();
end

areaName = upper(string(areaName));
assert(any(areaName == ["V1","V4"]), 'areaName must be ''V1'' or ''V4''.');

P = struct();
P.onlyCorrect = true;
P.correctCol = 9;
P.correctVal = 1;
P.dayCol = 11;
P.sessions = [1 2];
P.expectedMaxStim = 384;
P.chunkTrials = 200;
P.SNRthr = 0.7;
P.MinObjectStim = 1;
P.useCache = true;
P.saveResult = true;
P.plotFigure = true;

if ~isempty(Puser)
    fn = fieldnames(Puser);
    for i = 1:numel(fn)
        P.(fn{i}) = Puser.(fn{i});
    end
end

cfg = config();

if Monkey == 1
    monkeySuffix = "N";
elseif Monkey == 2
    monkeySuffix = "F";
else
    error('Monkey must be 1 (Nilson) or 2 (Figaro).');
end

switch char(areaName)
    case 'V1'
        tallFile = sprintf('Tall_V1_lines_%s.mat', char(monkeySuffix));
        tallField = 'Tall_V1';
        subsetLabel = sprintf('bestSNR > %.2f', P.SNRthr);
    case 'V4'
        tallFile = sprintf('Tall_V4_lines_%s.mat', char(monkeySuffix));
        tallField = 'Tall_V4';
        subsetLabel = sprintf('hasObjectRF & bestSNR > %.2f', P.SNRthr);
end

resp3binFile = sprintf('SNR_capsules_%s_d12.mat', char(monkeySuffix));
tallPath = fullfile(cfg.matDir, tallFile);
resp3binPath = fullfile(cfg.matDir, resp3binFile);
outPath = fullfile(cfg.resultsDir, sprintf('Compare_%s_SessionResponses_3bin_%s.mat', char(areaName), char(monkeySuffix)));

assert(exist(tallPath, 'file') == 2, 'Missing %s.', tallPath);
assert(exist(resp3binPath, 'file') == 2, 'Missing %s.', resp3binPath);

varsTall = who('-file', tallPath);
if any(strcmp(varsTall, 'RFrange'))
    Sgeo = load(tallPath, tallField, 'RFrange');
else
    Sgeo = load(tallPath, tallField);
end
assert(isfield(Sgeo, tallField) && isstruct(Sgeo.(tallField)), '%s must contain %s.', tallPath, tallField);
Tall = Sgeo.(tallField);
if isfield(Sgeo, 'RFrange') && ~isempty(Sgeo.RFrange)
    RFrange = Sgeo.RFrange(:);
elseif areaName == "V1"
    % The canonical V1 line-task scripts use the first 512 channels
    % directly rather than storing an explicit RFrange in Tall_V1_lines_*.mat.
    RFrange = (1:512).';
else
    error('%s must contain RFrange.', tallPath);
end
nArea = numel(RFrange);
siteRows = (1:nArea).';

stimNums = arrayfun(@(x) x.stimNum, Tall(:));
[stimNumsSorted, ordStim] = sort(stimNums(:));
assert(all(stimNumsSorted(:).' == 1:numel(Tall)), ...
    '%s.stimNum must cover 1..%d exactly.', tallField, numel(Tall));
Tall = Tall(ordStim);
nStim = numel(Tall);

R3full = load_capsules_struct_exclusion_aware(resp3binPath, monkeySuffix, 'cfg', cfg);
assert(size(R3full.meanAct,1) >= max(RFrange), ...
    'Response file %s has only %d rows; need at least %d for %s.', ...
    resp3binPath, size(R3full.meanAct,1), max(RFrange), char(areaName));
R3 = localize_response_rows_local(R3full, RFrange);
assert(size(R3.meanAct,1) == nArea && size(R3.meanAct,2) == nStim, ...
    'Localized %s response rows do not match geometry.', char(areaName));
timeWindows = double(R3.timeWindows);
assert(size(timeWindows,1) == 3, 'Expected exactly 3 coarse time windows.');

SNR = compute_snr_per_color_sites(R3, Tall, siteRows, 'Verbose', false);
SNRmat = [SNR.yellowEarly(siteRows), SNR.yellowLate(siteRows), ...
          SNR.purpleEarly(siteRows), SNR.purpleLate(siteRows)];
bestSNR = max(SNRmat, [], 2, 'omitnan');

nObjectStim = zeros(nArea,1);
for stim = 1:nStim
    T = Tall(stim).T;
    assign = string(T.assignment(siteRows));
    nObjectStim = nObjectStim + (assign == "target") + (assign == "distractor");
end
hasObjectRF = nObjectStim >= P.MinObjectStim;

switch char(areaName)
    case 'V1'
        keepResponsive = isfinite(bestSNR) & (bestSNR > P.SNRthr);
    case 'V4'
        keepResponsive = hasObjectRF & isfinite(bestSNR) & (bestSNR > P.SNRthr);
end

[Split, ~, cachePath] = load_capsules_day_split(resp3binPath, monkeySuffix, ...
    'cfg', cfg, ...
    'sessions', P.sessions, ...
    'onlyCorrect', P.onlyCorrect, ...
    'correctCol', P.correctCol, ...
    'correctVal', P.correctVal, ...
    'dayCol', P.dayCol, ...
    'chunkTrials', P.chunkTrials, ...
    'expectedMaxStim', P.expectedMaxStim, ...
    'useCache', P.useCache);

muBySess = nan(nArea, size(timeWindows,1), numel(P.sessions));
for iSess = 1:numel(P.sessions)
    muBySess(:,:,iSess) = weighted_site_window_means(Split(iSess).meanAct(RFrange,:,:), Split(iSess).nTrials);
end

preBySess = squeeze(muBySess(:,1,:));
earlyBySess = squeeze(muBySess(:,2,:));
lateBySess = squeeze(muBySess(:,3,:));
evokedEarly = earlyBySess - preBySess;
evokedLate = lateBySess - preBySess;

lateDrop = evokedLate(:,2) - evokedLate(:,1);
preShift = preBySess(:,2) - preBySess(:,1);
[~, exLateLocal] = min(lateDrop);
[~, exPreLocal] = min(preShift);
if exPreLocal == exLateLocal
    tmp = preShift;
    tmp(exLateLocal) = inf;
    [~, exPreLocal] = min(tmp);
end

OUT = struct();
OUT.P = P;
OUT.areaName = areaName;
OUT.monkeySuffix = monkeySuffix;
OUT.siteGlobal = RFrange(:);
OUT.keepResponsive = keepResponsive(:);
OUT.hasObjectRF = hasObjectRF(:);
OUT.bestSNR = bestSNR(:);
OUT.nObjectStim = nObjectStim(:);
OUT.timeWindows = timeWindows;
OUT.sessionIds = P.sessions(:);
OUT.muBySess = muBySess;
OUT.preBySess = preBySess;
OUT.earlyBySess = earlyBySess;
OUT.lateBySess = lateBySess;
OUT.evokedEarly = evokedEarly;
OUT.evokedLate = evokedLate;
OUT.exampleLateDropSiteGlobal = RFrange(exLateLocal);
OUT.exampleLateDropSiteLocal = exLateLocal;
OUT.examplePreShiftSiteGlobal = RFrange(exPreLocal);
OUT.examplePreShiftSiteLocal = exPreLocal;
OUT.daySplitCachePath = cachePath;
OUT.subsetLabel = subsetLabel;
OUT.summary = build_summary(preBySess, earlyBySess, lateBySess, evokedEarly, evokedLate, keepResponsive);

fprintf('%s session 3-bin comparison for monkey %s\n', char(areaName), char(monkeySuffix));
fprintf('All %s sites: %d | Responsive subset: %d (%s)\n', ...
    char(areaName), nArea, nnz(keepResponsive), subsetLabel);
fprintf('Late evoked median (all): sess%d=%.4f, sess%d=%.4f\n', ...
    P.sessions(1), OUT.summary.all.lateEvokedMedian(1), P.sessions(2), OUT.summary.all.lateEvokedMedian(2));
fprintf('Late evoked median (responsive): sess%d=%.4f, sess%d=%.4f\n', ...
    P.sessions(1), OUT.summary.responsive.lateEvokedMedian(1), P.sessions(2), OUT.summary.responsive.lateEvokedMedian(2));
fprintf('Responsive sites with sess%d late-evoked < sess%d late-evoked: %.1f%%\n', ...
    P.sessions(2), P.sessions(1), 100 * OUT.summary.responsive.fracLateDrop);
fprintf('Responsive sites with sess%d late-evoked <= 0 and sess%d late-evoked > 0: %.1f%%\n', ...
    P.sessions(2), P.sessions(1), 100 * OUT.summary.responsive.fracLateGone);
fprintf('Largest late drop site: global %d (local %d)\n', OUT.exampleLateDropSiteGlobal, OUT.exampleLateDropSiteLocal);
fprintf('Largest negative pre shift site: global %d (local %d)\n', OUT.examplePreShiftSiteGlobal, OUT.examplePreShiftSiteLocal);

if P.saveResult
    save(outPath, 'OUT', '-v7.3');
    fprintf('Saved %s session comparison to %s\n', char(areaName), outPath);
end

if P.plotFigure
    make_summary_figure(OUT);
end
end

function R_loc = localize_response_rows_local(R_full, siteGlobal)
R_loc = R_full;
R_loc.meanAct = R_full.meanAct(siteGlobal, :, :);
R_loc.meanSqAct = R_full.meanSqAct(siteGlobal, :, :);
if ismatrix(R_full.nTrials) && size(R_full.nTrials,1) >= max(siteGlobal)
    R_loc.nTrials = R_full.nTrials(siteGlobal, :);
else
    R_loc.nTrials = R_full.nTrials;
end
end

function mu = weighted_site_window_means(meanAct, nTrials)
nArea = size(meanAct, 1);
nWin = size(meanAct, 3);
w = double(nTrials(:));
mu = nan(nArea, nWin);
for s = 1:nArea
    for wIdx = 1:nWin
        x = squeeze(meanAct(s,:,wIdx)).';
        idx = isfinite(x) & isfinite(w) & (w > 0);
        if any(idx)
            mu(s, wIdx) = sum(w(idx) .* x(idx)) / sum(w(idx));
        end
    end
end
end

function summary = build_summary(preBySess, earlyBySess, lateBySess, evokedEarly, evokedLate, keepResponsive)
summary = struct();
summary.all = one_summary(preBySess, earlyBySess, lateBySess, evokedEarly, evokedLate, true(size(keepResponsive)));
summary.responsive = one_summary(preBySess, earlyBySess, lateBySess, evokedEarly, evokedLate, keepResponsive);
end

function S = one_summary(preBySess, earlyBySess, lateBySess, evokedEarly, evokedLate, mask)
mask = mask(:) & all(isfinite(evokedLate), 2);
S = struct();
S.nSites = nnz(mask);
S.preMedian = median(preBySess(mask,:), 1, 'omitnan');
S.earlyMedian = median(earlyBySess(mask,:), 1, 'omitnan');
S.lateMedian = median(lateBySess(mask,:), 1, 'omitnan');
S.earlyEvokedMedian = median(evokedEarly(mask,:), 1, 'omitnan');
S.lateEvokedMedian = median(evokedLate(mask,:), 1, 'omitnan');
S.fracLateDrop = mean(evokedLate(mask,2) < evokedLate(mask,1), 'omitnan');
S.fracLateGone = mean(evokedLate(mask,2) <= 0 & evokedLate(mask,1) > 0, 'omitnan');
end

function make_summary_figure(OUT)
sess1 = OUT.sessionIds(1);
sess2 = OUT.sessionIds(2);
resp = OUT.keepResponsive(:);
lateLocal = OUT.exampleLateDropSiteLocal;
preLocal = OUT.examplePreShiftSiteLocal;

fig = figure('Color', 'w', ...
    'Name', sprintf('%s session 3-bin comparison %s', char(OUT.areaName), char(OUT.monkeySuffix)));
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile;
scatter_with_subset(ax1, OUT.preBySess(:,1), OUT.preBySess(:,2), resp, lateLocal, preLocal);
xlabel(ax1, sprintf('sess%d pre', sess1));
ylabel(ax1, sprintf('sess%d pre', sess2));
title(ax1, 'Pre window');

ax2 = nexttile;
scatter_with_subset(ax2, OUT.evokedEarly(:,1), OUT.evokedEarly(:,2), resp, lateLocal, preLocal);
xlabel(ax2, sprintf('sess%d early-pre', sess1));
ylabel(ax2, sprintf('sess%d early-pre', sess2));
title(ax2, 'Early evoked');

ax3 = nexttile;
scatter_with_subset(ax3, OUT.evokedLate(:,1), OUT.evokedLate(:,2), resp, lateLocal, preLocal);
xlabel(ax3, sprintf('sess%d late-pre', sess1));
ylabel(ax3, sprintf('sess%d late-pre', sess2));
title(ax3, 'Late evoked');

ax4 = nexttile;
hold(ax4, 'on');
deltaAll = OUT.evokedLate(:,2) - OUT.evokedLate(:,1);
histogram(ax4, deltaAll, 24, 'FaceColor', [0.80 0.80 0.80], 'EdgeColor', 'none', 'FaceAlpha', 0.9);
histogram(ax4, deltaAll(resp), 20, 'FaceColor', [0.20 0.45 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.65);
line(ax4, [0 0], ylim(ax4), 'Color', [0 0 0], 'LineStyle', '--', 'LineWidth', 1.0);
xlabel(ax4, sprintf('late-pre: sess%d - sess%d', sess2, sess1));
ylabel(ax4, 'N sites');
title(ax4, 'Late evoked session change');
legend(ax4, {'all sites', 'responsive subset'}, 'Location', 'best', 'Box', 'off');
grid(ax4, 'on');

annotation(fig, 'textbox', [0.05 0.955 0.90 0.04], ...
    'String', sprintf(['%s session split using the standard 3 windows [%d %d], [%d %d], [%d %d] ms. ' ...
    'Blue = %s. Red marker = largest late drop site %d; magenta marker = strongest negative pre-shift site %d.'], ...
    char(OUT.areaName), ...
    OUT.timeWindows(1,1), OUT.timeWindows(1,2), ...
    OUT.timeWindows(2,1), OUT.timeWindows(2,2), ...
    OUT.timeWindows(3,1), OUT.timeWindows(3,2), ...
    OUT.subsetLabel, OUT.exampleLateDropSiteGlobal, OUT.examplePreShiftSiteGlobal), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'FontWeight', 'bold', 'FontSize', 10);
end

function scatter_with_subset(ax, x, y, keepMask, lateLocal, preLocal)
hold(ax, 'on');
scatter(ax, x, y, 22, [0.75 0.75 0.75], 'filled', ...
    'MarkerFaceAlpha', 0.75, 'MarkerEdgeAlpha', 0.75);
scatter(ax, x(keepMask), y(keepMask), 28, [0.20 0.45 0.85], 'filled', ...
    'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0.9);
if isfinite(lateLocal)
    scatter(ax, x(lateLocal), y(lateLocal), 72, [0.85 0.20 0.20], ...
        'p', 'filled', 'MarkerEdgeColor', [1 1 1], 'LineWidth', 0.8);
end
if isfinite(preLocal)
    scatter(ax, x(preLocal), y(preLocal), 72, [0.75 0.10 0.75], ...
        'd', 'filled', 'MarkerEdgeColor', [1 1 1], 'LineWidth', 0.8);
end
lims = finite_lims([x(:); y(:)]);
if all(isfinite(lims))
    pad = 0.05 * max(1e-6, diff(lims));
    lims = [lims(1)-pad, lims(2)+pad];
    xlim(ax, lims);
    ylim(ax, lims);
    line(ax, lims, lims, 'Color', [0.4 0.4 0.4], 'LineStyle', '--', 'LineWidth', 1.0);
end
grid(ax, 'on');
axis(ax, 'square');
end

function lims = finite_lims(x)
x = x(isfinite(x));
if isempty(x)
    lims = [NaN NaN];
else
    lims = [min(x) max(x)];
end
end
