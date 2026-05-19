function OUT = Compare_IT_SessionResponses_3bin(Puser)
% COMPARE_IT_SESSIONRESPONSES_3BIN
% Compare day-1 vs day-2 IT responses using the same coarse 3 windows as
% the standard SNR_capsules_*_d12.mat analysis. This is intended as a
% lightweight session-drift check that stays aligned to the existing
% 3-window pipeline.

P = struct();
P.Monkey = 1;                  % 1 = Nilson, 2 = Figaro
P.onlyCorrect = true;
P.correctCol = 9;
P.correctVal = 1;
P.dayCol = 11;
P.sessions = [1 2];
P.expectedMaxStim = 384;
P.chunkTrials = 200;
P.MinObjectStim = 1;
P.TopKQuartets = 5;
P.RespThr = 0.7;
P.useCache = true;
P.saveResult = true;
P.plotFigure = true;
P.exampleSiteGlobal = 868;     % current suspicious example site

if nargin >= 1 && ~isempty(Puser)
    fn = fieldnames(Puser);
    for i = 1:numel(fn)
        P.(fn{i}) = Puser.(fn{i});
    end
end

cfg = config();

if P.Monkey == 1
    monkeySuffix = "N";
    monkeyFolder = 'Mr Nilson';
    tallFile = 'Tall_IT_lines_N.mat';
    resp3binFile = 'SNR_capsules_N_d12.mat';
elseif P.Monkey == 2
    monkeySuffix = "F";
    monkeyFolder = 'Figaro';
    tallFile = 'Tall_IT_lines_F.mat';
    resp3binFile = 'SNR_capsules_F_d12.mat';
else
    error('Compare_IT_SessionResponses_3bin:InvalidMonkey', ...
        'Monkey must be 1 (Nilson) or 2 (Figaro).');
end

tallPath = fullfile(cfg.matDir, tallFile);
resp3binPath = fullfile(cfg.matDir, resp3binFile);
trialRespPath = fullfile(cfg.dataRoot, monkeyFolder, 'ObjAtt_lines_normMUA.mat');
trialInfoPath = fullfile(cfg.dataRoot, monkeyFolder, 'ObjAtt_lines_MUA_trials.mat');
cachePath = fullfile(cfg.resultsDir, sprintf('IT_session_split_3bin_%s.mat', char(monkeySuffix)));
outPath = fullfile(cfg.resultsDir, sprintf('Compare_IT_SessionResponses_3bin_%s.mat', char(monkeySuffix)));

assert(exist(tallPath, 'file') == 2, 'Missing %s.', tallPath);
assert(exist(resp3binPath, 'file') == 2, 'Missing %s.', resp3binPath);
assert(exist(trialRespPath, 'file') == 2, 'Missing %s.', trialRespPath);
assert(exist(trialInfoPath, 'file') == 2, 'Missing %s.', trialInfoPath);
assert(numel(P.sessions) == 2, 'This summary currently expects exactly two sessions.');

Sgeo = load(tallPath, 'Tall_IT', 'RFrange');
Tall_IT = Sgeo.Tall_IT;
RFrange = Sgeo.RFrange(:);
Rcombined = load_capsules_struct_exclusion_aware(resp3binPath, monkeySuffix, 'cfg', cfg);
timeWindows = double(Rcombined.timeWindows);
assert(size(timeWindows,1) == 3, 'Expected exactly 3 coarse time windows.');

stimNums = arrayfun(@(x) x.stimNum, Tall_IT(:));
[stimNumsSorted, ordStim] = sort(stimNums(:));
assert(all(stimNumsSorted(:).' == 1:numel(Tall_IT)), ...
    'Tall_IT.stimNum must cover 1..%d exactly.', numel(Tall_IT));
Tall_IT = Tall_IT(ordStim);
nStim = numel(Tall_IT);
nIT = numel(RFrange);
siteRows = (1:nIT).';

% Object-assignment gate from the approved IT RF.m workflow.
nObjectStim = zeros(nIT,1);
for stim = 1:nStim
    T = Tall_IT(stim).T;
    assign = string(T.assignment(siteRows));
    nObjectStim = nObjectStim + (assign == "target") + (assign == "distractor");
end
hasObjectRF = nObjectStim >= P.MinObjectStim;

% Standard IT responsiveness score from the existing 3-bin summary.
Rloc = localize_resp_to_it(Rcombined, RFrange);
SNR = compute_snr_per_color_sites(Rloc, Tall_IT, siteRows, 'Verbose', false);
topKAbsQuartetEarly = compute_topk_quartet_response(Rloc, SNR, P.TopKQuartets);
keepResponsive = hasObjectRF & isfinite(topKAbsQuartetEarly) & (topKAbsQuartetEarly > P.RespThr);

% Day-specific coarse summaries are cached after the first raw-data pass.
Split = load_or_compute_session_split(cachePath, trialRespPath, trialInfoPath, ...
    timeWindows, P, nStim);

muBySess = nan(nIT, size(timeWindows,1), numel(P.sessions));
for iSess = 1:numel(P.sessions)
    muBySess(:,:,iSess) = weighted_site_window_means(Split(iSess).meanAct(RFrange,:,:), Split(iSess).nTrials);
end

preBySess = squeeze(muBySess(:,1,:));
earlyBySess = squeeze(muBySess(:,2,:));
lateBySess = squeeze(muBySess(:,3,:));
evokedEarly = earlyBySess - preBySess;
evokedLate = lateBySess - preBySess;

exampleLocal = find(RFrange == P.exampleSiteGlobal, 1, 'first');
if isempty(exampleLocal)
    exampleLocal = NaN;
end
preShift = preBySess(:,2) - preBySess(:,1);
secondaryLocal = find_secondary_example(preShift, exampleLocal);
if isfinite(secondaryLocal)
    secondaryGlobal = RFrange(secondaryLocal);
else
    secondaryGlobal = NaN;
end

OUT = struct();
OUT.P = P;
OUT.monkeySuffix = monkeySuffix;
OUT.siteGlobal = RFrange(:);
OUT.keepResponsive = keepResponsive(:);
OUT.hasObjectRF = hasObjectRF(:);
OUT.topKAbsQuartetEarly = topKAbsQuartetEarly(:);
OUT.timeWindows = timeWindows;
OUT.sessionIds = P.sessions(:);
OUT.muBySess = muBySess;
OUT.preBySess = preBySess;
OUT.earlyBySess = earlyBySess;
OUT.lateBySess = lateBySess;
OUT.evokedEarly = evokedEarly;
OUT.evokedLate = evokedLate;
OUT.exampleSiteGlobal = P.exampleSiteGlobal;
OUT.exampleSiteLocal = exampleLocal;
OUT.secondarySiteGlobal = secondaryGlobal;
OUT.secondarySiteLocal = secondaryLocal;
OUT.summary = build_summary(preBySess, earlyBySess, lateBySess, evokedEarly, evokedLate, keepResponsive);

fprintf('IT session 3-bin comparison for monkey %s\n', char(monkeySuffix));
fprintf('All IT sites: %d | Responsive subset: %d\n', nIT, nnz(keepResponsive));
fprintf('Late evoked median (all): sess%d=%.4f, sess%d=%.4f\n', ...
    P.sessions(1), OUT.summary.all.lateEvokedMedian(1), ...
    P.sessions(2), OUT.summary.all.lateEvokedMedian(2));
fprintf('Late evoked median (responsive): sess%d=%.4f, sess%d=%.4f\n', ...
    P.sessions(1), OUT.summary.responsive.lateEvokedMedian(1), ...
    P.sessions(2), OUT.summary.responsive.lateEvokedMedian(2));
fprintf('Responsive sites with sess%d late-evoked < sess%d late-evoked: %.1f%%\n', ...
    P.sessions(2), P.sessions(1), 100 * OUT.summary.responsive.fracLateDrop);
fprintf('Responsive sites with sess%d late-evoked <= 0 and sess%d late-evoked > 0: %.1f%%\n', ...
    P.sessions(2), P.sessions(1), 100 * OUT.summary.responsive.fracLateGone);
if isfinite(exampleLocal)
    fprintf('Example site %d (local %d): late-evoked sess%d=%.4f, sess%d=%.4f\n', ...
        P.exampleSiteGlobal, exampleLocal, P.sessions(1), evokedLate(exampleLocal,1), ...
        P.sessions(2), evokedLate(exampleLocal,2));
end
if isfinite(secondaryLocal)
    fprintf('Negative pre-shift site %d (local %d): pre sess%d=%.4f, sess%d=%.4f\n', ...
        secondaryGlobal, secondaryLocal, P.sessions(1), preBySess(secondaryLocal,1), ...
        P.sessions(2), preBySess(secondaryLocal,2));
end

if P.saveResult
    save(outPath, 'OUT', '-v7.3');
    fprintf('Saved IT session comparison to %s\n', outPath);
end

if P.plotFigure
    make_summary_figure(OUT);
end
end

function Split = load_or_compute_session_split(cachePath, trialRespPath, trialInfoPath, timeWindows, P, nStim)
wantSessions = P.sessions(:).';

if P.useCache && exist(cachePath, 'file') == 2
    S = load(cachePath, 'Split', 'meta');
    if isfield(S, 'Split') && isfield(S, 'meta') ...
            && isequal(double(S.meta.timeWindows), double(timeWindows)) ...
            && isequal(double(S.meta.sessions(:).'), double(wantSessions))
        Split = S.Split;
        return;
    end
end

m1 = matfile(trialRespPath);
m2 = matfile(trialInfoPath);
Split = repmat(struct('sessionId', [], 'meanAct', [], 'meanSqAct', [], 'nTrials', [], 'stimList', []), ...
    numel(wantSessions), 1);

for iSess = 1:numel(wantSessions)
    sess = wantSessions(iSess);
    fprintf('Computing IT 3-bin session split for day %d...\n', sess);
    [meanAct, meanSqAct, nTrials, stimList] = avg_byStim( ...
        m1, m2, timeWindows, ...
        'onlyCorrect', P.onlyCorrect, ...
        'correctCol', P.correctCol, ...
        'correctVal', P.correctVal, ...
        'days', sess, ...
        'dayCol', P.dayCol, ...
        'chunkTrials', P.chunkTrials, ...
        'expectedMaxStim', P.expectedMaxStim, ...
        'verbose', true);
    assert(numel(stimList) == nStim, 'Stimulus count mismatch after day split.');
    Split(iSess).sessionId = sess;
    Split(iSess).meanAct = meanAct;
    Split(iSess).meanSqAct = meanSqAct;
    Split(iSess).nTrials = nTrials;
    Split(iSess).stimList = stimList;
end

if P.useCache
    meta = struct();
    meta.timeWindows = timeWindows;
    meta.sessions = wantSessions;
    save(cachePath, 'Split', 'meta', '-v7.3');
    fprintf('Saved IT session split cache to %s\n', cachePath);
end
end

function Rloc = localize_resp_to_it(R, RFrange)
Rloc = R;
Rloc.meanAct = R.meanAct(RFrange,:,:);
Rloc.meanSqAct = R.meanSqAct(RFrange,:,:);
if ismatrix(R.nTrials) && size(R.nTrials,1) >= max(RFrange)
    Rloc.nTrials = R.nTrials(RFrange, :);
else
    Rloc.nTrials = R.nTrials;
end
end

function topKAbsQuartetEarly = compute_topk_quartet_response(Rloc, SNR, topK)
nIT = size(Rloc.meanAct, 1);
nStim = size(Rloc.meanAct, 2);
muSpont = SNR.muSpont(:);
sdSpont = SNR.sdSpont(:);
badNoise = ~isfinite(sdSpont) | (sdSpont <= 0);

nQuartets = floor(nStim / 8) * 2;
quartetMembers = zeros(nQuartets, 4);
q = 0;
for base = 0:8:(nStim - 8)
    q = q + 1;
    quartetMembers(q,:) = base + [1 2 5 6];
    q = q + 1;
    quartetMembers(q,:) = base + [3 4 7 8];
end

if isvector(Rloc.nTrials)
    nTrialsAll = double(Rloc.nTrials(:));
    perSiteTrials = false;
else
    perSiteTrials = true;
end

muQuartetEarly = nan(nIT, nQuartets);
for iSite = 1:nIT
    if perSiteTrials
        nTrSite = double(Rloc.nTrials(iSite, :)).';
    else
        nTrSite = nTrialsAll;
    end
    rEarly = squeeze(Rloc.meanAct(iSite,:,2)).';
    for qIdx = 1:nQuartets
        stimQ = quartetMembers(qIdx,:);
        nTrQ = nTrSite(stimQ);
        rQ = rEarly(stimQ);
        idx = isfinite(rQ) & isfinite(nTrQ) & (nTrQ > 0);
        if any(idx)
            muQuartetEarly(iSite, qIdx) = sum(nTrQ(idx) .* rQ(idx)) / sum(nTrQ(idx));
        end
    end
end

signedQuartetEarly = bsxfun(@rdivide, bsxfun(@minus, muQuartetEarly, muSpont), sdSpont);
signedQuartetEarly(badNoise, :) = NaN;

topKAbsQuartetEarly = nan(nIT,1);
for iSite = 1:nIT
    vals = abs(signedQuartetEarly(iSite,:));
    vals = vals(isfinite(vals));
    if isempty(vals)
        continue;
    end
    vals = sort(vals, 'descend');
    k = min(topK, numel(vals));
    topKAbsQuartetEarly(iSite) = mean(vals(1:k));
end
end

function mu = weighted_site_window_means(meanAct, nTrials)
nIT = size(meanAct, 1);
nStim = size(meanAct, 2);
nWin = size(meanAct, 3);
assert(numel(nTrials) == nStim, 'nTrials must match stimulus dimension.');
w = double(nTrials(:));
mu = nan(nIT, nWin);
for iSite = 1:nIT
    for wIdx = 1:nWin
        x = squeeze(meanAct(iSite,:,wIdx)).';
        idx = isfinite(x) & isfinite(w) & (w > 0);
        if any(idx)
            mu(iSite, wIdx) = sum(w(idx) .* x(idx)) / sum(w(idx));
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
exLocal = OUT.exampleSiteLocal;
secLocal = OUT.secondarySiteLocal;

fig = figure('Color', 'w', ...
    'Name', sprintf('IT session 3-bin comparison %s', char(OUT.monkeySuffix)));
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile;
scatter_with_subset(ax1, OUT.preBySess(:,1), OUT.preBySess(:,2), resp, exLocal, secLocal);
xlabel(ax1, sprintf('sess%d pre', sess1));
ylabel(ax1, sprintf('sess%d pre', sess2));
title(ax1, 'Pre window');

ax2 = nexttile;
scatter_with_subset(ax2, OUT.evokedEarly(:,1), OUT.evokedEarly(:,2), resp, exLocal, secLocal);
xlabel(ax2, sprintf('sess%d early-pre', sess1));
ylabel(ax2, sprintf('sess%d early-pre', sess2));
title(ax2, 'Early evoked');

ax3 = nexttile;
scatter_with_subset(ax3, OUT.evokedLate(:,1), OUT.evokedLate(:,2), resp, exLocal, secLocal);
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
legend(ax4, {'all IT', 'responsive IT'}, 'Location', 'best', 'Box', 'off');
grid(ax4, 'on');

annotation(fig, 'textbox', [0.05 0.955 0.90 0.04], ...
    'String', sprintf(['IT session split using the standard 3 windows [%d %d], [%d %d], [%d %d] ms. ' ...
    'Blue = standard responsive/object-related IT subset. Red marker = site %d; magenta marker = strongest negative pre-shift site %d.'], ...
    OUT.timeWindows(1,1), OUT.timeWindows(1,2), ...
    OUT.timeWindows(2,1), OUT.timeWindows(2,2), ...
    OUT.timeWindows(3,1), OUT.timeWindows(3,2), ...
    OUT.exampleSiteGlobal, OUT.secondarySiteGlobal), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'FontWeight', 'bold', 'FontSize', 10);
end

function scatter_with_subset(ax, x, y, keepMask, exampleLocal, secondaryLocal)
hold(ax, 'on');
scatter(ax, x, y, 22, [0.75 0.75 0.75], 'filled', ...
    'MarkerFaceAlpha', 0.75, 'MarkerEdgeAlpha', 0.75);
scatter(ax, x(keepMask), y(keepMask), 28, [0.20 0.45 0.85], 'filled', ...
    'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0.9);
if isfinite(exampleLocal)
    scatter(ax, x(exampleLocal), y(exampleLocal), 72, [0.85 0.20 0.20], ...
        'p', 'filled', 'MarkerEdgeColor', [1 1 1], 'LineWidth', 0.8);
end
if isfinite(secondaryLocal)
    scatter(ax, x(secondaryLocal), y(secondaryLocal), 72, [0.75 0.10 0.75], ...
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

function siteLocal = find_secondary_example(preShift, exampleLocal)
siteLocal = NaN;
ok = isfinite(preShift);
if isfinite(exampleLocal) && exampleLocal >= 1 && exampleLocal <= numel(preShift)
    ok(exampleLocal) = false;
end
if ~any(ok)
    return;
end
[~, ix] = min(preShift(ok));
idx = find(ok);
siteLocal = idx(ix);
end

function lims = finite_lims(x)
x = x(isfinite(x));
if isempty(x)
    lims = [NaN NaN];
else
    lims = [min(x) max(x)];
end
end
