function OUT = NoiseCorr_V1V4_Prestim_HistoryModel(Monkey, Puser)
% NOISECORR_V1V4_PRESTIM_HISTORYMODEL
% Test whether previous-trial shared activity explains the current
% prestimulus TT vs DD difference for V1-V4 pairs.
%
% For each V1-V4 pair (RF distance > threshold), this diagnostic:
%   1) computes stimulus-conditioned z-scored residuals for pre/early/late,
%   2) forms current prestimulus cofluctuation x_pre(t)*y_pre(t),
%   3) links each current included trial to the previous included trial,
%   4) fits per-pair models:
%        model 0: currPre ~ 1 + currDD
%        model 1: currPre ~ 1 + currDD + prevPre + prevEarly + prevLate
%   5) compares the TT-vs-DD coefficient before and after history terms.

if nargin < 1 || isempty(Monkey)
    Monkey = 1;
end
if nargin < 2 || isempty(Puser)
    Puser = struct();
end

P = struct();
P.cacheTag = "all";
P.rfMinPx = 50;
P.minTrialsPerStimSess = 2;
P.minObsPerPair = 20;
P.minObsPerCond = 10;
P.sessions = [1 2];
P.onlyCorrect = true;
P.correctCol = 9;
P.correctVal = 1;
P.stimCol = 1;
P.dayCol = 11;
P.preWindowMs = [-200 0];
P.earlyWindowMs = [40 240];
P.lateWindowMs = [300 500];
P.useCache = true;
P.saveResult = true;
P.plotFigure = true;
P.verbose = true;

if ~isempty(Puser)
    fn = fieldnames(Puser);
    for i = 1:numel(fn)
        P.(fn{i}) = Puser.(fn{i});
    end
end

cfg = config();
[monkeySuffix, monkeyFolder] = local_monkey_info(Monkey);
outPath = fullfile(cfg.resultsDir, sprintf('NoiseCorr_V1V4_Prestim_HistoryModel_%s_%s.mat', char(monkeySuffix), char(P.cacheTag)));
cacheParams = local_cache_params(Monkey, P);

if P.useCache && exist(outPath, 'file') == 2
    S = load(outPath, 'OUT');
    if isfield(S, 'OUT') && local_cache_matches(S.OUT, cacheParams) && ...
            session_exclusion_cache_matches(S.OUT, monkeySuffix)
        OUT = S.OUT;
        OUT.P = P;
        if P.verbose
            fprintf('Loaded cached V1-V4 prestim history model from %s\n', outPath);
        end
        if P.plotFigure
            local_plot_summary(OUT);
        end
        return;
    end
end

Main = NoiseCorr_V1V4_RFDistance_TargetDistractor(Monkey, struct('windowIdx', [1 2 3], ...
    'plotFigure', false, 'verbose', false, 'useCache', true));

trialRespPath = fullfile(cfg.dataRoot, monkeyFolder, 'ObjAtt_lines_normMUA.mat');
trialInfoPath = fullfile(cfg.dataRoot, monkeyFolder, 'ObjAtt_lines_MUA_trials.mat');
assert(exist(trialRespPath, 'file') == 2, 'Missing %s.', trialRespPath);
assert(exist(trialInfoPath, 'file') == 2, 'Missing %s.', trialInfoPath);

m1 = matfile(trialRespPath);
m2 = matfile(trialInfoPath);
ALLMAT = double(m2.ALLMAT);
tb = double(m2.tb);
tb = tb(:);

stimAll = double(ALLMAT(:, P.stimCol));
dayAll = double(ALLMAT(:, P.dayCol));
trialInclude = isfinite(stimAll) & stimAll >= 1 & stimAll <= 384 & ...
    (floor(stimAll) == stimAll);
if P.onlyCorrect
    trialInclude = trialInclude & (double(ALLMAT(:, P.correctCol)) == P.correctVal);
end
trialInclude = trialInclude & ismember(dayAll, P.sessions(:));

trialIdx = find(trialInclude);
stimIncl = stimAll(trialIdx);
dayIncl = dayAll(trialIdx);
nTrialsIncl = numel(trialIdx);

winDefs = [P.preWindowMs; P.earlyWindowMs; P.lateWindowMs];
winNames = ["pre"; "early"; "late"];
nWin = size(winDefs, 1);
sampleIdx = cell(nWin, 1);
for w = 1:nWin
    sampleIdx{w} = find(tb >= winDefs(w, 1) & tb < winDefs(w, 2));
    assert(~isempty(sampleIdx{w}), 'No samples found for window %s.', winNames(w));
end

if P.verbose
    fprintf('V1-V4 prestim history model (%s)\n', char(monkeySuffix));
    fprintf('Using %d correct trials from sessions %s\n', nTrialsIncl, mat2str(P.sessions));
end

[respV1, siteKeepSessV1] = local_load_window_sites(m1, Main.SiteTableV1, Main.SiteTableV1.globalSite, ...
    sampleIdx, trialIdx, monkeySuffix, P.sessions);
[respV4, siteKeepSessV4] = local_load_window_sites(m1, Main.SiteTableV4, Main.SiteTableV4.globalSite, ...
    sampleIdx, trialIdx, monkeySuffix, P.sessions);

[resZV1, nTrialsV1] = local_compute_residuals_windows(respV1, stimIncl, dayIncl, P.sessions, 384, ...
    P.minTrialsPerStimSess, siteKeepSessV1);
[resZV4, nTrialsV4] = local_compute_residuals_windows(respV4, stimIncl, dayIncl, P.sessions, 384, ...
    P.minTrialsPerStimSess, siteKeepSessV4);

prevIncluded = nan(nTrialsIncl, 1);
for k = 1:numel(P.sessions)
    rows = find(dayIncl == P.sessions(k));
    if numel(rows) >= 2
        prevIncluded(rows(2:end)) = rows(1:end-1);
    end
end

siteMapV1 = nan(max(double(Main.SiteTableV1.siteIdx)), 1);
siteMapV1(double(Main.SiteTableV1.siteIdx)) = 1:height(Main.SiteTableV1);
siteMapV4 = nan(max(double(Main.SiteTableV4.siteIdx)), 1);
siteMapV4(double(Main.SiteTableV4.siteIdx)) = 1:height(Main.SiteTableV4);

pairMask = isfinite(Main.rfDistancePx) & (Main.rfDistancePx > P.rfMinPx);
pairIdxUse = find(pairMask);
nPairsUse = numel(pairIdxUse);

betaSimple = nan(nPairsUse, 2);   % intercept, currDD
betaHistory = nan(nPairsUse, 5);  % intercept, currDD, prevPre, prevEarly, prevLate
r2Simple = nan(nPairsUse, 1);
r2History = nan(nPairsUse, 1);
histCorr = nan(nPairsUse, 3);     % corr(currPre, prevPre/early/late)
nObsPair = zeros(nPairsUse, 1);
nObsTT = zeros(nPairsUse, 1);
nObsDD = zeros(nPairsUse, 1);

for kp = 1:nPairsUse
    p = pairIdxUse(kp);
    rowV1 = siteMapV1(double(Main.PairTable.site1Idx(p)));
    rowV4 = siteMapV4(double(Main.PairTable.site2Idx(p)));
    assert(isfinite(rowV1) && isfinite(rowV4), 'Could not map pair %d to V1/V4 rows.', p);

    currCondByStim = local_pair_curr_ttdd(Main.assignCodeV1(rowV1, :), Main.assignCodeV4(rowV4, :));
    currCond = currCondByStim(stimIncl);
    prevIdx = prevIncluded;
    hasPrev = isfinite(prevIdx) & prevIdx >= 1;
    sel = hasPrev & ismember(currCond, [1 2]);
    if ~any(sel)
        continue;
    end

    idxCurr = find(sel);
    idxPrev = prevIdx(idxCurr);

    currPre = reshape(resZV1(rowV1, idxCurr, 1), [], 1) .* reshape(resZV4(rowV4, idxCurr, 1), [], 1);
    prevPre = reshape(resZV1(rowV1, idxPrev, 1), [], 1) .* reshape(resZV4(rowV4, idxPrev, 1), [], 1);
    prevEarly = reshape(resZV1(rowV1, idxPrev, 2), [], 1) .* reshape(resZV4(rowV4, idxPrev, 2), [], 1);
    prevLate = reshape(resZV1(rowV1, idxPrev, 3), [], 1) .* reshape(resZV4(rowV4, idxPrev, 3), [], 1);
    currDD = reshape(double(currCond(idxCurr) == 2), [], 1);

    good = isfinite(currPre) & isfinite(prevPre) & isfinite(prevEarly) & isfinite(prevLate) & isfinite(currDD);
    if nnz(good) < P.minObsPerPair
        continue;
    end

    currPre = currPre(good);
    prevPre = prevPre(good);
    prevEarly = prevEarly(good);
    prevLate = prevLate(good);
    currDD = currDD(good);
    nObsPair(kp) = numel(currPre);
    nObsTT(kp) = nnz(currDD == 0);
    nObsDD(kp) = nnz(currDD == 1);
    if nObsTT(kp) < P.minObsPerCond || nObsDD(kp) < P.minObsPerCond
        continue;
    end

    X0 = [ones(numel(currPre), 1), currDD];
    b0 = X0 \ currPre;
    yhat0 = X0 * b0;
    betaSimple(kp, :) = b0(:)';
    r2Simple(kp) = local_r2(currPre, yhat0);

    X1 = [ones(numel(currPre), 1), currDD, prevPre, prevEarly, prevLate];
    b1 = X1 \ currPre;
    yhat1 = X1 * b1;
    betaHistory(kp, :) = b1(:)';
    r2History(kp) = local_r2(currPre, yhat1);

    histCorr(kp, 1) = local_safe_corr(currPre, prevPre);
    histCorr(kp, 2) = local_safe_corr(currPre, prevEarly);
    histCorr(kp, 3) = local_safe_corr(currPre, prevLate);

    if P.verbose && (mod(kp, 2000) == 0 || kp == nPairsUse)
        fprintf('  History-model pairs: %d / %d\n', kp, nPairsUse);
    end
end

validPair = isfinite(betaSimple(:, 2)) & isfinite(betaHistory(:, 2));
nValidPairs = nnz(validPair);

[meanBetaSimpleDD, semBetaSimpleDD] = local_mean_sem_col(betaSimple(validPair, 2));
[meanBetaHistoryDD, semBetaHistoryDD] = local_mean_sem_col(betaHistory(validPair, 2));
[meanBetaPrev, semBetaPrev] = local_mean_sem_cols(betaHistory(validPair, 3:5));
[meanHistCorr, semHistCorr] = local_mean_sem_cols(histCorr(validPair, :));
[meanR2, semR2] = local_mean_sem_cols([r2Simple(validPair), r2History(validPair)]);

OUT = struct();
OUT.P = P;
OUT.cacheParams = cacheParams;
OUT.siteSessionExclusions = site_session_exclusions(monkeySuffix);
OUT.monkeySuffix = monkeySuffix;
OUT.mainPath = fullfile(cfg.resultsDir, sprintf('NoiseCorr_V1V4_RFDistance_TargetDistractor_%s_%s.mat', char(monkeySuffix), char(P.cacheTag)));
OUT.rfMinPx = P.rfMinPx;
OUT.pairIdxUse = pairIdxUse;
OUT.validPair = validPair;
OUT.nValidPairs = nValidPairs;
OUT.nObsPair = nObsPair;
OUT.nObsTT = nObsTT;
OUT.nObsDD = nObsDD;
OUT.betaSimple = betaSimple;
OUT.betaHistory = betaHistory;
OUT.r2Simple = r2Simple;
OUT.r2History = r2History;
OUT.histCorr = histCorr;
OUT.predictorNames = ["currDD"; "prevPre"; "prevEarly"; "prevLate"];
OUT.meanBetaSimpleDD = meanBetaSimpleDD;
OUT.semBetaSimpleDD = semBetaSimpleDD;
OUT.meanBetaHistoryDD = meanBetaHistoryDD;
OUT.semBetaHistoryDD = semBetaHistoryDD;
OUT.meanBetaPrev = meanBetaPrev;
OUT.semBetaPrev = semBetaPrev;
OUT.meanHistCorr = meanHistCorr;
OUT.semHistCorr = semHistCorr;
OUT.meanR2 = meanR2;
OUT.semR2 = semR2;

if P.saveResult
    save(outPath, 'OUT', '-v7.3');
    if P.verbose
        fprintf('Saved V1-V4 prestim history model to %s\n', outPath);
    end
end

if P.plotFigure
    local_plot_summary(OUT);
end
end

function [respIncl, siteKeepSess] = local_load_window_sites(m1, SiteTable, siteGlobal, sampleIdx, trialIdx, monkeySuffix, sessions)
siteGlobal = double(siteGlobal(:));
siteBlock = min(siteGlobal):max(siteGlobal);
nSites = numel(siteGlobal);
nTrialsIncl = numel(trialIdx);
nWin = numel(sampleIdx);
respIncl = nan(nSites, nTrialsIncl, nWin);
for w = 1:nWin
    Xblock = double(m1.normMUA(siteBlock, :, sampleIdx{w}));
    Xblock = squeeze(mean(Xblock, 3, 'omitnan'));
    if isvector(Xblock)
        Xblock = Xblock(:);
    end
    respAll = Xblock(siteGlobal - siteBlock(1) + 1, :);
    respIncl(:, :, w) = respAll(:, trialIdx);
end

siteKeepSess = true(height(SiteTable), numel(sessions));
Texcl = site_session_exclusions(monkeySuffix);
for i = 1:height(Texcl)
    row = find(siteGlobal == double(Texcl.siteGlobal(i)), 1, 'first');
    k = find(sessions == double(Texcl.day(i)), 1, 'first');
    if ~isempty(row) && ~isempty(k)
        siteKeepSess(row, k) = false;
    end
end
end

function [ResZ, nTrialsSiteStimSess] = local_compute_residuals_windows(respRaw, stimIncl, dayIncl, sessions, nStim, minTrialsPerStimSess, siteKeepSess)
nSites = size(respRaw, 1);
nTrialsIncl = size(respRaw, 2);
nWin = size(respRaw, 3);
nSess = numel(sessions);
Res = nan(nSites, nTrialsIncl, nWin);
ResZ = nan(nSites, nTrialsIncl, nWin);
nTrialsSiteStimSess = zeros(nSites, nStim, nSess, nWin, 'uint16');

for s = 1:nSites
    for k = 1:nSess
        if ~siteKeepSess(s, k)
            continue;
        end
        mSess = (dayIncl == sessions(k));
        idxSess = find(mSess);
        for w = 1:nWin
            for stim = 1:nStim
                idxStimLocal = find(mSess & (stimIncl == stim));
                if numel(idxStimLocal) < minTrialsPerStimSess
                    continue;
                end
                x = reshape(respRaw(s, idxStimLocal, w), [], 1);
                good = isfinite(x);
                if nnz(good) < minTrialsPerStimSess
                    continue;
                end
                mu = mean(x(good));
                Res(s, idxStimLocal(good), w) = x(good) - mu;
                nTrialsSiteStimSess(s, stim, k, w) = uint16(nnz(good));
            end

            r = reshape(Res(s, idxSess, w), [], 1);
            good = isfinite(r);
            if nnz(good) < 2
                continue;
            end
            sig = std(r(good), 0);
            if ~isfinite(sig) || sig <= 0
                continue;
            end
            z = reshape(Res(s, idxSess, w), [], 1);
            goodZ = isfinite(z);
            z(goodZ) = z(goodZ) ./ sig;
            ResZ(s, idxSess, w) = z;
        end
    end
end
end

function currCondByStim = local_pair_curr_ttdd(code1, code2)
currCondByStim = zeros(384, 1, 'uint8');
currCondByStim((code1 == 1) & (code2 == 1)) = 1; % TT
currCondByStim((code1 == 2) & (code2 == 2)) = 2; % DD
end

function r2 = local_r2(y, yhat)
y = y(:);
yhat = yhat(:);
good = isfinite(y) & isfinite(yhat);
if nnz(good) < 2
    r2 = nan;
    return;
end
y = y(good);
yhat = yhat(good);
ssRes = sum((y - yhat).^2);
ssTot = sum((y - mean(y)).^2);
if ssTot <= 0
    r2 = nan;
else
    r2 = 1 - ssRes / ssTot;
end
end

function r = local_safe_corr(x, y)
x = double(x(:));
y = double(y(:));
good = isfinite(x) & isfinite(y);
if nnz(good) < 2
    r = nan;
    return;
end
R = corrcoef(x(good), y(good));
r = R(1, 2);
end

function [mu, sem] = local_mean_sem_col(x)
x = x(isfinite(x));
if isempty(x)
    mu = nan;
    sem = nan;
    return;
end
mu = mean(x, 'omitnan');
if numel(x) >= 2
    sem = std(x, 0, 'omitnan') / sqrt(numel(x));
else
    sem = 0;
end
end

function [mu, sem] = local_mean_sem_cols(X)
mu = mean(X, 1, 'omitnan');
sem = std(X, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(X), 1));
end

function [monkeySuffix, monkeyFolder] = local_monkey_info(monkeyId)
switch monkeyId
    case 1
        monkeySuffix = "N";
        monkeyFolder = 'Mr Nilson';
    case 2
        monkeySuffix = "F";
        monkeyFolder = 'Figaro';
    otherwise
        error('Monkey must be 1 (Nilson) or 2 (Figaro).');
end
end

function S = local_cache_params(Monkey, P)
S = struct();
S.cacheVersion = 1;
S.Monkey = double(Monkey);
S.cacheTag = string(P.cacheTag);
S.rfMinPx = double(P.rfMinPx);
S.minTrialsPerStimSess = double(P.minTrialsPerStimSess);
S.minObsPerPair = double(P.minObsPerPair);
S.minObsPerCond = double(P.minObsPerCond);
S.sessions = double(P.sessions(:)');
S.onlyCorrect = logical(P.onlyCorrect);
S.correctCol = double(P.correctCol);
S.correctVal = double(P.correctVal);
S.stimCol = double(P.stimCol);
S.dayCol = double(P.dayCol);
S.preWindowMs = double(P.preWindowMs);
S.earlyWindowMs = double(P.earlyWindowMs);
S.lateWindowMs = double(P.lateWindowMs);
end

function tf = local_cache_matches(OUT, cacheParams)
tf = isstruct(OUT) && isfield(OUT, 'cacheParams') && isequaln(OUT.cacheParams, cacheParams);
end

function local_plot_summary(OUT)
figure('Color', 'w', 'Name', sprintf('V1-V4 prestim history model (%s)', char(OUT.monkeySuffix)));
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

ax1 = nexttile;
hold(ax1, 'on');
y = [OUT.meanBetaSimpleDD, OUT.meanBetaHistoryDD];
e = [OUT.semBetaSimpleDD, OUT.semBetaHistoryDD];
x = 1:2;
bar(ax1, x, y, 'FaceColor', [0.25 0.45 0.80], 'EdgeColor', 'none');
errorbar(ax1, x, y, e, 'k.', 'LineWidth', 1.0, 'HandleVisibility', 'off');
set(ax1, 'XTick', x, 'XTickLabel', {'currDD only', 'currDD + history'});
ylabel(ax1, 'Mean per-pair \beta_{DD}');
title(ax1, sprintf('DD effect on current prestim | n=%d pairs', OUT.nValidPairs));
grid(ax1, 'on');

ax2 = nexttile;
hold(ax2, 'on');
x = 1:3;
bar(ax2, x, OUT.meanBetaPrev, 'FaceColor', [0.80 0.45 0.20], 'EdgeColor', 'none');
errorbar(ax2, x, OUT.meanBetaPrev, OUT.semBetaPrev, 'k.', 'LineWidth', 1.0, 'HandleVisibility', 'off');
set(ax2, 'XTick', x, 'XTickLabel', {'prev pre', 'prev early', 'prev late'});
ylabel(ax2, 'Mean history coefficient');
title(ax2, 'History predictors');
grid(ax2, 'on');

ax3 = nexttile;
hold(ax3, 'on');
x = 1:3;
bar(ax3, x, OUT.meanHistCorr, 'FaceColor', [0.30 0.65 0.35], 'EdgeColor', 'none');
errorbar(ax3, x, OUT.meanHistCorr, OUT.semHistCorr, 'k.', 'LineWidth', 1.0, 'HandleVisibility', 'off');
plot(ax3, x, zeros(size(x)), 'k:', 'HandleVisibility', 'off');
set(ax3, 'XTick', x, 'XTickLabel', {'prev pre', 'prev early', 'prev late'});
ylabel(ax3, 'Corr(curr pre, prev window)');
title(ax3, 'Current-pre vs previous-window coupling');
grid(ax3, 'on');
end
