function OUT = Inspect_V1V4_Cofluctuation_TrialOrder(Monkey, Puser)
% INSPECT_V1V4_COFLUCTUATION_TRIALORDER
% Show running-average V1-V4 prestimulus cofluctuation across trial order
% for a few example pairs, before and after the LOESS baseline correction.

if nargin < 1 || isempty(Monkey)
    Monkey = 1;
end
if nargin < 2 || isempty(Puser)
    Puser = struct();
end

P = struct();
P.cacheTag = "all";
P.rfMinPx = 50;
P.runningWin = 80;
P.nExamples = 4;
P.effectQuantiles = [0.15 0.45 0.75 0.95];
P.loessFrac = 0.4;
P.preWindowMs = [-200 0];
P.sessions = [1 2];
P.minTrialsPerStimSess = 2;
P.onlyCorrect = true;
P.correctCol = 9;
P.correctVal = 1;
P.stimCol = 1;
P.dayCol = 11;
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
outPath = fullfile(cfg.resultsDir, sprintf('Inspect_V1V4_Cofluctuation_TrialOrder_%s_%s.mat', char(monkeySuffix), char(P.cacheTag)));
cacheParams = local_cache_params(Monkey, P);

if P.useCache && exist(outPath, 'file') == 2
    S = load(outPath, 'OUT');
    if isfield(S, 'OUT') && local_cache_matches(S.OUT, cacheParams) && ...
            session_exclusion_cache_matches(S.OUT, monkeySuffix)
        OUT = S.OUT;
        OUT.P = P;
        if P.verbose
            fprintf('Loaded cached V1-V4 trial-order cofluctuation diagnostic from %s\n', outPath);
        end
        if P.plotFigure
            local_plot_summary(OUT);
        end
        return;
    end
end

CorrLO = NoiseCorr_V1V4_RFDistance_TargetDistractor_LOESS(Monkey, struct( ...
    'plotFigure', false, ...
    'verbose', false, ...
    'useCache', true, ...
    'loessFrac', P.loessFrac));
Main = NoiseCorr_V1V4_RFDistance_TargetDistractor(Monkey, struct( ...
    'windowIdx', [1 2 3], ...
    'plotFigure', false, ...
    'verbose', false, ...
    'useCache', true));

prePosCorr = find(CorrLO.windowIdx == 1, 1, 'first');
assert(~isempty(prePosCorr), 'Could not locate prestim window in LOESS output.');
validPairs = isfinite(CorrLO.rfDistancePx) & (CorrLO.rfDistancePx > P.rfMinPx) & ...
    isfinite(CorrLO.pairCorrCond(:, 1, prePosCorr)) & isfinite(CorrLO.pairCorrCond(:, 2, prePosCorr));
pairEffect = CorrLO.pairCorrCond(:, 2, prePosCorr) - CorrLO.pairCorrCond(:, 1, prePosCorr); % DD - TT
pairIdxSel = local_pick_quantile_pairs(find(validPairs), pairEffect(validPairs), P.effectQuantiles, P.nExamples);
assert(~isempty(pairIdxSel), 'No valid V1-V4 pairs found for the requested selection.');

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

preIdx = find(tb >= P.preWindowMs(1) & tb < P.preWindowMs(2));
assert(~isempty(preIdx), 'No prestim samples found in [%d %d).', P.preWindowMs(1), P.preWindowMs(2));

siteMapV1 = nan(max(double(Main.SiteTableV1.siteIdx)), 1);
siteMapV1(double(Main.SiteTableV1.siteIdx)) = 1:height(Main.SiteTableV1);
siteMapV4 = nan(max(double(Main.SiteTableV4.siteIdx)), 1);
siteMapV4(double(Main.SiteTableV4.siteIdx)) = 1:height(Main.SiteTableV4);

pairRowsV1 = siteMapV1(double(Main.PairTable.site1Idx(pairIdxSel)));
pairRowsV4 = siteMapV4(double(Main.PairTable.site2Idx(pairIdxSel)));
assert(all(isfinite(pairRowsV1)) && all(isfinite(pairRowsV4)), 'Could not map example pairs to site rows.');

siteRowsV1 = unique(pairRowsV1(:)');
siteRowsV4 = unique(pairRowsV4(:)');
respV1Raw = local_load_prestim_subset(m1, Main.SiteTableV1.globalSite(siteRowsV1), preIdx, trialIdx);
respV4Raw = local_load_prestim_subset(m1, Main.SiteTableV4.globalSite(siteRowsV4), preIdx, trialIdx);

[respV1Corr, driftV1] = local_apply_loess(respV1Raw, dayIncl, P.sessions, P.loessFrac);
[respV4Corr, driftV4] = local_apply_loess(respV4Raw, dayIncl, P.sessions, P.loessFrac);

siteKeepSessV1 = local_site_keep(Main.SiteTableV1.globalSite(siteRowsV1), monkeySuffix, P.sessions);
siteKeepSessV4 = local_site_keep(Main.SiteTableV4.globalSite(siteRowsV4), monkeySuffix, P.sessions);
[resRawV1, ~] = local_compute_residuals(respV1Raw, stimIncl, dayIncl, P.sessions, 384, P.minTrialsPerStimSess, siteKeepSessV1);
[resRawV4, ~] = local_compute_residuals(respV4Raw, stimIncl, dayIncl, P.sessions, 384, P.minTrialsPerStimSess, siteKeepSessV4);
[resCorV1, ~] = local_compute_residuals(respV1Corr, stimIncl, dayIncl, P.sessions, 384, P.minTrialsPerStimSess, siteKeepSessV1);
[resCorV4, ~] = local_compute_residuals(respV4Corr, stimIncl, dayIncl, P.sessions, 384, P.minTrialsPerStimSess, siteKeepSessV4);

sessionCounts = arrayfun(@(s) nnz(dayIncl == s), P.sessions);
sessionEnds = cumsum(sessionCounts);

Examples = repmat(struct( ...
    'pairIdx', [], 'pairGlobalV1', [], 'pairGlobalV4', [], 'effectDDminusTT', [], ...
    'cofluctRaw', [], 'cofluctCorr', [], 'runRaw', [], 'runCorr', [], ...
    'stimType', [], 'driftSdV1', [], 'driftSdV4', []), numel(pairIdxSel), 1);

for i = 1:numel(pairIdxSel)
    p = pairIdxSel(i);
    rowV1 = pairRowsV1(i);
    rowV4 = pairRowsV4(i);
    locV1 = find(siteRowsV1 == rowV1, 1, 'first');
    locV4 = find(siteRowsV4 == rowV4, 1, 'first');

    code1 = Main.assignCodeV1(rowV1, :);
    code2 = Main.assignCodeV4(rowV4, :);
    stimType = zeros(numel(stimIncl), 1, 'uint8');
    stimType((code1(stimIncl) == 1) & (code2(stimIncl) == 1)) = 1; % TT
    stimType((code1(stimIncl) == 2) & (code2(stimIncl) == 2)) = 2; % DD
    stimType(((code1(stimIncl) == 1) & (code2(stimIncl) == 2)) | ((code1(stimIncl) == 2) & (code2(stimIncl) == 1))) = 3; % TD

    cofluctRaw = reshape(resRawV1(locV1, :), [], 1) .* reshape(resRawV4(locV4, :), [], 1);
    cofluctCorr = reshape(resCorV1(locV1, :), [], 1) .* reshape(resCorV4(locV4, :), [], 1);
    runRaw = local_running_by_session(cofluctRaw, dayIncl, P.sessions, P.runningWin);
    runCorr = local_running_by_session(cofluctCorr, dayIncl, P.sessions, P.runningWin);

    Examples(i).pairIdx = p;
    Examples(i).pairGlobalV1 = double(Main.SiteTableV1.globalSite(rowV1));
    Examples(i).pairGlobalV4 = double(Main.SiteTableV4.globalSite(rowV4));
    Examples(i).effectDDminusTT = double(pairEffect(p));
    Examples(i).cofluctRaw = cofluctRaw;
    Examples(i).cofluctCorr = cofluctCorr;
    Examples(i).runRaw = runRaw;
    Examples(i).runCorr = runCorr;
    Examples(i).stimType = stimType;
    Examples(i).driftSdV1 = std(reshape(driftV1(locV1, :), [], 1), 0, 'omitnan');
    Examples(i).driftSdV4 = std(reshape(driftV4(locV4, :), [], 1), 0, 'omitnan');
end

OUT = struct();
OUT.P = P;
OUT.cacheParams = cacheParams;
OUT.siteSessionExclusions = site_session_exclusions(monkeySuffix);
OUT.monkeySuffix = monkeySuffix;
OUT.dayIncl = dayIncl;
OUT.stimIncl = stimIncl;
OUT.sessionEnds = sessionEnds;
OUT.examplePairs = Examples;

if P.saveResult
    save(outPath, 'OUT', '-v7.3');
    if P.verbose
        fprintf('Saved V1-V4 trial-order cofluctuation diagnostic to %s\n', outPath);
    end
end

if P.plotFigure
    local_plot_summary(OUT);
end
end

function resp = local_load_prestim_subset(m1, siteGlobalVec, preIdx, trialIdx)
siteGlobalVec = double(siteGlobalVec(:));
siteBlock = min(siteGlobalVec):max(siteGlobalVec);
X = double(m1.normMUA(siteBlock, :, preIdx));
X = squeeze(mean(X, 3, 'omitnan'));
if isvector(X)
    X = X(:);
end
respAll = X(siteGlobalVec - siteBlock(1) + 1, :);
resp = respAll(:, trialIdx);
end

function siteKeepSess = local_site_keep(siteGlobalVec, monkeySuffix, sessions)
siteGlobalVec = double(siteGlobalVec(:));
siteKeepSess = true(numel(siteGlobalVec), numel(sessions));
Texcl = site_session_exclusions(monkeySuffix);
for i = 1:height(Texcl)
    row = find(siteGlobalVec == double(Texcl.siteGlobal(i)), 1, 'first');
    k = find(sessions == double(Texcl.day(i)), 1, 'first');
    if ~isempty(row) && ~isempty(k)
        siteKeepSess(row, k) = false;
    end
end
end

function [respCorr, drift] = local_apply_loess(respRaw, dayIncl, sessions, frac)
respCorr = respRaw;
nSites = size(respRaw, 1);
nTrials = size(respRaw, 2);
drift = nan(nSites, nTrials);
for s = 1:nSites
    for k = 1:numel(sessions)
        idx = find(dayIncl == sessions(k));
        if numel(idx) < 5
            continue;
        end
        y = reshape(respRaw(s, idx), [], 1);
        good = isfinite(y);
        if nnz(good) < 5
            continue;
        end
        spanN = min(nnz(good), max(5, ceil(frac * nnz(good))));
        if spanN >= nnz(good)
            fitVals = repmat(mean(y(good), 'omitnan'), nnz(good), 1);
        else
            fitVals = smoothdata(y(good), 'loess', spanN);
        end
        d = nan(size(y));
        d(good) = fitVals;
        drift(s, idx) = d;
        y(good) = y(good) - d(good);
        respCorr(s, idx) = y;
    end
end
end

function [ResZ, nTrialsSiteStimSess] = local_compute_residuals(respRaw, stimIncl, dayIncl, sessions, nStim, minTrialsPerStimSess, siteKeepSess)
nSites = size(respRaw, 1);
nTrialsIncl = size(respRaw, 2);
nSess = numel(sessions);
Res = nan(nSites, nTrialsIncl);
ResZ = nan(nSites, nTrialsIncl);
nTrialsSiteStimSess = zeros(nSites, nStim, nSess, 'uint16');

for s = 1:nSites
    for k = 1:nSess
        if ~siteKeepSess(s, k)
            continue;
        end
        mSess = (dayIncl == sessions(k));
        idxSess = find(mSess);
        for stim = 1:nStim
            idxStim = find(mSess & (stimIncl == stim));
            if numel(idxStim) < minTrialsPerStimSess
                continue;
            end
            x = reshape(respRaw(s, idxStim), [], 1);
            good = isfinite(x);
            if nnz(good) < minTrialsPerStimSess
                continue;
            end
            mu = mean(x(good));
            Res(s, idxStim(good)) = x(good) - mu;
            nTrialsSiteStimSess(s, stim, k) = uint16(nnz(good));
        end
        r = reshape(Res(s, idxSess), [], 1);
        good = isfinite(r);
        if nnz(good) < 2
            continue;
        end
        sig = std(r(good), 0);
        if ~isfinite(sig) || sig <= 0
            continue;
        end
        z = reshape(Res(s, idxSess), [], 1);
        goodZ = isfinite(z);
        z(goodZ) = z(goodZ) ./ sig;
        ResZ(s, idxSess) = z;
    end
end
end

function run = local_running_by_session(x, dayIncl, sessions, win)
run = nan(size(x));
for k = 1:numel(sessions)
    idx = find(dayIncl == sessions(k));
    if isempty(idx)
        continue;
    end
    y = reshape(x(idx), [], 1);
    run(idx) = movmean(y, win, 'omitnan');
end
end

function picked = local_pick_quantile_pairs(pairIdx, effectVals, quantiles, nExamples)
effectVals = effectVals(:);
pairIdx = pairIdx(:);
[valsSorted, ord] = sort(effectVals, 'ascend');
pairSorted = pairIdx(ord);
q = quantiles(:)';
picked = zeros(0, 1);
for i = 1:numel(q)
    pos = max(1, min(numel(pairSorted), round(q(i) * numel(pairSorted))));
    candidate = pairSorted(pos);
    if ~ismember(candidate, picked)
        picked(end+1,1) = candidate; %#ok<AGROW>
    end
end
if numel(picked) < nExamples
    for i = 1:numel(pairSorted)
        if ~ismember(pairSorted(i), picked)
            picked(end+1,1) = pairSorted(i); %#ok<AGROW>
        end
        if numel(picked) >= nExamples
            break;
        end
    end
end
picked = picked(1:min(nExamples, numel(picked)));
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
S.runningWin = double(P.runningWin);
S.nExamples = double(P.nExamples);
S.effectQuantiles = double(P.effectQuantiles(:)');
S.loessFrac = double(P.loessFrac);
S.preWindowMs = double(P.preWindowMs);
S.sessions = double(P.sessions(:)');
S.minTrialsPerStimSess = double(P.minTrialsPerStimSess);
S.onlyCorrect = logical(P.onlyCorrect);
end

function tf = local_cache_matches(OUT, cacheParams)
tf = isstruct(OUT) && isfield(OUT, 'cacheParams') && isequaln(OUT.cacheParams, cacheParams);
end

function local_plot_summary(OUT)
figure('Color', 'w', 'Name', sprintf('V1-V4 cofluctuation across trial order (%s)', char(OUT.monkeySuffix)));
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

for i = 1:numel(OUT.examplePairs)
    ex = OUT.examplePairs(i);
    ax = nexttile;
    hold(ax, 'on');
    x = 1:numel(ex.runRaw);
    plot(ax, x, ex.runRaw, '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.2, 'DisplayName', 'raw');
    plot(ax, x, ex.runCorr, '-', 'Color', [0.15 0.45 0.85], 'LineWidth', 1.8, 'DisplayName', 'LOESS 0.4');
    yMin = min([ex.runRaw(:); ex.runCorr(:)], [], 'omitnan');
    yMax = max([ex.runRaw(:); ex.runCorr(:)], [], 'omitnan');
    if ~isfinite(yMin) || ~isfinite(yMax) || yMin == yMax
        yMin = -0.1;
        yMax = 0.1;
    end
    yPad = 0.08 * (yMax - yMin + eps);
    tt = ex.stimType == 1;
    dd = ex.stimType == 2;
    plot(ax, find(tt), repmat(yMin - yPad, nnz(tt), 1), '.', 'Color', [0.20 0.60 0.30], 'MarkerSize', 5, 'DisplayName', 'TT trials');
    plot(ax, find(dd), repmat(yMin - 2*yPad, nnz(dd), 1), '.', 'Color', [0.85 0.40 0.15], 'MarkerSize', 5, 'DisplayName', 'DD trials');
    for b = 1:numel(OUT.sessionEnds)-1
        xline(ax, OUT.sessionEnds(b), 'k:', 'HandleVisibility', 'off');
    end
    xlabel(ax, sprintf('Included trial order (running mean %d)', OUT.P.runningWin));
    ylabel(ax, 'Prestim cofluctuation');
    title(ax, sprintf('V1 %d - V4 %d | DD-TT=%.3f', ...
        ex.pairGlobalV1, ex.pairGlobalV4, ex.effectDDminusTT));
    subtitle(ax, sprintf('drift SD V1=%.3f | V4=%.3f', ex.driftSdV1, ex.driftSdV4));
    grid(ax, 'on');
    if i == 1
        legend(ax, 'Location', 'best', 'Box', 'off');
    end
end
end
