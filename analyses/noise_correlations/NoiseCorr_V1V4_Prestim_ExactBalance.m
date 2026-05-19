function OUT = NoiseCorr_V1V4_Prestim_ExactBalance(Monkey, Puser)
% NOISECORR_V1V4_PRESTIM_EXACTBALANCE
% Raw-trial prestimulus control for the V1-V4 TT vs DD comparison.
%
% This diagnostic revisits the prestimulus window using the original
% normMUA trials. For each V1-V4 pair, session, and 8-stimulus block, it:
%   1) computes stimulus-conditioned residuals,
%   2) z-scores residuals within site x session,
%   3) subsamples TT and DD to the same number of valid paired trials,
%   4) computes separate TT and DD noise correlations,
%   5) averages equally across blocks and then equally across sessions.
%
% The goal is to test whether the prestim TT vs DD difference survives
% exact within-session trial balancing.

if nargin < 1 || isempty(Monkey)
    Monkey = 1;
end
if nargin < 2 || isempty(Puser)
    Puser = struct();
end

P = struct();
P.cacheTag = "all";
P.preWindowMs = [-200 0];
P.rfMinPx = 50;
P.minTrialsPerStimSess = 2;
P.minTrialsPerBlockCondSess = 2;
P.sessions = [1 2];
P.onlyCorrect = true;
P.correctCol = 9;
P.correctVal = 1;
P.stimCol = 1;
P.dayCol = 11;
P.balanceSeed = 1;
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
outPath = fullfile(cfg.resultsDir, sprintf('NoiseCorr_V1V4_Prestim_ExactBalance_%s_%s.mat', char(monkeySuffix), char(P.cacheTag)));
cacheParams = local_cache_params(Monkey, P);

if P.useCache && exist(outPath, 'file') == 2
    S = load(outPath, 'OUT');
    if isfield(S, 'OUT') && local_cache_matches(S.OUT, cacheParams) && ...
            session_exclusion_cache_matches(S.OUT, monkeySuffix)
        OUT = S.OUT;
        OUT.P = P;
        if P.verbose
            fprintf('Loaded cached V1-V4 exact-balance prestim control from %s\n', outPath);
        end
        if P.plotFigure
            local_plot_summary(OUT);
        end
        return;
    end
end

Main = NoiseCorr_V1V4_RFDistance_TargetDistractor(Monkey, struct('windowIdx', [1 2 3], 'plotFigure', false, 'verbose', false, 'useCache', true));

trialRespPath = fullfile(cfg.dataRoot, monkeyFolder, 'ObjAtt_lines_normMUA.mat');
trialInfoPath = fullfile(cfg.dataRoot, monkeyFolder, 'ObjAtt_lines_MUA_trials.mat');
assert(exist(trialRespPath, 'file') == 2, 'Missing %s.', trialRespPath);
assert(exist(trialInfoPath, 'file') == 2, 'Missing %s.', trialInfoPath);

m1 = matfile(trialRespPath);
m2 = matfile(trialInfoPath);
ALLMAT = m2.ALLMAT;
tb = double(m2.tb);
tb = tb(:);

stimPerTrial = double(ALLMAT(:, P.stimCol));
trialInclude = isfinite(stimPerTrial) & stimPerTrial >= 1 & stimPerTrial <= 384 & ...
    (floor(stimPerTrial) == stimPerTrial);
if P.onlyCorrect
    trialInclude = trialInclude & (double(ALLMAT(:, P.correctCol)) == P.correctVal);
end
trialInclude = trialInclude & ismember(double(ALLMAT(:, P.dayCol)), P.sessions(:));

trialIdx = find(trialInclude);
stimIncl = double(stimPerTrial(trialIdx));
dayIncl = double(ALLMAT(trialIdx, P.dayCol));
nTrialsIncl = numel(trialIdx);

preIdx = find(tb >= P.preWindowMs(1) & tb < P.preWindowMs(2));
assert(~isempty(preIdx), 'No prestim samples found in [%d %d).', P.preWindowMs(1), P.preWindowMs(2));

if P.verbose
    fprintf('V1-V4 exact-balance prestim control (%s)\n', char(monkeySuffix));
    fprintf('Using %d correct trials from sessions %s\n', nTrialsIncl, mat2str(P.sessions));
end

[respV1, siteKeepSessV1] = local_load_prestim_sites(m1, Main.SiteTableV1, Main.SiteTableV1.globalSite, ...
    Main.SiteTableV1.areaLocalSite, preIdx, trialIdx, monkeySuffix, P.sessions);
[respV4, siteKeepSessV4] = local_load_prestim_sites(m1, Main.SiteTableV4, Main.SiteTableV4.globalSite, ...
    Main.SiteTableV4.areaLocalSite, preIdx, trialIdx, monkeySuffix, P.sessions);

trialIdxByStimSess = cell(384, numel(P.sessions));
for k = 1:numel(P.sessions)
    mSess = dayIncl == P.sessions(k);
    for stim = 1:384
        trialIdxByStimSess{stim, k} = find(mSess & (stimIncl == stim));
    end
end

[resZV1, nTrialsV1] = local_compute_residuals(respV1, stimIncl, dayIncl, P.sessions, 384, P.minTrialsPerStimSess, siteKeepSessV1);
[resZV4, nTrialsV4] = local_compute_residuals(respV4, stimIncl, dayIncl, P.sessions, 384, P.minTrialsPerStimSess, siteKeepSessV4);

pairMask = isfinite(Main.rfDistancePx) & (Main.rfDistancePx > P.rfMinPx);
pairIdxUse = find(pairMask);
nPairsUse = numel(pairIdxUse);

pairCorrEqual = nan(nPairsUse, 2);      % TT, DD combined across sessions
pairCorrSess = nan(nPairsUse, 2, numel(P.sessions));
pairEffTrials = zeros(nPairsUse, 2);
pairEffTrialsSess = zeros(nPairsUse, 2, numel(P.sessions));

stream = RandStream('mt19937ar', 'Seed', P.balanceSeed);

for kp = 1:nPairsUse
    p = pairIdxUse(kp);
    r1 = Main.PairTable.site1Idx(p);
    r2 = Main.PairTable.site2Idx(p);

    rowV1 = find(double(Main.SiteTableV1.siteIdx) == double(r1), 1, 'first');
    rowV4 = find(double(Main.SiteTableV4.siteIdx) == double(r2), 1, 'first');
    assert(~isempty(rowV1) && ~isempty(rowV4), 'Could not map pair %d to V1/V4 rows.', p);

    code1 = Main.assignCodeV1(rowV1, :);
    code2 = Main.assignCodeV4(rowV4, :);

    zSess = nan(2, numel(P.sessions));   % TT, DD
    nSess = zeros(2, numel(P.sessions));

    for k = 1:numel(P.sessions)
        if ~siteKeepSessV1(rowV1, k) || ~siteKeepSessV4(rowV4, k)
            continue;
        end

        zBlocksTT = [];
        zBlocksDD = [];
        nEff = 0;
        for b = 1:48
            stimBlock = (b - 1) * 8 + (1:8);
            ttStim = stimBlock((code1(stimBlock) == 1) & (code2(stimBlock) == 1));
            ddStim = stimBlock((code1(stimBlock) == 2) & (code2(stimBlock) == 2));
            if isempty(ttStim) || isempty(ddStim)
                continue;
            end

            idxTT = local_collect_trials(ttStim, trialIdxByStimSess, k);
            idxDD = local_collect_trials(ddStim, trialIdxByStimSess, k);
            goodTT = isfinite(resZV1(rowV1, idxTT)) & isfinite(resZV4(rowV4, idxTT));
            goodDD = isfinite(resZV1(rowV1, idxDD)) & isfinite(resZV4(rowV4, idxDD));
            idxTT = idxTT(goodTT);
            idxDD = idxDD(goodDD);
            nUse = min(numel(idxTT), numel(idxDD));
            if nUse < P.minTrialsPerBlockCondSess
                continue;
            end

            if numel(idxTT) > nUse
                idxTT = idxTT(randperm(stream, numel(idxTT), nUse));
            end
            if numel(idxDD) > nUse
                idxDD = idxDD(randperm(stream, numel(idxDD), nUse));
            end

            rTT = local_safe_corr(resZV1(rowV1, idxTT), resZV4(rowV4, idxTT));
            rDD = local_safe_corr(resZV1(rowV1, idxDD), resZV4(rowV4, idxDD));
            if isfinite(rTT) && isfinite(rDD)
                zBlocksTT(end+1,1) = atanh(max(min(rTT, 0.999999), -0.999999)); %#ok<AGROW>
                zBlocksDD(end+1,1) = atanh(max(min(rDD, 0.999999), -0.999999)); %#ok<AGROW>
                nEff = nEff + nUse;
            end
        end

        if ~isempty(zBlocksTT) && ~isempty(zBlocksDD)
            zSess(1, k) = mean(zBlocksTT, 'omitnan');
            zSess(2, k) = mean(zBlocksDD, 'omitnan');
            nSess(:, k) = nEff;
            pairCorrSess(kp, :, k) = tanh([zSess(1, k), zSess(2, k)]);
            pairEffTrialsSess(kp, :, k) = nEff;
        end
    end

    for c = 1:2
        validSess = isfinite(zSess(c, :)) & (nSess(c, :) > 0);
        if any(validSess)
            pairCorrEqual(kp, c) = tanh(mean(zSess(c, validSess), 2, 'omitnan'));
            pairEffTrials(kp, c) = sum(nSess(c, validSess), 'all', 'omitnan');
        end
    end

    if P.verbose && (mod(kp, 500) == 0 || kp == nPairsUse)
        fprintf('  Balanced prestim pairs: %d / %d\n', kp, nPairsUse);
    end
end

meanExact = nan(2, 1);
semExact = nan(2, 1);
nPairsExact = zeros(2, 1);
for c = 1:2
    v = isfinite(pairCorrEqual(:, c));
    vals = pairCorrEqual(v, c);
    nPairsExact(c) = numel(vals);
    if isempty(vals)
        continue;
    end
    meanExact(c) = mean(vals, 'omitnan');
    if numel(vals) >= 2
        semExact(c) = std(vals, 0, 'omitnan') / sqrt(numel(vals));
    else
        semExact(c) = 0;
    end
end

meanExactSess = nan(2, numel(P.sessions));
semExactSess = nan(2, numel(P.sessions));
nPairsExactSess = zeros(2, numel(P.sessions));
for k = 1:numel(P.sessions)
    for c = 1:2
        vals = pairCorrSess(isfinite(pairCorrSess(:, c, k)), c, k);
        vals = vals(:);
        nPairsExactSess(c, k) = numel(vals);
        if isempty(vals)
            continue;
        end
        meanExactSess(c, k) = mean(vals, 'omitnan');
        if numel(vals) >= 2
            semExactSess(c, k) = std(vals, 0, 'omitnan') / sqrt(numel(vals));
        else
            semExactSess(c, k) = 0;
        end
    end
end

OUT = struct();
OUT.P = P;
OUT.cacheParams = cacheParams;
OUT.siteSessionExclusions = site_session_exclusions(monkeySuffix);
OUT.monkeySuffix = monkeySuffix;
OUT.mainPath = fullfile(cfg.resultsDir, sprintf('NoiseCorr_V1V4_RFDistance_TargetDistractor_%s_%s.mat', char(monkeySuffix), char(P.cacheTag)));
OUT.mainSummary = struct();
OUT.mainSummary.condNames = ["T-T"; "D-D"];
OUT.mainSummary.meanAboveMinRf = Main.meanAboveMinRf(1:2, 1);
OUT.mainSummary.semAboveMinRf = Main.semAboveMinRf(1:2, 1);
OUT.mainSummary.nPairsAboveMinRf = Main.nPairsAboveMinRf(1:2, 1);
OUT.mainSummary.trialWeightedMeanAboveMinRf = [ ...
    mean(Main.pairCorrCondTrialWeighted(isfinite(Main.rfDistancePx) & Main.rfDistancePx > P.rfMinPx & isfinite(Main.pairCorrCondTrialWeighted(:,1,1)), 1, 1), 'omitnan'); ...
    mean(Main.pairCorrCondTrialWeighted(isfinite(Main.rfDistancePx) & Main.rfDistancePx > P.rfMinPx & isfinite(Main.pairCorrCondTrialWeighted(:,2,1)), 2, 1), 'omitnan')];
OUT.rfMinPx = P.rfMinPx;
OUT.pairIdxUse = pairIdxUse;
OUT.pairCorrEqual = pairCorrEqual;
OUT.pairCorrSess = pairCorrSess;
OUT.pairEffTrials = pairEffTrials;
OUT.pairEffTrialsSess = pairEffTrialsSess;
OUT.meanExact = meanExact;
OUT.semExact = semExact;
OUT.nPairsExact = nPairsExact;
OUT.meanExactSess = meanExactSess;
OUT.semExactSess = semExactSess;
OUT.nPairsExactSess = nPairsExactSess;
OUT.condNames = ["T-T"; "D-D"];

if P.saveResult
    save(outPath, 'OUT', '-v7.3');
    if P.verbose
        fprintf('Saved V1-V4 exact-balance prestim control to %s\n', outPath);
    end
end

if P.plotFigure
    local_plot_summary(OUT);
end
end

function [respIncl, siteKeepSess] = local_load_prestim_sites(m1, SiteTable, siteGlobal, areaLocalSite, preIdx, trialIdx, monkeySuffix, sessions)
siteGlobal = double(siteGlobal(:));
siteBlock = min(siteGlobal):max(siteGlobal);
Xblock = double(m1.normMUA(siteBlock, :, preIdx));
Xblock = squeeze(mean(Xblock, 3, 'omitnan'));
if isvector(Xblock)
    Xblock = Xblock(:);
end
respAll = Xblock(siteGlobal - siteBlock(1) + 1, :);
respIncl = respAll(:, trialIdx);

siteKeepSess = true(height(SiteTable), numel(sessions));
Texcl = site_session_exclusions(monkeySuffix);
for i = 1:height(Texcl)
    row = find(siteGlobal == double(Texcl.siteGlobal(i)), 1, 'first');
    k = find(sessions == double(Texcl.day(i)), 1, 'first');
    if ~isempty(row) && ~isempty(k)
        siteKeepSess(row, k) = false;
    end
end

assert(numel(areaLocalSite) == height(SiteTable), 'Site table/local-site mismatch.');
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

        r = reshape(Res(s, mSess), [], 1);
        good = isfinite(r);
        if nnz(good) < 2
            continue;
        end
        sig = std(r(good), 0);
        if ~isfinite(sig) || sig <= 0
            continue;
        end
        idxUse = find(mSess);
        z = reshape(Res(s, idxUse), [], 1);
        goodZ = isfinite(z);
        z(goodZ) = z(goodZ) ./ sig;
        ResZ(s, idxUse) = z;
    end
end
end

function idx = local_collect_trials(stimList, trialIdxByStimSess, k)
idx = zeros(0, 1);
for s = 1:numel(stimList)
    idx = [idx; trialIdxByStimSess{stimList(s), k}(:)]; %#ok<AGROW>
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
r = R(1,2);
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
S.preWindowMs = double(P.preWindowMs);
S.rfMinPx = double(P.rfMinPx);
S.minTrialsPerStimSess = double(P.minTrialsPerStimSess);
S.minTrialsPerBlockCondSess = double(P.minTrialsPerBlockCondSess);
S.sessions = double(P.sessions(:)');
S.onlyCorrect = logical(P.onlyCorrect);
S.correctCol = double(P.correctCol);
S.correctVal = double(P.correctVal);
S.stimCol = double(P.stimCol);
S.dayCol = double(P.dayCol);
S.balanceSeed = double(P.balanceSeed);
end

function tf = local_cache_matches(OUT, cacheParams)
tf = isstruct(OUT) && isfield(OUT, 'cacheParams') && isequaln(OUT.cacheParams, cacheParams);
end

function local_plot_summary(OUT)
figure('Color', 'w', 'Name', sprintf('V1-V4 exact-balance prestim control (%s)', char(OUT.monkeySuffix)));
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

colors = [0.15 0.45 0.85; 0.85 0.40 0.15];
labels = cellstr(OUT.condNames);

ax1 = nexttile;
hold(ax1, 'on');
x = 1:2;
y = OUT.mainSummary.meanAboveMinRf(:);
e = OUT.mainSummary.semAboveMinRf(:);
b = bar(ax1, x, y, 'FaceColor', 'flat', 'EdgeColor', 'none');
b.CData = colors;
errorbar(ax1, x, y, e, 'k.', 'LineWidth', 1.0, 'HandleVisibility', 'off');
for i = 1:2
    text(ax1, x(i), y(i), sprintf('  n=%d', OUT.mainSummary.nPairsAboveMinRf(i)), ...
        'Rotation', 90, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', ...
        'FontSize', 8, 'Color', [0.2 0.2 0.2]);
end
set(ax1, 'XTick', x, 'XTickLabel', labels);
ylabel(ax1, 'Mean noise correlation');
title(ax1, sprintf('Current main estimate | RF distance > %d px', OUT.rfMinPx));
grid(ax1, 'on');

ax2 = nexttile;
hold(ax2, 'on');
y = OUT.meanExact(:);
e = OUT.semExact(:);
b = bar(ax2, x, y, 'FaceColor', 'flat', 'EdgeColor', 'none');
b.CData = colors;
errorbar(ax2, x, y, e, 'k.', 'LineWidth', 1.0, 'HandleVisibility', 'off');
for i = 1:2
    text(ax2, x(i), y(i), sprintf('  n=%d', OUT.nPairsExact(i)), ...
        'Rotation', 90, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', ...
        'FontSize', 8, 'Color', [0.2 0.2 0.2]);
end
set(ax2, 'XTick', x, 'XTickLabel', labels);
ylabel(ax2, 'Mean noise correlation');
title(ax2, sprintf('Exact balanced raw-trial control | RF distance > %d px', OUT.rfMinPx));
grid(ax2, 'on');
end
