function OUT = NoiseCorr_V1V4_RFDistance_TargetDistractor_LOESS(Monkey, Puser)
% NOISECORR_V1V4_RFDISTANCE_TARGETDISTRACTOR_LOESS
% V1-V4 noise-correlation summary after subtracting a slow prestimulus
% baseline drift estimated with a LOESS curve fit within each session.
%
% The drift estimate is fit to the prestimulus baseline [-200 0] ms only,
% then subtracted from pre/early/late window responses on the same trials.

if nargin < 1 || isempty(Monkey)
    Monkey = 1;
end
if nargin < 2 || isempty(Puser)
    Puser = struct();
end

P = struct();
P.cacheTag = "all";
P.windowIdx = [1 2 3];
P.binWidthPx = 10;
P.minTrials = 15;
P.summaryMinRfPx = 50;
P.sessionCombineMode = "equal";
P.preWindowMs = [-200 0];
P.earlyWindowMs = [40 240];
P.lateWindowMs = [300 500];
P.loessFrac = 0.4;
P.loessMethod = "loess";
P.minTrialsPerStimSess = 2;
P.sessions = [1 2];
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
outPath = fullfile(cfg.resultsDir, sprintf('NoiseCorr_V1V4_RFDistance_TargetDistractor_LOESS_%s_%s.mat', char(monkeySuffix), char(P.cacheTag)));
cacheParams = local_cache_params(Monkey, P);

if P.useCache && exist(outPath, 'file') == 2
    S = load(outPath, 'OUT');
    if isfield(S, 'OUT') && local_cache_matches(S.OUT, cacheParams) && ...
            session_exclusion_cache_matches(S.OUT, monkeySuffix)
        OUT = S.OUT;
        OUT.P = P;
        if P.verbose
            fprintf('Loaded cached V1-V4 LOESS drift-corrected analysis from %s\n', outPath);
        end
        if P.plotFigure
            local_plot_summary(OUT);
            local_plot_difference_summary(OUT);
        end
        return;
    end
end

Main = NoiseCorr_V1V4_RFDistance_TargetDistractor(Monkey, struct( ...
    'windowIdx', P.windowIdx, ...
    'binWidthPx', P.binWidthPx, ...
    'minTrials', P.minTrials, ...
    'summaryMinRfPx', P.summaryMinRfPx, ...
    'sessionCombineMode', P.sessionCombineMode, ...
    'plotFigure', false, ...
    'verbose', false, ...
    'useCache', true));

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
sampleIdx = cell(3, 1);
for w = 1:3
    sampleIdx{w} = find(tb >= winDefs(w, 1) & tb < winDefs(w, 2));
    assert(~isempty(sampleIdx{w}), 'No samples found for analysis window %d.', w);
end

if P.verbose
    fprintf('V1-V4 LOESS drift-corrected analysis (%s)\n', char(monkeySuffix));
    fprintf('Using %d correct trials from sessions %s\n', nTrialsIncl, mat2str(P.sessions));
end

[respV1, siteKeepSessV1] = local_load_window_sites(m1, Main.SiteTableV1, Main.SiteTableV1.globalSite, ...
    sampleIdx, trialIdx, monkeySuffix, P.sessions);
[respV4, siteKeepSessV4] = local_load_window_sites(m1, Main.SiteTableV4, Main.SiteTableV4.globalSite, ...
    sampleIdx, trialIdx, monkeySuffix, P.sessions);

[respV1Det, driftV1] = local_apply_loess(respV1, dayIncl, P.sessions, P.loessFrac, P.loessMethod);
[respV4Det, driftV4] = local_apply_loess(respV4, dayIncl, P.sessions, P.loessFrac, P.loessMethod);

[resZV1, nTrialsV1] = local_compute_residuals_windows(respV1Det, stimIncl, dayIncl, P.sessions, 384, ...
    P.minTrialsPerStimSess, siteKeepSessV1);
[resZV4, nTrialsV4] = local_compute_residuals_windows(respV4Det, stimIncl, dayIncl, P.sessions, 384, ...
    P.minTrialsPerStimSess, siteKeepSessV4);

trialIdxByStimSess = cell(384, numel(P.sessions));
for k = 1:numel(P.sessions)
    mSess = dayIncl == P.sessions(k);
    for stim = 1:384
        trialIdxByStimSess{stim, k} = find(mSess & (stimIncl == stim));
    end
end

siteMapV1 = nan(max(double(Main.SiteTableV1.siteIdx)), 1);
siteMapV1(double(Main.SiteTableV1.siteIdx)) = 1:height(Main.SiteTableV1);
siteMapV4 = nan(max(double(Main.SiteTableV4.siteIdx)), 1);
siteMapV4(double(Main.SiteTableV4.siteIdx)) = 1:height(Main.SiteTableV4);

nPairs = height(Main.PairTable);
nWinSel = numel(Main.windowIdx);
nCond = numel(Main.condNames);
pairSite1V1 = siteMapV1(double(Main.PairTable.site1Idx));
pairSite2V4 = siteMapV4(double(Main.PairTable.site2Idx));
assert(all(isfinite(pairSite1V1)) && all(isfinite(pairSite2V4)), ...
    'Could not map all V1-V4 pairs to area-specific rows.');

pairLin = sub2ind([height(Main.SiteTableV1), height(Main.SiteTableV4)], pairSite1V1, pairSite2V4);

pairCondCode = zeros(nPairs, 384, 'uint8');
for p = 1:nPairs
    i1 = pairSite1V1(p);
    i2 = pairSite2V4(p);
    code1 = Main.assignCodeV1(i1, :);
    code2 = Main.assignCodeV4(i2, :);
    condCode = zeros(1, 384, 'uint8');
    condCode((code1 == 1) & (code2 == 1)) = 1; % T-T
    condCode((code1 == 1) & (code2 == 2)) = 2; % T-D
    condCode((code1 == 2) & (code2 == 1)) = 3; % D-T
    condCode((code1 == 2) & (code2 == 2)) = 4; % D-D
    condCode((code1 == 3) & (code2 == 3)) = 5; % B-B
    condCode((code1 == 3) & (code2 == 1)) = 6; % B-T
    condCode((code1 == 1) & (code2 == 3)) = 7; % T-B
    condCode((code1 == 3) & (code2 == 2)) = 8; % B-D
    condCode((code1 == 2) & (code2 == 3)) = 9; % D-B
    pairCondCode(p, :) = condCode;
end

sumX2V1 = zeros(height(Main.SiteTableV1), 384, numel(P.sessions), nWinSel, 'single');
sumX2V4 = zeros(height(Main.SiteTableV4), 384, numel(P.sessions), nWinSel, 'single');
sumXY = zeros(nPairs, 384, numel(P.sessions), nWinSel, 'single');

for k = 1:numel(P.sessions)
    for stim = 1:384
        idx = trialIdxByStimSess{stim, k};
        if isempty(idx)
            continue;
        end
        for wSel = 1:nWinSel
            wRaw = Main.windowIdx(wSel);
            Z1 = double(resZV1(:, idx, wRaw));
            Z4 = double(resZV4(:, idx, wRaw));
            sumX2V1(:, stim, k, wSel) = single(sum(Z1.^2, 2, 'omitnan'));
            sumX2V4(:, stim, k, wSel) = single(sum(Z4.^2, 2, 'omitnan'));
            Z1(~isfinite(Z1)) = 0;
            Z4(~isfinite(Z4)) = 0;
            C = Z1 * Z4.';
            sumXY(:, stim, k, wSel) = single(C(pairLin));
        end
    end
    if P.verbose
        fprintf('  Accumulated detrended moments for session %d / %d\n', k, numel(P.sessions));
    end
end

pairCorrCond = nan(nPairs, nCond, nWinSel);
pairCorrCondSess = nan(nPairs, nCond, nWinSel, numel(P.sessions));
pairTrialCond = zeros(nPairs, nCond, nWinSel);
pairTrialCondSess = zeros(nPairs, nCond, nWinSel, numel(P.sessions));

for p = 1:nPairs
    i1 = pairSite1V1(p);
    i2 = pairSite2V4(p);
    for c = 1:nCond
        stimMask = (pairCondCode(p, :) == c);
        if ~any(stimMask)
            continue;
        end
        for wSel = 1:nWinSel
            rSess = nan(1, numel(P.sessions));
            wSess = zeros(1, numel(P.sessions));
            for k = 1:numel(P.sessions)
                nUseSess = sum(min( ...
                    double(squeeze(nTrialsV1(i1, stimMask, k, Main.windowIdx(wSel)))), ...
                    double(squeeze(nTrialsV4(i2, stimMask, k, Main.windowIdx(wSel))))), ...
                    'all', 'omitnan');
                pairTrialCondSess(p, c, wSel, k) = nUseSess;
                if ~(isfinite(nUseSess) && nUseSess > 0)
                    continue;
                end
                num = sum(double(squeeze(sumXY(p, stimMask, k, wSel))), 'all', 'omitnan');
                den1 = sum(double(squeeze(sumX2V1(i1, stimMask, k, wSel))), 'all', 'omitnan');
                den2 = sum(double(squeeze(sumX2V4(i2, stimMask, k, wSel))), 'all', 'omitnan');
                denom = sqrt(den1 * den2);
                if isfinite(denom) && denom > 0
                    rSess(k) = num / denom;
                    pairCorrCondSess(p, c, wSel, k) = rSess(k);
                    wSess(k) = nUseSess;
                end
            end

            pairTrialCond(p, c, wSel) = sum(wSess, 'all', 'omitnan');
            validSess = isfinite(rSess) & isfinite(wSess) & (wSess > 0);
            if any(validSess) && pairTrialCond(p, c, wSel) >= P.minTrials
                zSess = atanh(max(min(rSess(validSess), 0.999999), -0.999999));
                switch string(P.sessionCombineMode)
                    case "equal"
                        pairCorrCond(p, c, wSel) = tanh(mean(zSess));
                    case "trial_weighted"
                        pairCorrCond(p, c, wSel) = tanh(sum(wSess(validSess) .* zSess) ./ sum(wSess(validSess)));
                    otherwise
                        error('Unsupported sessionCombineMode: %s', char(P.sessionCombineMode));
                end
            end
        end
    end
    if P.verbose && (mod(p, 2000) == 0 || p == nPairs)
        fprintf('  Pair summaries: %d / %d\n', p, nPairs);
    end
end

[meanByBin, semByBin, nPairsByBin, meanAboveMinRf, semAboveMinRf, nPairsAboveMinRf, ...
    meanLateMinusPre, semLateMinusPre, nPairsLateMinusPre] = ...
    local_summarize_pairs(pairCorrCond, Main.rfDistancePx, Main.condNames, Main.windowIdx, P);

exampleV1 = local_pick_example_site(driftV1, Main.SiteTableV1.globalSite);
exampleV4 = local_pick_example_site(driftV4, Main.SiteTableV4.globalSite);

OUT = struct();
OUT.P = P;
OUT.cacheParams = cacheParams;
OUT.siteSessionExclusions = site_session_exclusions(monkeySuffix);
OUT.monkeySuffix = monkeySuffix;
OUT.mainPath = Main.masterPath;
OUT.mainSummary = struct();
OUT.mainSummary.meanAboveMinRf = Main.meanAboveMinRf;
OUT.mainSummary.semAboveMinRf = Main.semAboveMinRf;
OUT.mainSummary.nPairsAboveMinRf = Main.nPairsAboveMinRf;
OUT.mainSummary.meanLateMinusPre = Main.meanLateMinusPre;
OUT.mainSummary.semLateMinusPre = Main.semLateMinusPre;
OUT.mainSummary.nPairsLateMinusPre = Main.nPairsLateMinusPre;
OUT.windowIdx = Main.windowIdx;
OUT.timeWindowsSelected = Main.timeWindowsSelected;
OUT.condNames = Main.condNames;
OUT.rfDistancePx = Main.rfDistancePx;
OUT.binEdges = Main.binEdges;
OUT.binCenters = Main.binCenters;
OUT.pairCorrCond = pairCorrCond;
OUT.pairTrialCond = pairTrialCond;
OUT.pairCorrCondSess = pairCorrCondSess;
OUT.pairTrialCondSess = pairTrialCondSess;
OUT.meanByBin = meanByBin;
OUT.semByBin = semByBin;
OUT.nPairsByBin = nPairsByBin;
OUT.meanAboveMinRf = meanAboveMinRf;
OUT.semAboveMinRf = semAboveMinRf;
OUT.nPairsAboveMinRf = nPairsAboveMinRf;
OUT.meanLateMinusPre = meanLateMinusPre;
OUT.semLateMinusPre = semLateMinusPre;
OUT.nPairsLateMinusPre = nPairsLateMinusPre;
OUT.loess = struct();
OUT.loess.frac = P.loessFrac;
OUT.loess.method = P.loessMethod;
OUT.loess.exampleV1 = exampleV1;
OUT.loess.exampleV4 = exampleV4;

if P.saveResult
    save(outPath, 'OUT', '-v7.3');
    if P.verbose
        fprintf('Saved V1-V4 LOESS drift-corrected analysis to %s\n', outPath);
    end
end

if P.plotFigure
    local_plot_summary(OUT);
    local_plot_difference_summary(OUT);
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

function [respDet, drift] = local_apply_loess(respRaw, dayIncl, sessions, loessFrac, loessMethod)
respDet = respRaw;
nSites = size(respRaw, 1);
nTrialsIncl = size(respRaw, 2);
nWin = size(respRaw, 3);
drift = nan(nSites, nTrialsIncl);
for s = 1:nSites
    for k = 1:numel(sessions)
        idxSess = find(dayIncl == sessions(k));
        if numel(idxSess) < 5
            continue;
        end
        xPre = reshape(respRaw(s, idxSess, 1), [], 1);
        good = isfinite(xPre);
        if nnz(good) < 5
            continue;
        end
        spanN = min(nnz(good), max(5, ceil(loessFrac * nnz(good))));
        if spanN >= nnz(good)
            fitVals = repmat(mean(xPre(good), 'omitnan'), nnz(good), 1);
        else
            fitVals = smoothdata(xPre(good), char(loessMethod), spanN);
        end
        d = nan(size(xPre));
        d(good) = fitVals;
        drift(s, idxSess) = d;
        for w = 1:nWin
            x = reshape(respRaw(s, idxSess, w), [], 1);
            goodW = isfinite(x) & isfinite(d);
            x(goodW) = x(goodW) - d(goodW);
            respDet(s, idxSess, w) = x;
        end
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

function [meanByBin, semByBin, nPairsByBin, meanAboveMinRf, semAboveMinRf, nPairsAboveMinRf, ...
    meanLateMinusPre, semLateMinusPre, nPairsLateMinusPre] = local_summarize_pairs(pairCorrCond, rfDistancePx, condNames, windowIdx, P)
nCond = numel(condNames);
nWinSel = numel(windowIdx);
maxDist = max(rfDistancePx(isfinite(rfDistancePx)));
binEdges = 0:P.binWidthPx:(ceil(maxDist / P.binWidthPx) * P.binWidthPx + P.binWidthPx);
if numel(binEdges) < 2
    binEdges = [0 P.binWidthPx];
end
nBins = numel(binEdges) - 1;

meanByBin = nan(nBins, nCond, nWinSel);
semByBin = nan(nBins, nCond, nWinSel);
nPairsByBin = zeros(nBins, nCond, nWinSel);
meanAboveMinRf = nan(nCond, nWinSel);
semAboveMinRf = nan(nCond, nWinSel);
nPairsAboveMinRf = zeros(nCond, nWinSel);
for wSel = 1:nWinSel
    for c = 1:nCond
        valid = isfinite(pairCorrCond(:, c, wSel)) & isfinite(rfDistancePx);
        binIdx = discretize(rfDistancePx(valid), binEdges);
        vals = pairCorrCond(valid, c, wSel);
        for b = 1:nBins
            xb = vals(binIdx == b);
            nPairsByBin(b, c, wSel) = numel(xb);
            if isempty(xb)
                continue;
            end
            meanByBin(b, c, wSel) = mean(xb, 'omitnan');
            if numel(xb) >= 2
                semByBin(b, c, wSel) = std(xb, 0, 'omitnan') / sqrt(numel(xb));
            else
                semByBin(b, c, wSel) = 0;
            end
        end

        validFar = valid & (rfDistancePx > P.summaryMinRfPx);
        valsFar = pairCorrCond(validFar, c, wSel);
        nPairsAboveMinRf(c, wSel) = numel(valsFar);
        if isempty(valsFar)
            continue;
        end
        meanAboveMinRf(c, wSel) = mean(valsFar, 'omitnan');
        if numel(valsFar) >= 2
            semAboveMinRf(c, wSel) = std(valsFar, 0, 'omitnan') / sqrt(numel(valsFar));
        else
            semAboveMinRf(c, wSel) = 0;
        end
    end
end

prePos = find(windowIdx == 1, 1, 'first');
latePos = find(windowIdx == 3, 1, 'first');
meanLateMinusPre = nan(nCond, 1);
semLateMinusPre = nan(nCond, 1);
nPairsLateMinusPre = zeros(nCond, 1);
if ~isempty(prePos) && ~isempty(latePos)
    for c = 1:nCond
        validDiff = isfinite(rfDistancePx) & (rfDistancePx > P.summaryMinRfPx) & ...
            isfinite(pairCorrCond(:, c, prePos)) & isfinite(pairCorrCond(:, c, latePos));
        delta = pairCorrCond(validDiff, c, latePos) - pairCorrCond(validDiff, c, prePos);
        nPairsLateMinusPre(c) = numel(delta);
        if isempty(delta)
            continue;
        end
        meanLateMinusPre(c) = mean(delta, 'omitnan');
        if numel(delta) >= 2
            semLateMinusPre(c) = std(delta, 0, 'omitnan') / sqrt(numel(delta));
        else
            semLateMinusPre(c) = 0;
        end
    end
end
end

function ex = local_pick_example_site(drift, siteGlobal)
meanDrift = mean(drift, 2, 'omitnan');
sdDrift = std(drift, 0, 2, 'omitnan');
[~, idx] = max(sdDrift);
ex = struct();
ex.globalSite = double(siteGlobal(idx));
ex.meanDrift = meanDrift(idx);
ex.sdDrift = sdDrift(idx);
ex.driftTrace = drift(idx, :);
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
S.cacheVersion = 3;
S.Monkey = double(Monkey);
S.cacheTag = string(P.cacheTag);
S.windowIdx = double(unique(P.windowIdx(:)'));
S.binWidthPx = double(P.binWidthPx);
S.minTrials = double(P.minTrials);
S.summaryMinRfPx = double(P.summaryMinRfPx);
S.sessionCombineMode = string(P.sessionCombineMode);
S.preWindowMs = double(P.preWindowMs);
S.earlyWindowMs = double(P.earlyWindowMs);
S.lateWindowMs = double(P.lateWindowMs);
S.loessFrac = double(P.loessFrac);
S.loessMethod = string(P.loessMethod);
S.minTrialsPerStimSess = double(P.minTrialsPerStimSess);
S.sessions = double(P.sessions(:)');
S.onlyCorrect = logical(P.onlyCorrect);
end

function tf = local_cache_matches(OUT, cacheParams)
tf = isstruct(OUT) && isfield(OUT, 'cacheParams') && isequaln(OUT.cacheParams, cacheParams);
end

function local_plot_difference_summary(OUT)
colors = [ ...
    0.15 0.45 0.85; ...
    0.20 0.60 0.30; ...
    0.10 0.70 0.55; ...
    0.85 0.40 0.15; ...
    0.15 0.15 0.15; ...
    0.45 0.45 0.45; ...
    0.58 0.58 0.58; ...
    0.72 0.72 0.72; ...
    0.84 0.84 0.84];

figure('Color', 'w', 'Name', sprintf('V1-V4 LOESS late-minus-pre (%s)', char(OUT.monkeySuffix)));
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

ax1 = nexttile;
hold(ax1, 'on');
x = 1:numel(OUT.condNames);
y = OUT.mainSummary.meanLateMinusPre(:);
e = OUT.mainSummary.semLateMinusPre(:);
b = bar(ax1, x, y, 'FaceColor', 'flat', 'EdgeColor', 'none');
b.CData = colors;
errorbar(ax1, x(isfinite(y)), y(isfinite(y)), e(isfinite(y)), 'k.', 'LineWidth', 1.0, 'HandleVisibility', 'off');
plot(ax1, x, zeros(size(x)), 'k:', 'HandleVisibility', 'off');
set(ax1, 'XTick', x, 'XTickLabel', cellstr(OUT.condNames));
xtickangle(ax1, 30);
ylabel(ax1, '\Delta noise correlation (late - pre)');
title(ax1, 'Current main');
grid(ax1, 'on');

ax2 = nexttile;
hold(ax2, 'on');
y = OUT.meanLateMinusPre(:);
e = OUT.semLateMinusPre(:);
b = bar(ax2, x, y, 'FaceColor', 'flat', 'EdgeColor', 'none');
b.CData = colors;
errorbar(ax2, x(isfinite(y)), y(isfinite(y)), e(isfinite(y)), 'k.', 'LineWidth', 1.0, 'HandleVisibility', 'off');
plot(ax2, x, zeros(size(x)), 'k:', 'HandleVisibility', 'off');
set(ax2, 'XTick', x, 'XTickLabel', cellstr(OUT.condNames));
xtickangle(ax2, 30);
ylabel(ax2, '\Delta noise correlation (late - pre)');
title(ax2, sprintf('LOESS corrected (span=%.2f)', OUT.loess.frac));
grid(ax2, 'on');
end

function local_plot_summary(OUT)
colors = [ ...
    0.15 0.45 0.85; ...
    0.85 0.40 0.15; ...
    0.20 0.60 0.30; ...
    0.15 0.15 0.15; ...
    0.45 0.45 0.45; ...
    0.72 0.72 0.72];
lineCondIdx = [1 2 3 4];

figure('Color', 'w', 'Name', sprintf('V1-V4 LOESS drift-corrected noise corr (%s)', char(OUT.monkeySuffix)));
tiledlayout(numel(OUT.windowIdx), 2, 'Padding', 'compact', 'TileSpacing', 'compact');

for wSel = 1:numel(OUT.windowIdx)
    tw = OUT.timeWindowsSelected(wSel, :);

    ax1 = nexttile;
    hold(ax1, 'on');
    for c = lineCondIdx
        x = OUT.binCenters;
        y = OUT.meanByBin(:, c, wSel);
        e = OUT.semByBin(:, c, wSel);
        good = isfinite(y);
        if any(good)
            patch(ax1, [x(good) fliplr(x(good))], ...
                [transpose(y(good)-e(good)) fliplr(transpose(y(good)+e(good)))], ...
                colors(c,:), 'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            plot(ax1, x(good), y(good), '-', 'Color', colors(c,:), 'LineWidth', 1.8, ...
                'DisplayName', sprintf('%s', OUT.condNames(c)));
        end
    end
    xlabel(ax1, 'RF distance (px)');
    ylabel(ax1, 'Mean noise correlation');
    title(ax1, sprintf('LOESS corrected | %d-%d ms', tw(1), tw(2)));
    grid(ax1, 'on');
    legend(ax1, 'Location', 'best', 'Box', 'off');

    ax2 = nexttile;
    hold(ax2, 'on');
    xBar = 1:numel(OUT.condNames);
    yBar = OUT.meanAboveMinRf(:, wSel);
    eBar = OUT.semAboveMinRf(:, wSel);
    b = bar(ax2, xBar, yBar, 'FaceColor', 'flat', 'EdgeColor', 'none');
    b.CData = colors;
    good = isfinite(yBar);
    errorbar(ax2, xBar(good), yBar(good), eBar(good), 'k.', 'LineWidth', 1.0, ...
        'HandleVisibility', 'off');
    set(ax2, 'XTick', xBar, 'XTickLabel', cellstr(OUT.condNames));
    xtickangle(ax2, 30);
    ylabel(ax2, 'Mean noise correlation');
    title(ax2, sprintf('RF distance > %d px', OUT.P.summaryMinRfPx));
    grid(ax2, 'on');
end
end
