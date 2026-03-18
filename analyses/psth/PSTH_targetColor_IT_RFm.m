function OUT = PSTH_targetColor_IT_RFm(Monkey, Opts)
% PSTH_TARGETCOLOR_IT_RFM
% Balanced IT time course for best vs worst target-capsule color,
% irrespective of where the RF falls.
%
% The key IT-specific choice is to work at the stimulus level:
%   - each stimulus is labeled by the color of the target capsule
%   - complementary stimulus pairs are used so geometry is matched while
%     the target/distractor colors swap
%   - site selection uses the paired precision-weighted target-color metric
%     because this behaved better than pooled metrics in IT
%
% By default, sites are first gated to be object-related and early
% responsive, then target-color significance is evaluated in the requested
% 3-bin window ('early' by default).

if nargin < 1 || isempty(Monkey)
    Monkey = 1; % 1 = Nilson, 2 = Figaro
end
if nargin < 2 || isempty(Opts)
    Opts = struct();
end

Opts = normalize_opts_local(Opts);
cfg = config();

if Monkey == 1
    monkeySuffix = "N";
    tallFile = 'Tall_IT_lines_N.mat';
    respFile = 'Resp_capsules_N_d12.mat';
    resp3binFile = 'SNR_capsules_N_d12.mat';
elseif Monkey == 2
    monkeySuffix = "F";
    tallFile = 'Tall_IT_lines_F.mat';
    respFile = 'Resp_capsules_F_d12.mat';
    resp3binFile = 'SNR_capsules_F_d12.mat';
else
    error('PSTH_targetColor_IT_RFm:InvalidMonkey', ...
        'Monkey must be 1 (Nilson) or 2 (Figaro).');
end

tallPath = fullfile(cfg.matDir, tallFile);
respPath = fullfile(cfg.matDir, respFile);
resp3binPath = fullfile(cfg.matDir, resp3binFile);

assert(exist(tallPath, 'file') == 2, ...
    'Missing %s. Run Line_Stimuli_IT.m first.', tallPath);
assert(exist(respPath, 'file') == 2, ...
    'Missing %s. Create the high-resolution IT response summary first.', respPath);
assert(exist(resp3binPath, 'file') == 2, ...
    'Missing %s. Need the 3-bin IT response summary for gating/selection.', resp3binPath);

Sgeo = load(tallPath);
Sresp = load(respPath, 'R');
Sresp3 = load(resp3binPath, 'R');

assert(isfield(Sgeo, 'Tall_IT') && isstruct(Sgeo.Tall_IT) && ...
       isfield(Sgeo, 'RFrange') && ~isempty(Sgeo.RFrange) && ...
       isfield(Sgeo, 'RTAB384'), ...
    '%s must contain Tall_IT, RFrange, and RTAB384.', tallPath);
assert(isfield(Sresp, 'R') && isstruct(Sresp.R), ...
    '%s must contain struct R.', respPath);
assert(isfield(Sresp3, 'R') && isstruct(Sresp3.R), ...
    '%s must contain struct R.', resp3binPath);

Tall_IT = Sgeo.Tall_IT;
RFrange = Sgeo.RFrange(:);
RTAB384 = Sgeo.RTAB384;
nIT = numel(RFrange);
siteRows = (1:nIT).';

R_full = Sresp.R;
R3_full = Sresp3.R;
R_resp = localize_response_rows_local(R_full, RFrange);
R3 = localize_response_rows_local(R3_full, RFrange);

[nSelResp, nStim, nBins] = size(R_resp.meanAct);
assert(nSelResp == nIT, 'Localized IT rows do not match RFrange.');
assert(nStim == 384, 'Expected 384 stimuli in %s.', respFile);
assert(size(R_resp.timeWindows, 1) == nBins, ...
    'R_resp.timeWindows rows must equal the number of bins.');
assert(size(R3.meanAct, 1) == nIT && size(R3.meanAct, 2) == nStim, ...
    'Localized 3-bin IT responses do not match geometry.');

tCenters = mean(double(R_resp.timeWindows), 2);

[TallSorted, targetColorByStim] = build_target_color_labels_local(Tall_IT, RTAB384, nStim);
[pairsA, pairsB] = build_complementary_pairs_local(nStim);
assert(numel(pairsA) == 192, 'Expected 192 complementary pairs.');

pairOpposite = targetColorByStim(pairsA) ~= targetColorByStim(pairsB);
assert(all(pairOpposite), ...
    'Complementary pair target colors must swap in every pair.');

SNR = compute_snr_per_color_sites(R3, TallSorted, siteRows, 'Verbose', false);
[hasObjectRF, nObjectStim, topKAbsQuartetEarly] = compute_it_responsive_gate_local( ...
    TallSorted, R3, SNR, Opts.MinObjectStim, Opts.TopKQuartets);

keepResponsive = hasObjectRF & isfinite(topKAbsQuartetEarly) & ...
    (topKAbsQuartetEarly > Opts.RespThr);
if Opts.RequireResponsive
    keepBase = keepResponsive;
    baseLabel = sprintf('object RF + top-%d early quartet > %.2f', ...
        Opts.TopKQuartets, Opts.RespThr);
else
    keepBase = true(nIT, 1);
    baseLabel = 'all IT sites';
end

TargetColor = compute_target_color_metrics_local(R3, targetColorByStim, pairsA, pairsB);
assert(isfield(TargetColor, Opts.Window), ...
    'TargetColor missing window "%s".', Opts.Window);
TC = TargetColor.(Opts.Window);

validTarget = keepBase & isfinite(TC.pairedP) & isfinite(TC.pairedWeightedDiff) & ...
    (TC.pairedWeightedDiff ~= 0);
if Opts.RequirePairedSig
    keep = validTarget & (TC.pairedP < Opts.PtargetThresh);
    selectionLabel = sprintf('%s, paired target-color p < %.3f (%s window)', ...
        baseLabel, Opts.PtargetThresh, Opts.Window);
else
    keep = validTarget;
    selectionLabel = sprintf('%s, valid paired target-color metric (%s window)', ...
        baseLabel, Opts.Window);
end

prefIsTargetYellow = TC.pairedWeightedDiff > 0;

fprintf('IT target-color PSTH (%s, %s window)\n', char(monkeySuffix), Opts.Window);
fprintf('Base site pool: %s\n', baseLabel);
fprintf('Object-related + early responsive IT sites: %d / %d\n', nnz(keepResponsive), nIT);
fprintf('Valid paired target-color metric in base pool: %d / %d\n', nnz(validTarget), nIT);
fprintf('Target-color pairedP < 0.05 in base pool: %d / %d\n', ...
    nnz(keepBase & isfinite(TC.pairedP) & (TC.pairedP < 0.05)), nIT);
fprintf('Using selection: %s\n', selectionLabel);
fprintf('Selected %d / %d IT sites for PSTH\n', nnz(keep), nIT);

assert(any(keep), ...
    'No IT sites passed the requested target-color selection criteria.');

siteLocalSel = find(keep);
siteGlobalSel = RFrange(siteLocalSel);
prefIsTargetYellowSel = prefIsTargetYellow(siteLocalSel);

nTrialsRaw = R_resp.nTrials;
if isvector(nTrialsRaw)
    nTrialsByStim = double(nTrialsRaw(:).');
    perSiteTrials = false;
    assert(numel(nTrialsByStim) == nStim, 'R_resp.nTrials vector must have %d elements.', nStim);
elseif ismatrix(nTrialsRaw) && (size(nTrialsRaw, 2) == nStim)
    perSiteTrials = true;
    nTrialsByStim = [];
    assert(size(nTrialsRaw, 1) >= nIT, ...
        'R_resp.nTrials has %d rows; need at least %d.', size(nTrialsRaw, 1), nIT);
else
    error('R_resp.nTrials must be a vector(384) or matrix(nSites x 384).');
end

if Opts.NormalizeResponses
    req = {'muSpont','muYellowEarly','muYellowLate','muPurpleEarly','muPurpleLate'};
    assert(all(isfield(SNR, req)), ...
        'SNR must contain muSpont/muYellowEarly/muYellowLate/muPurpleEarly/muPurpleLate.');
    bAll = double(SNR.muSpont(siteRows));
    topMat = [double(SNR.muYellowEarly(siteRows)), ...
              double(SNR.muYellowLate(siteRows)), ...
              double(SNR.muPurpleEarly(siteRows)), ...
              double(SNR.muPurpleLate(siteRows))];
    scaleAll = max(topMat, [], 2) - bAll(:);
    scaleAll(~isfinite(scaleAll) | scaleAll <= 0) = NaN;
else
    bAll = zeros(nIT, 1);
    scaleAll = ones(nIT, 1);
end

muBest = nan(numel(siteLocalSel), nBins);
muWorst = nan(numel(siteLocalSel), nBins);
muTargetYellow = nan(numel(siteLocalSel), nBins);
muTargetPurple = nan(numel(siteLocalSel), nBins);
nBestTrials = zeros(numel(siteLocalSel), 1);
nWorstTrials = zeros(numel(siteLocalSel), 1);
nPairsUsed = zeros(numel(siteLocalSel), 1);

for ii = 1:numel(siteLocalSel)
    iSite = siteLocalSel(ii);

    if perSiteTrials
        nTr = double(nTrialsRaw(iSite, :));
    else
        nTr = nTrialsByStim;
    end
    nTr(~isfinite(nTr) | nTr < 0) = 0;

    sumY = zeros(nBins, 1);
    sumP = zeros(nBins, 1);
    NY = 0;
    NP = 0;
    nPairValid = 0;

    for ip = 1:numel(pairsA)
        a = pairsA(ip);
        b = pairsB(ip);
        nEff = min(nTr(a), nTr(b));
        if ~(isfinite(nEff) && nEff > 0)
            continue;
        end

        ma = squeeze(double(R_resp.meanAct(iSite, a, :)));
        mb = squeeze(double(R_resp.meanAct(iSite, b, :)));
        if ~all(isfinite(ma)) || ~all(isfinite(mb))
            continue;
        end

        if targetColorByStim(a) == "yellowArm" && targetColorByStim(b) == "purple"
            sumY = sumY + nEff * ma;
            sumP = sumP + nEff * mb;
        elseif targetColorByStim(a) == "purple" && targetColorByStim(b) == "yellowArm"
            sumY = sumY + nEff * mb;
            sumP = sumP + nEff * ma;
        else
            continue;
        end

        NY = NY + nEff;
        NP = NP + nEff;
        nPairValid = nPairValid + 1;
    end

    if (NY < 1) || (NP < 1)
        continue;
    end

    muY = (sumY / NY).';
    muP = (sumP / NP).';

    if Opts.NormalizeResponses
        b = bAll(iSite);
        sc = scaleAll(iSite);
        if isfinite(sc) && (sc > 0)
            muY = (muY - b) ./ sc;
            muP = (muP - b) ./ sc;
        else
            muY = nan(1, nBins);
            muP = nan(1, nBins);
            NY = 0;
            NP = 0;
            nPairValid = 0;
        end
    end

    muTargetYellow(ii, :) = muY;
    muTargetPurple(ii, :) = muP;

    if prefIsTargetYellowSel(ii)
        muBest(ii, :) = muY;
        muWorst(ii, :) = muP;
        nBestTrials(ii) = NY;
        nWorstTrials(ii) = NP;
    else
        muBest(ii, :) = muP;
        muWorst(ii, :) = muY;
        nBestTrials(ii) = NP;
        nWorstTrials(ii) = NY;
    end
    nPairsUsed(ii) = nPairValid;
end

good = any(isfinite(muBest), 2) & any(isfinite(muWorst), 2);
siteLocalUsed = siteLocalSel(good);
siteGlobalUsed = siteGlobalSel(good);
prefIsTargetYellowUsed = prefIsTargetYellowSel(good);
muBest = muBest(good, :);
muWorst = muWorst(good, :);
muTargetYellow = muTargetYellow(good, :);
muTargetPurple = muTargetPurple(good, :);
nBestTrials = nBestTrials(good);
nWorstTrials = nWorstTrials(good);
nPairsUsed = nPairsUsed(good);

fprintf('After requiring usable balanced target-color trials: N=%d IT sites\n', numel(siteLocalUsed));
if ~isempty(nPairsUsed)
    fprintf('Median balanced pair count per used site: %.1f\n', median(nPairsUsed));
    fprintf('Median balanced trial totals per used site: best=%g, worst=%g\n', ...
        median(nBestTrials), median(nWorstTrials));
end

assert(~isempty(siteLocalUsed), ...
    'No selected IT sites had usable balanced target-color trials.');

mBest = mean(muBest, 1, 'omitnan');
mWorst = mean(muWorst, 1, 'omitnan');
mTargetYellow = mean(muTargetYellow, 1, 'omitnan');
mTargetPurple = mean(muTargetPurple, 1, 'omitnan');

semBest = std(muBest, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(muBest), 1));
semWorst = std(muWorst, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(muWorst), 1));
semTargetYellow = std(muTargetYellow, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(muTargetYellow), 1));
semTargetPurple = std(muTargetPurple, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(muTargetPurple), 1));

figTitle = sprintf(['IT target-color coding (%s) ' ...
    '(%s, N=%d)'], ...
    char(monkeySuffix), selectionLabel, numel(siteLocalUsed));

fig = [];
ax = [];
if Opts.PlotFigure
    fig = figure('Color', 'w', 'Name', figTitle, 'NumberTitle', 'off', ...
        'Tag', 'IT_targetColor_RFm');
    ax = axes('Parent', fig);
    hold(ax, 'on');

    cBest = [0.80 0.15 0.15];
    cWorst = [0.15 0.35 0.85];
    cTY = [0.85 0.60 0.10];
    cTP = [0.55 0.25 0.70];

    if Opts.PlotSem
        plot_sem_band_local(ax, tCenters, mBest, semBest, cBest, 0.12);
        plot_sem_band_local(ax, tCenters, mWorst, semWorst, cWorst, 0.10);
    end

    hBest = plot(ax, tCenters, mBest, 'Color', cBest, 'LineWidth', 2.2);
    hWorst = plot(ax, tCenters, mWorst, 'Color', cWorst, 'LineWidth', 2.2);
    hTY = plot(ax, tCenters, mTargetYellow, '--', 'Color', cTY, 'LineWidth', 1.6);
    hTP = plot(ax, tCenters, mTargetPurple, '--', 'Color', cTP, 'LineWidth', 1.6);

    xline(ax, 0, 'k-');
    xlabel(ax, 'Time from stimulus onset (ms)');
    ylabel(ax, 'Mean response (a.u.)');
    title(ax, figTitle);
    legend(ax, [hBest, hWorst, hTY, hTP], ...
        sprintf('Best target color (N=%d)', numel(siteLocalUsed)), ...
        sprintf('Worst target color (N=%d)', numel(siteLocalUsed)), ...
        'Unsorted: target yellow', ...
        'Unsorted: target purple', ...
        'Location', 'best');
    grid(ax, 'on');
end

OUT = struct();
OUT.Monkey = Monkey;
OUT.monkeySuffix = monkeySuffix;
OUT.tallPath = tallPath;
OUT.respPath = respPath;
OUT.resp3binPath = resp3binPath;
OUT.filters = Opts;
OUT.selectionLabel = selectionLabel;
OUT.targetColorByStim = targetColorByStim;
OUT.keepResponsive = keepResponsive;
OUT.keepBase = keepBase;
OUT.validTarget = validTarget;
OUT.keepMask = keep;
OUT.TargetColor = TargetColor;
OUT.siteLocalSelected = siteLocalSel;
OUT.siteGlobalSelected = siteGlobalSel;
OUT.siteLocalUsed = siteLocalUsed;
OUT.siteGlobalUsed = siteGlobalUsed;
OUT.prefIsTargetYellowSelected = prefIsTargetYellowSel;
OUT.prefIsTargetYellowUsed = prefIsTargetYellowUsed;
OUT.nObjectStim = nObjectStim;
OUT.topKAbsQuartetEarly = topKAbsQuartetEarly;
OUT.nBestTrialsUsed = nBestTrials;
OUT.nWorstTrialsUsed = nWorstTrials;
OUT.nPairsUsed = nPairsUsed;
OUT.timeMs = tCenters;
OUT.muBestBySite = muBest;
OUT.muWorstBySite = muWorst;
OUT.muTargetYellowBySite = muTargetYellow;
OUT.muTargetPurpleBySite = muTargetPurple;
OUT.meanBest = mBest;
OUT.meanWorst = mWorst;
OUT.meanTargetYellow = mTargetYellow;
OUT.meanTargetPurple = mTargetPurple;
OUT.semBest = semBest;
OUT.semWorst = semWorst;
OUT.semTargetYellow = semTargetYellow;
OUT.semTargetPurple = semTargetPurple;
OUT.figure = fig;
OUT.axes = ax;
end

function Opts = normalize_opts_local(Opts)
defaults = struct();
defaults.Window = 'early';
defaults.PtargetThresh = 0.05;
defaults.RequirePairedSig = true;
defaults.RequireResponsive = true;
defaults.MinObjectStim = 1;
defaults.RespThr = 0.7;
defaults.TopKQuartets = 5;
defaults.NormalizeResponses = true;
defaults.PlotFigure = true;
defaults.PlotSem = true;

fn = fieldnames(defaults);
for i = 1:numel(fn)
    if ~isfield(Opts, fn{i}) || isempty(Opts.(fn{i}))
        Opts.(fn{i}) = defaults.(fn{i});
    end
end
Opts.Window = char(string(Opts.Window));
end

function R_loc = localize_response_rows_local(R_full, siteGlobal)
R_loc = R_full;
R_loc.meanAct = R_full.meanAct(siteGlobal, :, :);
R_loc.meanSqAct = R_full.meanSqAct(siteGlobal, :, :);

if ismatrix(R_full.nTrials) && size(R_full.nTrials, 1) >= max(siteGlobal)
    R_loc.nTrials = R_full.nTrials(siteGlobal, :);
else
    R_loc.nTrials = R_full.nTrials;
end
end

function [TallSorted, targetColorByStim] = build_target_color_labels_local(Tall_IT, RTAB384, nStim)
stimNums = arrayfun(@(x) x.stimNum, Tall_IT(:));
[stimNumsSorted, ord] = sort(stimNums(:));
assert(numel(stimNumsSorted) == nStim && all(stimNumsSorted(:).' == 1:nStim), ...
    'Tall_IT.stimNum must cover 1..%d exactly.', nStim);
TallSorted = Tall_IT(ord);

assert(size(RTAB384, 1) >= nStim && size(RTAB384, 2) >= 8, ...
    'RTAB384 must have at least %d rows and 8 columns.', nStim);

targetColorByStim = strings(nStim, 1);
for stim = 1:nStim
    colIdx = double(RTAB384(stim, 8));
    if colIdx == 1
        targetColorByStim(stim) = "purple";
    elseif colIdx == 2
        targetColorByStim(stim) = "yellowArm";
    elseif mod(colIdx, 2) == 1
        targetColorByStim(stim) = "purple";
    else
        targetColorByStim(stim) = "yellowArm";
    end
end
end

function [pairsA, pairsB] = build_complementary_pairs_local(nStim)
pairsA = zeros(nStim / 2, 1);
pairsB = zeros(nStim / 2, 1);
k = 0;
for a = 1:nStim
    pos = mod(a - 1, 8) + 1;
    if pos <= 4
        k = k + 1;
        pairsA(k) = a;
        pairsB(k) = a + 4;
    end
end
pairsA = pairsA(1:k);
pairsB = pairsB(1:k);
end

function TargetColor = compute_target_color_metrics_local(DATA, targetColorByStim, pairsA, pairsB)
WIN_EARLY = 2;
WIN_LATE = 3;
useBesselCorrection = true;
varFloorFrac = 1e-3;

[nSitesTotal, nStim, nWin] = size(DATA.meanAct);
assert(all(size(DATA.meanSqAct) == size(DATA.meanAct)), ...
    'DATA.meanSqAct must match DATA.meanAct.');
assert(nStim == numel(targetColorByStim), ...
    'Stimulus count mismatch between responses and targetColorByStim.');
assert(nWin >= WIN_LATE, 'DATA must contain at least 3 time windows.');

nTrials = DATA.nTrials;
if isvector(nTrials)
    nTrials = double(nTrials(:));
    assert(numel(nTrials) == nStim, 'DATA.nTrials vector must have %d elements.', nStim);
    perSiteTrials = false;
elseif ismatrix(nTrials) && size(nTrials,1) == nSitesTotal && size(nTrials,2) == nStim
    nTrials = double(nTrials);
    perSiteTrials = true;
else
    error('DATA.nTrials must be a vector or [nSites x nStim] matrix.');
end

fields = ["pairedMuTargetYellow","pairedMuTargetPurple","pairedMeanDiff", ...
    "pairedWeightedDiff","pairedWeightSum","pairedEffN","pairedT","pairedP","pairedNPairs"];
for f = fields
    TargetColor.early.(f) = nan(nSitesTotal,1);
    TargetColor.late.(f)  = nan(nSitesTotal,1);
end

for iSite = 1:nSitesTotal
    mEarly = squeeze(DATA.meanAct(iSite, :, WIN_EARLY));   mEarly = mEarly(:);
    msqEarly = squeeze(DATA.meanSqAct(iSite, :, WIN_EARLY)); msqEarly = msqEarly(:);
    mLate = squeeze(DATA.meanAct(iSite, :, WIN_LATE));     mLate = mLate(:);
    msqLate = squeeze(DATA.meanSqAct(iSite, :, WIN_LATE)); msqLate = msqLate(:);

    pairDiffEarly = nan(numel(pairsA), 1);
    pairVarEarly = nan(numel(pairsA), 1);
    pairYEarly = nan(numel(pairsA), 1);
    pairPEarly = nan(numel(pairsA), 1);
    pairDiffLate = nan(numel(pairsA), 1);
    pairVarLate = nan(numel(pairsA), 1);
    pairYLate = nan(numel(pairsA), 1);
    pairPLate = nan(numel(pairsA), 1);
    nPairEarly = 0;
    nPairLate = 0;

    for k = 1:numel(pairsA)
        a = pairsA(k);
        b = pairsB(k);

        if perSiteTrials
            na = nTrials(iSite, a);
            nb = nTrials(iSite, b);
        else
            na = nTrials(a);
            nb = nTrials(b);
        end

        nEff = min(na, nb);
        if ~(isfinite(nEff) && nEff > 0)
            continue;
        end

        ca = targetColorByStim(a);
        cb = targetColorByStim(b);

        if ca == "yellowArm" && cb == "purple"
            [dEarly, vEarly, muYEarly, muPEarly] = build_paired_diff_local( ...
                mEarly(a), msqEarly(a), na, mEarly(b), msqEarly(b), nb, nEff, useBesselCorrection);
            [dLate, vLate, muYLate, muPLate] = build_paired_diff_local( ...
                mLate(a), msqLate(a), na, mLate(b), msqLate(b), nb, nEff, useBesselCorrection);
        elseif ca == "purple" && cb == "yellowArm"
            [dEarly, vEarly, muYEarly, muPEarly] = build_paired_diff_local( ...
                mEarly(b), msqEarly(b), nb, mEarly(a), msqEarly(a), na, nEff, useBesselCorrection);
            [dLate, vLate, muYLate, muPLate] = build_paired_diff_local( ...
                mLate(b), msqLate(b), nb, mLate(a), msqLate(a), na, nEff, useBesselCorrection);
        else
            dEarly = NaN; vEarly = NaN; muYEarly = NaN; muPEarly = NaN;
            dLate = NaN; vLate = NaN; muYLate = NaN; muPLate = NaN;
        end

        if isfinite(dEarly) && isfinite(vEarly)
            nPairEarly = nPairEarly + 1;
            pairDiffEarly(nPairEarly) = dEarly;
            pairVarEarly(nPairEarly) = vEarly;
            pairYEarly(nPairEarly) = muYEarly;
            pairPEarly(nPairEarly) = muPEarly;
        end
        if isfinite(dLate) && isfinite(vLate)
            nPairLate = nPairLate + 1;
            pairDiffLate(nPairLate) = dLate;
            pairVarLate(nPairLate) = vLate;
            pairYLate(nPairLate) = muYLate;
            pairPLate(nPairLate) = muPLate;
        end
    end

    [pairMuY, pairMuP, pairMeanDiff, pairWeightedDiff, pairWeightSum, pairEffN, pairT, pairP] = ...
        compute_paired_stats_local(pairDiffEarly(1:nPairEarly), pairVarEarly(1:nPairEarly), ...
        pairYEarly(1:nPairEarly), pairPEarly(1:nPairEarly), varFloorFrac);
    TargetColor.early.pairedMuTargetYellow(iSite) = pairMuY;
    TargetColor.early.pairedMuTargetPurple(iSite) = pairMuP;
    TargetColor.early.pairedMeanDiff(iSite) = pairMeanDiff;
    TargetColor.early.pairedWeightedDiff(iSite) = pairWeightedDiff;
    TargetColor.early.pairedWeightSum(iSite) = pairWeightSum;
    TargetColor.early.pairedEffN(iSite) = pairEffN;
    TargetColor.early.pairedT(iSite) = pairT;
    TargetColor.early.pairedP(iSite) = pairP;
    TargetColor.early.pairedNPairs(iSite) = nPairEarly;

    [pairMuY, pairMuP, pairMeanDiff, pairWeightedDiff, pairWeightSum, pairEffN, pairT, pairP] = ...
        compute_paired_stats_local(pairDiffLate(1:nPairLate), pairVarLate(1:nPairLate), ...
        pairYLate(1:nPairLate), pairPLate(1:nPairLate), varFloorFrac);
    TargetColor.late.pairedMuTargetYellow(iSite) = pairMuY;
    TargetColor.late.pairedMuTargetPurple(iSite) = pairMuP;
    TargetColor.late.pairedMeanDiff(iSite) = pairMeanDiff;
    TargetColor.late.pairedWeightedDiff(iSite) = pairWeightedDiff;
    TargetColor.late.pairedWeightSum(iSite) = pairWeightSum;
    TargetColor.late.pairedEffN(iSite) = pairEffN;
    TargetColor.late.pairedT(iSite) = pairT;
    TargetColor.late.pairedP(iSite) = pairP;
    TargetColor.late.pairedNPairs(iSite) = nPairLate;
end
end

function [pairDiff, pairVar, muY, muP] = build_paired_diff_local(muY, msqY, nOrigY, muP, msqP, nOrigP, nEff, useBessel)
pairDiff = NaN;
pairVar = NaN;
if ~(isfinite(muY) && isfinite(msqY) && isfinite(nOrigY) && nOrigY > 1 && ...
        isfinite(muP) && isfinite(msqP) && isfinite(nOrigP) && nOrigP > 1 && ...
        isfinite(nEff) && nEff > 0)
    return;
end

varMeanY = variance_of_balanced_mean_local(muY, msqY, nOrigY, nEff, useBessel);
varMeanP = variance_of_balanced_mean_local(muP, msqP, nOrigP, nEff, useBessel);
if ~(isfinite(varMeanY) && isfinite(varMeanP))
    return;
end

pairDiff = muY - muP;
pairVar = varMeanY + varMeanP;
end

function varMean = variance_of_balanced_mean_local(mu, msq, nOrig, nEff, useBessel)
varMean = NaN;
if ~(isfinite(mu) && isfinite(msq) && isfinite(nOrig) && nOrig > 1 && isfinite(nEff) && nEff > 0)
    return;
end
sampleVar = max(0, msq - mu^2);
if useBessel
    sampleVar = sampleVar * (nOrig / (nOrig - 1));
end
varMean = sampleVar / nEff;
end

function [muYw, muPw, meanDiff, weightedDiff, sumW, effN, tStat, pVal] = ...
        compute_paired_stats_local(pairDiff, pairVar, pairY, pairP, varFloorFrac)
muYw = NaN; muPw = NaN; meanDiff = NaN; weightedDiff = NaN;
sumW = NaN; effN = NaN; tStat = NaN; pVal = NaN;

valid = isfinite(pairDiff) & isfinite(pairVar) & (pairVar >= 0) & isfinite(pairY) & isfinite(pairP);
pairDiff = pairDiff(valid);
pairVar = pairVar(valid);
pairY = pairY(valid);
pairP = pairP(valid);
if numel(pairDiff) < 2
    return;
end

posVar = pairVar(pairVar > 0);
if isempty(posVar)
    varFloor = 1;
else
    varFloor = varFloorFrac * median(posVar);
end
varAdj = max(pairVar, varFloor);
weights = 1 ./ varAdj;
sumW = sum(weights);
if ~(isfinite(sumW) && sumW > 0)
    return;
end

weightedDiff = sum(weights .* pairDiff) / sumW;
muYw = sum(weights .* pairY) / sumW;
muPw = sum(weights .* pairP) / sumW;
meanDiff = mean(pairDiff);

effN = (sumW^2) / sum(weights.^2);
seWeighted = sqrt(1 / sumW);
if ~(isfinite(seWeighted) && seWeighted > 0)
    return;
end
tStat = weightedDiff / seWeighted;

if isfinite(effN) && (effN > 1)
    pVal = 2 * tcdf(-abs(tStat), effN - 1);
else
    pVal = NaN;
end
end

function [hasObjectRF, nObjectStim, topKAbsQuartetEarly] = compute_it_responsive_gate_local( ...
        Tall_IT, R3, SNR, MinObjectStim, TopKQuartets)
nIT = size(R3.meanAct, 1);
nStim = numel(Tall_IT);
siteRows = (1:nIT).';

nTargetStim = zeros(nIT, 1);
nDistrStim = zeros(nIT, 1);
for stim = 1:nStim
    T = Tall_IT(stim).T;
    assign = string(T.assignment(siteRows));
    nTargetStim = nTargetStim + (assign == "target");
    nDistrStim = nDistrStim + (assign == "distractor");
end
nObjectStim = nTargetStim + nDistrStim;
hasObjectRF = nObjectStim >= MinObjectStim;

if isvector(R3.nTrials)
    nTrialsAll = double(R3.nTrials(:));
    perSiteTrials = false;
else
    perSiteTrials = true;
end

muSpont = SNR.muSpont(siteRows);
sdSpont = SNR.sdSpont(siteRows);
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

muQuartetEarly = nan(nIT, nQuartets);
for iSite = 1:nIT
    if perSiteTrials
        nTrSite = double(R3.nTrials(iSite, :)).';
    else
        nTrSite = nTrialsAll;
    end

    rEarly = squeeze(R3.meanAct(iSite,:,2)).';
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

topKAbsQuartetEarly = nan(nIT, 1);
for iSite = 1:nIT
    vals = abs(signedQuartetEarly(iSite, :));
    vals = vals(isfinite(vals));
    if isempty(vals)
        continue;
    end
    vals = sort(vals, 'descend');
    kUse = min(TopKQuartets, numel(vals));
    topKAbsQuartetEarly(iSite) = mean(vals(1:kUse));
end
end

function plot_sem_band_local(ax, t, mu, sem, color, alphaVal)
if ~any(isfinite(mu))
    return;
end
lo = mu(:) - sem(:);
hi = mu(:) + sem(:);
fill(ax, [t(:); flipud(t(:))], [lo; flipud(hi)], color, ...
    'FaceAlpha', alphaVal, 'EdgeColor', 'none');
end
