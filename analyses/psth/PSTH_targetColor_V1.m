function OUT = PSTH_targetColor_V1(Opts)
% PSTH_TARGETCOLOR_V1
% Balanced V1 time course for best vs worst target-capsule color,
% irrespective of RF location.
%
% This follows the standard V1 visual gate from the line-stimulus analysis:
% bestSNR > threshold.

if nargin < 1 || isempty(Opts)
    Opts = struct();
end

Opts = normalize_opts_local(Opts);
cfg = config();

tallPath = fullfile(cfg.matDir, 'Tall_V1_lines_N.mat');
respPath = fullfile(cfg.matDir, 'Resp_capsules_N_d12.mat');
resp3binPath = fullfile(cfg.matDir, 'SNR_capsules_N_d12.mat');

assert(exist(tallPath, 'file') == 2, ...
    'Missing %s. Run Line_Stimuli.m first.', tallPath);
assert(exist(respPath, 'file') == 2, ...
    'Missing %s. Create the high-resolution response summary first.', respPath);
assert(exist(resp3binPath, 'file') == 2, ...
    'Missing %s. Need the 3-bin response summary for gating/selection.', resp3binPath);

Sgeo = load(tallPath);
Sresp = load(respPath, 'R');
Sresp3 = load(resp3binPath, 'R');

assert(isfield(Sgeo, 'Tall_V1') && isstruct(Sgeo.Tall_V1), ...
    '%s must contain Tall_V1.', tallPath);
assert(isfield(Sresp, 'R') && isstruct(Sresp.R), ...
    '%s must contain struct R.', respPath);
assert(isfield(Sresp3, 'R') && isstruct(Sresp3.R), ...
    '%s must contain struct R.', resp3binPath);

Tall_V1 = Sgeo.Tall_V1;
if isfield(Sgeo, 'RTAB384')
    RTAB384 = Sgeo.RTAB384;
else
    rtabPath = fullfile(cfg.logsDir, 'RTAB384.mat');
    assert(exist(rtabPath, 'file') == 2, ...
        'Missing %s. Need RTAB384 to label target color.', rtabPath);
    Srtab = load(rtabPath);
    assert(isfield(Srtab, 'RTAB384'), ...
        '%s must contain RTAB384.', rtabPath);
    RTAB384 = Srtab.RTAB384;
end
v1Sites = (1:512).';
nV1 = numel(v1Sites);

R_resp = localize_response_rows_local(Sresp.R, v1Sites);
R3 = localize_response_rows_local(Sresp3.R, v1Sites);

[nCh, nStim, nBins] = size(R_resp.meanAct);
assert(nCh == nV1, 'Localized V1 response struct should have 512 rows.');
assert(nStim == 384, 'Expected 384 stimuli in %s.', respPath);
assert(size(R_resp.timeWindows, 1) == nBins, ...
    'R_resp.timeWindows rows must equal the number of bins.');
tCenters = mean(double(R_resp.timeWindows), 2);

[TallSorted, targetColorByStim] = build_target_color_labels_local(Tall_V1, RTAB384, nStim);
[pairsA, pairsB] = build_complementary_pairs_local(nStim);
assert(all(targetColorByStim(pairsA) ~= targetColorByStim(pairsB)), ...
    'Complementary pair target colors must swap in every pair.');

SNR = compute_snr_per_color_sites(R3, TallSorted, v1Sites, 'Verbose', false);
SNRmat = [SNR.yellowEarly(v1Sites), SNR.yellowLate(v1Sites), ...
          SNR.purpleEarly(v1Sites), SNR.purpleLate(v1Sites)];
bestSNR = max(SNRmat, [], 2, 'omitnan');
bestSNR = bestSNR(:);

keepBase = isfinite(bestSNR) & (bestSNR > Opts.VisualSNRthr);
baseLabel = sprintf('bestSNR > %.2f', Opts.VisualSNRthr);

TargetColor = compute_paired_target_color_sites(R3, targetColorByStim, pairsA, pairsB);
TC = TargetColor.(Opts.Window);
pPaired = TC.pairedP(:);
dPaired = TC.pairedWeightedDiff(:);
validTarget = keepBase & isfinite(pPaired) & isfinite(dPaired) & (dPaired ~= 0);
if Opts.RequirePairedSig
    keep = validTarget & (pPaired < Opts.PtargetThresh);
    selectionLabel = sprintf('%s, paired target-color p < %.3f (%s window)', ...
        baseLabel, Opts.PtargetThresh, Opts.Window);
else
    keep = validTarget;
    selectionLabel = sprintf('%s, valid paired target-color metric (%s window)', ...
        baseLabel, Opts.Window);
end

prefIsTargetYellow = dPaired > 0;

fprintf('V1 target-color PSTH (N, %s window)\n', Opts.Window);
fprintf('Base site pool: %s\n', baseLabel);
fprintf('Visual V1 sites in base pool: %d / %d\n', nnz(keepBase), nV1);
fprintf('Valid paired target-color metric in base pool: %d / %d\n', nnz(validTarget), nV1);
fprintf('Target-color pairedP < 0.05 in base pool: %d / %d\n', ...
    nnz(keepBase & isfinite(pPaired) & (pPaired < 0.05)), nV1);
fprintf('Using selection: %s\n', selectionLabel);
fprintf('Selected %d / %d V1 sites for PSTH\n', nnz(keep), nV1);

assert(any(keep), 'No V1 sites passed the requested target-color selection criteria.');

siteLocalSel = find(keep);
siteGlobalSel = siteLocalSel;
prefIsTargetYellowSel = prefIsTargetYellow(siteLocalSel);

nTrialsRaw = R_resp.nTrials;
if isvector(nTrialsRaw)
    nTrialsByStim = double(nTrialsRaw(:).');
    perSiteTrials = false;
    assert(numel(nTrialsByStim) == nStim, 'R_resp.nTrials vector must have %d elements.', nStim);
elseif ismatrix(nTrialsRaw) && (size(nTrialsRaw, 2) == nStim)
    perSiteTrials = true;
    nTrialsByStim = [];
    assert(size(nTrialsRaw, 1) >= nV1, ...
        'R_resp.nTrials has %d rows; need at least %d.', size(nTrialsRaw, 1), nV1);
else
    error('R_resp.nTrials must be a vector(384) or matrix(nSites x 384).');
end

if Opts.NormalizeResponses
    req = {'muSpont','muYellowEarly','muYellowLate','muPurpleEarly','muPurpleLate'};
    assert(all(isfield(SNR, req)), ...
        'SNR must contain muSpont/muYellowEarly/muYellowLate/muPurpleEarly/muPurpleLate.');
    bAll = double(SNR.muSpont(v1Sites));
    topMat = [double(SNR.muYellowEarly(v1Sites)), ...
              double(SNR.muYellowLate(v1Sites)), ...
              double(SNR.muPurpleEarly(v1Sites)), ...
              double(SNR.muPurpleLate(v1Sites))];
    scaleAll = max(topMat, [], 2) - bAll(:);
    scaleAll(~isfinite(scaleAll) | scaleAll <= 0) = NaN;
else
    bAll = zeros(nV1, 1);
    scaleAll = ones(nV1, 1);
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

    [muY, muP, NY, NP, nPairValid] = paired_target_curves_local( ...
        squeeze(double(R_resp.meanAct(iSite, :, :))), nTr, targetColorByStim, pairsA, pairsB);

    if Opts.NormalizeResponses
        b = bAll(iSite);
        sc = scaleAll(iSite);
        if isfinite(sc) && (sc > 0)
            muY = (muY - b) ./ sc;
            muP = (muP - b) ./ sc;
        else
            muY = nan(1, nBins);
            muP = nan(1, nBins);
            NY = 0; NP = 0; nPairValid = 0;
        end
    end

    if ~any(isfinite(muY)) || ~any(isfinite(muP))
        continue;
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

OUT = finalize_output_local(muBest, muWorst, muTargetYellow, muTargetPurple, ...
    nBestTrials, nWorstTrials, nPairsUsed, siteLocalSel, siteGlobalSel, ...
    prefIsTargetYellowSel, tCenters, selectionLabel, 'N', 'V1', Opts);
OUT.bestSNR = bestSNR;
OUT.TargetColor = TargetColor;
OUT.keepBase = keepBase;
OUT.validTarget = validTarget;
end

function Opts = normalize_opts_local(Opts)
defaults = struct();
defaults.Window = 'early';
defaults.PtargetThresh = 0.05;
defaults.RequirePairedSig = false;
defaults.VisualSNRthr = 0.7;
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

function [TallSorted, targetColorByStim] = build_target_color_labels_local(Tall, RTAB384, nStim)
stimNums = arrayfun(@(x) x.stimNum, Tall(:));
[stimNumsSorted, ord] = sort(stimNums(:));
assert(numel(stimNumsSorted) == nStim && all(stimNumsSorted(:).' == 1:nStim), ...
    'Tall.stimNum must cover 1..%d exactly.', nStim);
TallSorted = Tall(ord);

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

function [muY, muP, NY, NP, nPairValid] = paired_target_curves_local(respByStim, nTr, targetColorByStim, pairsA, pairsB)
[nStim, nBins] = size(respByStim);
if nStim ~= numel(nTr) && size(respByStim, 2) == numel(nTr)
    respByStim = respByStim.';
    [nStim, nBins] = size(respByStim);
end
assert(nStim == numel(nTr), 'Response/trial stimulus count mismatch.');

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

    ma = respByStim(a, :).';
    mb = respByStim(b, :).';
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

if (NY > 0) && (NP > 0)
    muY = (sumY / NY).';
    muP = (sumP / NP).';
else
    muY = nan(1, nBins);
    muP = nan(1, nBins);
end
end

function OUT = finalize_output_local(muBest, muWorst, muTargetYellow, muTargetPurple, ...
        nBestTrials, nWorstTrials, nPairsUsed, siteLocalSel, siteGlobalSel, ...
        prefIsTargetYellowSel, tCenters, selectionLabel, monkeySuffix, areaLabel, Opts)
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

fprintf('After requiring usable balanced target-color trials: N=%d %s sites\n', ...
    numel(siteLocalUsed), areaLabel);
if ~isempty(nPairsUsed)
    fprintf('Median balanced pair count per used site: %.1f\n', median(nPairsUsed));
    fprintf('Median balanced trial totals per used site: best=%g, worst=%g\n', ...
        median(nBestTrials), median(nWorstTrials));
end
assert(~isempty(siteLocalUsed), ...
    'No selected sites had usable balanced target-color trials.');

mBest = mean(muBest, 1, 'omitnan');
mWorst = mean(muWorst, 1, 'omitnan');
mTargetYellow = mean(muTargetYellow, 1, 'omitnan');
mTargetPurple = mean(muTargetPurple, 1, 'omitnan');

semBest = std(muBest, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(muBest), 1));
semWorst = std(muWorst, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(muWorst), 1));
semTargetYellow = std(muTargetYellow, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(muTargetYellow), 1));
semTargetPurple = std(muTargetPurple, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(muTargetPurple), 1));

figTitle = sprintf('%s target-color coding (%s) (%s, N=%d)', ...
    areaLabel, char(monkeySuffix), selectionLabel, numel(siteLocalUsed));

fig = [];
ax = [];
if Opts.PlotFigure
    fig = figure('Color', 'w', 'Name', figTitle, 'NumberTitle', 'off', ...
        'Tag', sprintf('%s_targetColor', areaLabel));
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
OUT.filters = Opts;
OUT.selectionLabel = selectionLabel;
OUT.monkeySuffix = monkeySuffix;
OUT.siteLocalSelected = siteLocalSel;
OUT.siteGlobalSelected = siteGlobalSel;
OUT.siteLocalUsed = siteLocalUsed;
OUT.siteGlobalUsed = siteGlobalUsed;
OUT.prefIsTargetYellowSelected = prefIsTargetYellowSel;
OUT.prefIsTargetYellowUsed = prefIsTargetYellowUsed;
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

function plot_sem_band_local(ax, t, mu, sem, color, alphaVal)
if ~any(isfinite(mu))
    return;
end
lo = mu(:) - sem(:);
hi = mu(:) + sem(:);
fill(ax, [t(:); flipud(t(:))], [lo; flipud(hi)], color, ...
    'FaceAlpha', alphaVal, 'EdgeColor', 'none');
end
