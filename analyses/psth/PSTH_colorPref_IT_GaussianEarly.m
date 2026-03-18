function OUT = PSTH_colorPref_IT_GaussianEarly(Monkey, Opts)
% PSTH_COLORPREF_IT_GAUSSIANEARLY
% Balanced IT time course for best vs worst color plus background using
% early Gaussian RF centers and the early precision-weighted paired color
% significance from GaussianOccupancy_Tuning_IT.

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
    gaussFile = 'GaussianOccupancy_Tuning_IT_N.mat';
    tallBaseFile = 'Tall_IT_lines_N.mat';
    respFile = 'Resp_capsules_N_d12.mat';
elseif Monkey == 2
    monkeySuffix = "F";
    gaussFile = 'GaussianOccupancy_Tuning_IT_F.mat';
    tallBaseFile = 'Tall_IT_lines_F.mat';
    respFile = 'Resp_capsules_F_d12.mat';
else
    error('PSTH_colorPref_IT_GaussianEarly:InvalidMonkey', ...
        'Monkey must be 1 (Nilson) or 2 (Figaro).');
end

gaussPath = fullfile(cfg.matDir, gaussFile);
tallBasePath = fullfile(cfg.matDir, tallBaseFile);
respPath = fullfile(cfg.matDir, respFile);

assert(exist(gaussPath, 'file') == 2, ...
    'Missing %s. Run GaussianOccupancy_Tuning_IT.m first.', gaussPath);
assert(exist(tallBasePath, 'file') == 2, ...
    'Missing %s. Run Line_Stimuli_IT.m first.', tallBasePath);
assert(exist(respPath, 'file') == 2, ...
    'Missing %s. Create the high-resolution IT response summary first.', respPath);

Sg = load(gaussPath, 'OUT');
Sgeo = load(tallBasePath, 'ALLCOORDS', 'RTAB384');
Sresp = load(respPath, 'R');

assert(isfield(Sg, 'OUT') && isstruct(Sg.OUT), ...
    '%s must contain struct OUT.', gaussPath);
assert(isfield(Sgeo, 'ALLCOORDS') && isfield(Sgeo, 'RTAB384'), ...
    '%s must contain ALLCOORDS and RTAB384.', tallBasePath);
assert(isfield(Sresp, 'R') && isstruct(Sresp.R), ...
    '%s must contain struct R.', respPath);

GOUT = Sg.OUT;
ALLCOORDS = Sgeo.ALLCOORDS;
RTAB384 = Sgeo.RTAB384;
R_full = Sresp.R;

assert(isfield(GOUT, 'FitSpatialEarly') && isfield(GOUT, 'PairColorEarly') && isfield(GOUT, 'RFrange'), ...
    '%s must contain FitSpatialEarly, PairColorEarly, and RFrange.', gaussFile);

FitEarly = GOUT.FitSpatialEarly;
PairEarly = GOUT.PairColorEarly;
RFrange = GOUT.RFrange(:);
nIT = numel(RFrange);

assert(numel(FitEarly) == nIT && numel(PairEarly) == nIT, ...
    'Gaussian early IT outputs must match RFrange length.');

centerX = [FitEarly.centerX].';
centerY = [FitEarly.centerY].';
sigmaPx = [FitEarly.sigmaPx].';
bestGaussianIdx = [FitEarly.bestGaussianIdx].';
spatialP = [FitEarly.pValueApprox].';
pairP = [PairEarly.pValue].';
weightedPairDiff = [PairEarly.weightedPairDiff].';
pairNPairs = [PairEarly.nPairs].';
pairStatus = string({PairEarly.status}).';

keep = isfinite(centerX) & isfinite(centerY) & ...
       isfinite(weightedPairDiff) & (weightedPairDiff ~= 0) & ...
       isfinite(pairP) & (pairP < Opts.PcolorThresh) & ...
       (pairStatus == "ok");

if Opts.RequireFiniteSigma
    keep = keep & isfinite(sigmaPx) & (sigmaPx > 0);
end
if Opts.RequireSpatialSig
    keep = keep & isfinite(spatialP) & (spatialP < Opts.PspatialThresh);
end
keep = keep & isfinite(bestGaussianIdx) & (bestGaussianIdx >= 1);

siteGlobalSel = RFrange(keep);
x_rf = centerX(keep);
y_rf = centerY(keep);
bestGaussianIdxSel = bestGaussianIdx(keep);
prefIsYellow = weightedPairDiff(keep) > 0;
pairPSel = pairP(keep);
weightedPairDiffSel = weightedPairDiff(keep);
pairNPairsSel = pairNPairs(keep);

[overlapYellow, overlapPurple, dominanceRatio] = load_overlap_cache_local(GOUT, cfg, monkeySuffix);

fprintf('IT Gaussian-early color PSTH (%s)\n', char(monkeySuffix));
fprintf('Selecting sites with early weighted paired p < %.3f\n', Opts.PcolorThresh);
if Opts.RequireSpatialSig
    fprintf('Also requiring early spatial p < %.3f\n', Opts.PspatialThresh);
end
fprintf('Best/worst conditions use Gaussian overlap dominance ratio >= %.2f\n', dominanceRatio);
fprintf('Selected %d / %d IT sites with finite Gaussian-early centers\n', ...
    numel(siteGlobalSel), nIT);

assert(~isempty(siteGlobalSel), ...
    'No IT sites passed the requested Gaussian-early color-selection criteria.');

fprintf('Rebuilding 384-stimulus IT geometry for the selected Gaussian-early centers\n');
Tall_IT_GaussianEarly = build_all_stim_tables(ALLCOORDS, RTAB384, x_rf, y_rf);

R_resp = localize_response_rows_local(R_full, siteGlobalSel);
[nSel, nStim, nBins] = size(R_resp.meanAct);
assert(nSel == numel(siteGlobalSel), 'Localized response rows do not match selected IT sites.');
assert(nStim == 384, 'Expected 384 stimuli in %s.', respFile);
assert(size(R_resp.timeWindows, 1) == nBins, ...
    'R_resp.timeWindows rows must equal the number of bins.');

tCenters = mean(double(R_resp.timeWindows), 2);
nTrialsRaw = R_resp.nTrials;
if isvector(nTrialsRaw)
    nTrialsByStim = double(nTrialsRaw(:).');
    perSiteTrials = false;
    assert(numel(nTrialsByStim) == nStim, 'R_resp.nTrials vector must have %d elements.', nStim);
elseif ismatrix(nTrialsRaw) && (size(nTrialsRaw, 2) == nStim)
    perSiteTrials = true;
    nTrialsByStim = [];
    assert(size(nTrialsRaw, 1) >= nSel, ...
        'R_resp.nTrials has %d rows; need at least %d.', size(nTrialsRaw, 1), nSel);
else
    error('R_resp.nTrials must be a vector(384) or matrix(nSites x 384).');
end

[TallSorted, CC, Dist] = build_color_label_mats_local(Tall_IT_GaussianEarly, nSel, nStim); %#ok<ASGLU>
assert(all(bestGaussianIdxSel <= size(overlapYellow, 2)), ...
    'Selected Gaussian index exceeds cached overlap library size.');

COL_G = "gray";

qualGray = (CC == COL_G) & isfinite(Dist) & (Dist >= Opts.MinDistThr);
nQualGray = sum(qualGray, 2);
keepGraySite = (nQualGray > 0);

fprintf(['Gray selection (for gray curve only): %d / %d selected IT sites ' ...
         'have >=1 gray stimulus with dist>=%.1f px\n'], ...
    nnz(keepGraySite), nSel, Opts.MinDistThr);

[pairsA, pairsB] = build_complementary_pairs_local(nStim);
nPairs = numel(pairsA);

muBest = nan(nSel, nBins);
muWorst = nan(nSel, nBins);
muGray = nan(nSel, nBins);
nDominantPairsHiRes = zeros(nSel, 1);
nBestTrials = zeros(nSel, 1);
nWorstTrials = zeros(nSel, 1);

for ii = 1:nSel
    if perSiteTrials
        nTr = double(nTrialsRaw(ii, :));
    else
        nTr = nTrialsByStim;
    end
    nTr(~isfinite(nTr) | nTr < 0) = 0;

    sumY = zeros(nBins, 1); NY = 0;
    sumP = zeros(nBins, 1); NP = 0;
    sumG = zeros(nBins, 1); NG = 0;
    doGray = keepGraySite(ii);
    bestIdx = bestGaussianIdxSel(ii);
    xY = double(overlapYellow(:, bestIdx));
    xP = double(overlapPurple(:, bestIdx));

    for ip = 1:nPairs
        a = pairsA(ip);
        b = pairsB(ip);

        na = nTr(a);
        nb = nTr(b);

        ma = squeeze(R_resp.meanAct(ii, a, :));
        mb = squeeze(R_resp.meanAct(ii, b, :));

        if doGray
            if qualGray(ii, a) && (na > 0)
                sumG = sumG + na * ma;
                NG = NG + na;
            end
            if qualGray(ii, b) && (nb > 0)
                sumG = sumG + nb * mb;
                NG = NG + nb;
            end
        end

        nEff = min(na, nb);
        if nEff <= 0
            continue;
        end

        isYellowA = is_color_dominant_local(xY(a), xP(a), dominanceRatio);
        isPurpleA = is_color_dominant_local(xP(a), xY(a), dominanceRatio);
        isYellowB = is_color_dominant_local(xY(b), xP(b), dominanceRatio);
        isPurpleB = is_color_dominant_local(xP(b), xY(b), dominanceRatio);

        if isYellowA && isPurpleB
            yStim = a;
            pStim = b;
        elseif isYellowB && isPurpleA
            yStim = b;
            pStim = a;
        else
            continue;
        end

        nDominantPairsHiRes(ii) = nDominantPairsHiRes(ii) + 1;
        sumY = sumY + nEff * squeeze(R_resp.meanAct(ii, yStim, :));
        sumP = sumP + nEff * squeeze(R_resp.meanAct(ii, pStim, :));
        NY = NY + nEff;
        NP = NP + nEff;
    end

    if NG > 0
        muGray(ii, :) = (sumG / NG).';
    end

    if (NY < 1) || (NP < 1)
        continue;
    end

    muY = (sumY / NY).';
    muP = (sumP / NP).';

    if prefIsYellow(ii)
        muBest(ii, :) = muY;
        muWorst(ii, :) = muP;
    else
        muBest(ii, :) = muP;
        muWorst(ii, :) = muY;
    end
    nBestTrials(ii) = NY;
    nWorstTrials(ii) = NP;
end

good = any(isfinite(muBest), 2) & any(isfinite(muWorst), 2);
siteGlobalUsed = siteGlobalSel(good);
prefIsYellowUsed = prefIsYellow(good);
pairPUsed = pairPSel(good);
weightedPairDiffUsed = weightedPairDiffSel(good);
pairNPairsUsed = pairNPairsSel(good);
nDominantPairsUsed = nDominantPairsHiRes(good);
nBestTrialsUsed = nBestTrials(good);
nWorstTrialsUsed = nWorstTrials(good);

muBest = muBest(good, :);
muWorst = muWorst(good, :);
muGray = muGray(good, :);

nBestSites = size(muBest, 1);
nGraySites = sum(any(isfinite(muGray), 2));

fprintf('After requiring usable balanced best/worst trials: N=%d IT sites\n', nBestSites);
fprintf('Gray curve uses subset via NaNs: N=%d IT sites\n', nGraySites);
if ~isempty(nDominantPairsUsed)
    fprintf('Median dominant complementary pairs per used site (high-res PSTH): %.1f\n', ...
        median(nDominantPairsUsed));
end

assert(nBestSites > 0, ...
    'No selected IT sites had usable balanced best/worst trials.');

mBest = mean(muBest, 1, 'omitnan');
mWorst = mean(muWorst, 1, 'omitnan');
mGray = mean(muGray, 1, 'omitnan');

semBest = std(muBest, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(muBest), 1));
semWorst = std(muWorst, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(muWorst), 1));
semGray = std(muGray, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(muGray), 1));

figTitle = sprintf(['IT Gaussian-early best vs worst + gray (%s) ' ...
    '(early weighted paired p<%.2f, gray dist>=%.0f px, N best/worst=%d, N gray=%d)'], ...
    char(monkeySuffix), Opts.PcolorThresh, Opts.MinDistThr, nBestSites, nGraySites);

fig = [];
ax = [];
if Opts.PlotFigure
    fig = figure('Color', 'w', 'Name', figTitle, 'NumberTitle', 'off', ...
        'Tag', 'IT_colorPref_GaussianEarly');
    ax = axes('Parent', fig);
    hold(ax, 'on');

    cBest = [0.80 0.15 0.15];
    cWorst = [0.15 0.35 0.85];
    cGray = [0.35 0.35 0.35];

    if Opts.PlotSem
        plot_sem_band_local(ax, tCenters, mBest, semBest, cBest, 0.12);
        plot_sem_band_local(ax, tCenters, mWorst, semWorst, cWorst, 0.10);
        plot_sem_band_local(ax, tCenters, mGray, semGray, cGray, 0.08);
    end

    hBest = plot(ax, tCenters, mBest, 'Color', cBest, 'LineWidth', 2);
    hWorst = plot(ax, tCenters, mWorst, 'Color', cWorst, 'LineWidth', 2);
    hGray = plot(ax, tCenters, mGray, 'Color', cGray, 'LineWidth', 2);

    xline(ax, 0, 'k-');
    xlabel(ax, 'Time from stimulus onset (ms)');
    ylabel(ax, 'Mean response (a.u.)');
    title(ax, figTitle);
    legend(ax, [hBest, hWorst, hGray], ...
        sprintf('Best color (N=%d)', nBestSites), ...
        sprintf('Worst color (N=%d)', nBestSites), ...
        sprintf('Gray in RF (N=%d)', nGraySites), ...
        'Location', 'best');
    grid(ax, 'on');
end

OUT = struct();
OUT.Monkey = Monkey;
OUT.monkeySuffix = monkeySuffix;
OUT.gaussPath = gaussPath;
OUT.tallBasePath = tallBasePath;
OUT.respPath = respPath;
OUT.filters = Opts;
OUT.keepMaskGaussian = keep;
OUT.siteGlobalSelected = siteGlobalSel;
OUT.siteGlobalUsed = siteGlobalUsed;
OUT.prefIsYellowSelected = prefIsYellow;
OUT.prefIsYellowUsed = prefIsYellowUsed;
OUT.pairedPSelected = pairPSel;
OUT.pairedPUsed = pairPUsed;
OUT.weightedPairDiffSelected = weightedPairDiffSel;
OUT.weightedPairDiffUsed = weightedPairDiffUsed;
OUT.nDominantPairsSelected = pairNPairsSel;
OUT.nDominantPairsUsed = pairNPairsUsed;
OUT.nDominantPairsHiResSelected = nDominantPairsHiRes;
OUT.nDominantPairsHiResUsed = nDominantPairsUsed;
OUT.nBestTrialsUsed = nBestTrialsUsed;
OUT.nWorstTrialsUsed = nWorstTrialsUsed;
OUT.keepGraySiteSelected = keepGraySite;
OUT.keepGraySiteUsed = keepGraySite(good);
OUT.x_rf = x_rf;
OUT.y_rf = y_rf;
OUT.bestGaussianIdxSelected = bestGaussianIdxSel;
OUT.Tall_IT_GaussianEarly = Tall_IT_GaussianEarly;
OUT.timeMs = tCenters;
OUT.muBestBySite = muBest;
OUT.muWorstBySite = muWorst;
OUT.muGrayBySite = muGray;
OUT.meanBest = mBest;
OUT.meanWorst = mWorst;
OUT.meanGray = mGray;
OUT.semBest = semBest;
OUT.semWorst = semWorst;
OUT.semGray = semGray;
OUT.figure = fig;
OUT.axes = ax;
end

function Opts = normalize_opts_local(Opts)
defaults = struct();
defaults.PcolorThresh = 0.05;
defaults.PspatialThresh = 0.05;
defaults.RequireSpatialSig = false;
defaults.RequireFiniteSigma = true;
defaults.MinDistThr = 30;
defaults.PlotFigure = true;
defaults.PlotSem = true;

fn = fieldnames(defaults);
for i = 1:numel(fn)
    if ~isfield(Opts, fn{i}) || isempty(Opts.(fn{i}))
        Opts.(fn{i}) = defaults.(fn{i});
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

function [overlapYellow, overlapPurple, dominanceRatio] = load_overlap_cache_local(GOUT, cfg, monkeySuffix)
if isfield(GOUT, 'cachePath') && ~isempty(GOUT.cachePath)
    cachePath = GOUT.cachePath;
else
    cachePath = fullfile(cfg.matDir, sprintf('GaussianOccupancy_Library_IT_%s.mat', char(monkeySuffix)));
end

assert(exist(cachePath, 'file') == 2, ...
    'Missing %s. Rerun GaussianOccupancy_Tuning_IT.m to rebuild the Gaussian overlap cache.', cachePath);

S = load(cachePath, 'overlapYellow', 'overlapPurple');
assert(isfield(S, 'overlapYellow') && isfield(S, 'overlapPurple'), ...
    '%s must contain overlapYellow and overlapPurple.', cachePath);

overlapYellow = S.overlapYellow;
overlapPurple = S.overlapPurple;

dominanceRatio = 2;
if isfield(GOUT, 'P') && isfield(GOUT.P, 'colorDominanceRatio') && isfinite(GOUT.P.colorDominanceRatio)
    dominanceRatio = double(GOUT.P.colorDominanceRatio);
end
end

function [TallSorted, CC, Dist] = build_color_label_mats_local(Tall, nSites, nStim)
stimNums = arrayfun(@(x) x.stimNum, Tall(:));
[stimNumsSorted, ord] = sort(stimNums(:));
assert(numel(stimNumsSorted) == nStim && all(stimNumsSorted(:).' == 1:nStim), ...
    'Tall.stimNum must cover 1..%d exactly.', nStim);
TallSorted = Tall(ord);

T0 = TallSorted(1).T;
assert(istable(T0), 'Tall(stim).T must be a table.');
vn = string(T0.Properties.VariableNames);

ccIdx = find(vn == "center_color", 1);
if isempty(ccIdx)
    ccIdx = find(contains(lower(vn), "center") & contains(lower(vn), "color"), 1);
end
assert(~isempty(ccIdx), 'Could not find center_color column.');

distIdx = find(vn == "dist_to_nearest_color_px", 1);
assert(~isempty(distIdx), 'Could not find dist_to_nearest_color_px column.');

CC = strings(nSites, nStim);
Dist = nan(nSites, nStim);
for stim = 1:nStim
    Ti = TallSorted(stim).T;
    labs = string(Ti{:, ccIdx});
    CC(:, stim) = strtrim(labs);
    Dist(:, stim) = double(Ti{:, distIdx});
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

function tf = is_color_dominant_local(mainOcc, otherOcc, ratioThr)
occFloor = 1e-9;
tf = isfinite(mainOcc) && isfinite(otherOcc) && (mainOcc > 0) && ...
    (mainOcc >= ratioThr * max(otherOcc, occFloor));
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
