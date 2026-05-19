function OUT = PSTH_colorPref_IT_RFm(Monkey, Opts)
% PSTH_COLORPREF_IT_RFM
% Balanced IT time course for best vs worst color plus background using
% the original RF.m-center Tall_IT geometry and RF-center color tuning.

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
    colorTuneFile = 'ColorTune_balanced_IT_N.mat';
    respFile = 'Resp_capsules_N_d12.mat';
elseif Monkey == 2
    monkeySuffix = "F";
    tallFile = 'Tall_IT_lines_F.mat';
    colorTuneFile = 'ColorTune_balanced_IT_F.mat';
    respFile = 'Resp_capsules_F_d12.mat';
else
    error('PSTH_colorPref_IT_RFm:InvalidMonkey', ...
        'Monkey must be 1 (Nilson) or 2 (Figaro).');
end

tallPath = fullfile(cfg.matDir, tallFile);
colorTunePath = fullfile(cfg.matDir, colorTuneFile);
respPath = fullfile(cfg.matDir, respFile);

assert(exist(tallPath, 'file') == 2, ...
    'Missing %s. Run Line_Stimuli_IT.m first.', tallPath);
assert(exist(colorTunePath, 'file') == 2, ...
    'Missing %s. Run ColorTuning_Capsules_IT.m first.', colorTunePath);
assert(exist(respPath, 'file') == 2, ...
    'Missing %s. Create the high-resolution IT response summary first.', respPath);

Sgeo = load(tallPath);
Sct = load(colorTunePath, 'ColorTune');

assert(isfield(Sgeo, 'Tall_IT') && isstruct(Sgeo.Tall_IT) && isfield(Sgeo, 'RFrange'), ...
    '%s must contain Tall_IT and RFrange.', tallPath);
assert(isfield(Sct, 'ColorTune') && isstruct(Sct.ColorTune), ...
    '%s must contain ColorTune.', colorTunePath);
Tall_IT = Sgeo.Tall_IT;
RFrange = Sgeo.RFrange(:);
ColorTune = Sct.ColorTune;
R_full = load_capsules_struct_exclusion_aware(respPath, monkeySuffix, 'cfg', cfg);
nIT = numel(RFrange);

assert(isfield(ColorTune, Opts.Window), ...
    'ColorTune is missing window "%s".', Opts.Window);
assert(isfield(ColorTune, 'keepSites') && ~isempty(ColorTune.keepSites), ...
    'ColorTune must contain keepSites.');
assert(isfield(ColorTune, 'RFrange') && numel(ColorTune.RFrange) == nIT, ...
    'ColorTune.RFrange must match the localized IT RFrange.');
assert(isequal(ColorTune.RFrange(:), RFrange), ...
    'ColorTune.RFrange does not match Tall_IT RFrange. Rebuild ColorTune_balanced_IT first.');

CT = ColorTune.(Opts.Window);
assert(isfield(CT, 'pairedP') && isfield(CT, 'pairedWeightedDiff') && ...
       isfield(CT, 'p') && isfield(CT, 'colorIndex'), ...
    'ColorTune.%s must contain pairedP, pairedWeightedDiff, p, and colorIndex.', Opts.Window);

keepSNR = false(nIT, 1);
keepSNR(ColorTune.keepSites(:)) = true;

pairedP = CT.pairedP(:);
pairedDiff = CT.pairedWeightedDiff(:);
pooledP = CT.p(:);
pooledIndex = CT.colorIndex(:);

validPaired = keepSNR & isfinite(pairedP) & isfinite(pairedDiff) & (pairedDiff ~= 0);
validPooled = keepSNR & isfinite(pooledP) & isfinite(pooledIndex) & (pooledIndex ~= 0);

fprintf('IT RF.m color PSTH (%s, %s window)\n', char(monkeySuffix), Opts.Window);
fprintf('SNR-kept IT sites from ColorTune: %d / %d\n', nnz(keepSNR), nIT);
fprintf('Valid paired-weighted color effect: %d / %d\n', nnz(validPaired), nIT);
fprintf('Valid pooled color effect: %d / %d\n', nnz(validPooled), nIT);
fprintf('Counts among SNR-kept IT sites:\n');
fprintf('  pairedP < 0.05: %d | pairedP < 0.10: %d\n', ...
    nnz(validPaired & (pairedP < 0.05)), nnz(validPaired & (pairedP < 0.10)));
fprintf('  pooled p < 0.05: %d | pooled p < 0.10: %d\n', ...
    nnz(validPooled & (pooledP < 0.05)), nnz(validPooled & (pooledP < 0.10)));

switch lower(Opts.SelectionMode)
    case 'paired'
        keep = validPaired & (pairedP < Opts.PcolorThresh);
        prefIsYellow = pairedDiff > 0;
        pSelectedAll = pairedP;
        effectSelectedAll = pairedDiff;
        selectionLabel = sprintf('pairedP < %.3f', Opts.PcolorThresh);
    case 'pooled'
        keep = validPooled & (pooledP < Opts.PcolorThresh);
        prefIsYellow = pooledIndex > 0;
        pSelectedAll = pooledP;
        effectSelectedAll = pooledIndex;
        selectionLabel = sprintf('pooled p < %.3f', Opts.PcolorThresh);
    otherwise
        error('SelectionMode must be ''paired'' or ''pooled''.');
end

siteLocalSel = find(keep);
siteGlobalSel = RFrange(siteLocalSel);
prefIsYellowSel = prefIsYellow(siteLocalSel);
pSelected = pSelectedAll(siteLocalSel);
effectSelected = effectSelectedAll(siteLocalSel);

fprintf('Using RF.m-center selection: %s\n', selectionLabel);
fprintf('Selected %d / %d IT sites for PSTH\n', numel(siteLocalSel), nIT);

assert(~isempty(siteLocalSel), ...
    'No IT sites passed the requested RF.m-center color-selection criteria.');

R_resp = localize_response_rows_local(R_full, siteGlobalSel);
[nSel, nStim, nBins] = size(R_resp.meanAct);
assert(nSel == numel(siteLocalSel), 'Localized response rows do not match selected IT sites.');
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
    assert(size(nTrialsRaw, 1) >= numel(siteGlobalSel), ...
        'R_resp.nTrials has %d rows; need at least %d.', size(nTrialsRaw, 1), numel(siteGlobalSel));
else
    error('R_resp.nTrials must be a vector(384) or matrix(nSites x 384).');
end

[TallSorted, CC, Dist] = build_color_label_mats_local(Tall_IT, siteLocalSel, nStim);

COL_Y = "yellowArm";
COL_P = "purple";
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

    for ip = 1:nPairs
        a = pairsA(ip);
        b = pairsB(ip);

        na = nTr(a);
        nb = nTr(b);

        ma = squeeze(R_resp.meanAct(ii, a, :));
        mb = squeeze(R_resp.meanAct(ii, b, :));

        ca = CC(ii, a);
        cb = CC(ii, b);

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

        if ca == COL_Y
            sumY = sumY + nEff * ma;
            NY = NY + nEff;
        elseif ca == COL_P
            sumP = sumP + nEff * ma;
            NP = NP + nEff;
        end

        if cb == COL_Y
            sumY = sumY + nEff * mb;
            NY = NY + nEff;
        elseif cb == COL_P
            sumP = sumP + nEff * mb;
            NP = NP + nEff;
        end
    end

    if NG > 0
        muGray(ii, :) = (sumG / NG).';
    end

    if (NY < 1) || (NP < 1)
        continue;
    end

    muY = (sumY / NY).';
    muP = (sumP / NP).';

    if prefIsYellowSel(ii)
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
siteLocalUsed = siteLocalSel(good);
siteGlobalUsed = siteGlobalSel(good);
prefIsYellowUsed = prefIsYellowSel(good);
pUsed = pSelected(good);
effectUsed = effectSelected(good);
nBestTrialsUsed = nBestTrials(good);
nWorstTrialsUsed = nWorstTrials(good);

muBest = muBest(good, :);
muWorst = muWorst(good, :);
muGray = muGray(good, :);

nBestSites = size(muBest, 1);
nGraySites = sum(any(isfinite(muGray), 2));

fprintf('After requiring usable balanced best/worst trials: N=%d IT sites\n', nBestSites);
fprintf('Gray curve uses subset via NaNs: N=%d IT sites\n', nGraySites);
if ~isempty(nBestTrialsUsed)
    fprintf('Median balanced trial totals per used site: best=%g, worst=%g\n', ...
        median(nBestTrialsUsed), median(nWorstTrialsUsed));
end

assert(nBestSites > 0, ...
    'No selected IT sites had usable balanced best/worst trials.');

mBest = mean(muBest, 1, 'omitnan');
mWorst = mean(muWorst, 1, 'omitnan');
mGray = mean(muGray, 1, 'omitnan');

semBest = std(muBest, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(muBest), 1));
semWorst = std(muWorst, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(muWorst), 1));
semGray = std(muGray, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(muGray), 1));

figTitle = sprintf(['IT RF.m-center best vs worst + gray (%s) ' ...
    '(%s, gray dist>=%.0f px, N best/worst=%d, N gray=%d)'], ...
    char(monkeySuffix), selectionLabel, Opts.MinDistThr, nBestSites, nGraySites);

fig = [];
ax = [];
if Opts.PlotFigure
    fig = figure('Color', 'w', 'Name', figTitle, 'NumberTitle', 'off', ...
        'Tag', 'IT_colorPref_RFm');
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
OUT.tallPath = tallPath;
OUT.colorTunePath = colorTunePath;
OUT.respPath = respPath;
OUT.filters = Opts;
OUT.selectionLabel = selectionLabel;
OUT.keepSNR = keepSNR;
OUT.validPaired = validPaired;
OUT.validPooled = validPooled;
OUT.keepMask = keep;
OUT.siteLocalSelected = siteLocalSel;
OUT.siteLocalUsed = siteLocalUsed;
OUT.siteGlobalSelected = siteGlobalSel;
OUT.siteGlobalUsed = siteGlobalUsed;
OUT.prefIsYellowSelected = prefIsYellowSel;
OUT.prefIsYellowUsed = prefIsYellowUsed;
OUT.pSelected = pSelected;
OUT.pUsed = pUsed;
OUT.effectSelected = effectSelected;
OUT.effectUsed = effectUsed;
OUT.keepGraySiteSelected = keepGraySite;
OUT.keepGraySiteUsed = keepGraySite(good);
OUT.nBestTrialsUsed = nBestTrialsUsed;
OUT.nWorstTrialsUsed = nWorstTrialsUsed;
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
defaults.Window = 'early';
defaults.SelectionMode = 'paired';
defaults.PcolorThresh = 0.05;
defaults.MinDistThr = 30;
defaults.PlotFigure = true;
defaults.PlotSem = true;

fn = fieldnames(defaults);
for i = 1:numel(fn)
    if ~isfield(Opts, fn{i}) || isempty(Opts.(fn{i}))
        Opts.(fn{i}) = defaults.(fn{i});
    end
end
Opts.Window = char(string(Opts.Window));
Opts.SelectionMode = char(string(Opts.SelectionMode));
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

function [TallSorted, CC, Dist] = build_color_label_mats_local(Tall, siteRows, nStim)
stimNums = arrayfun(@(x) x.stimNum, Tall(:));
[stimNumsSorted, ord] = sort(stimNums(:));
assert(numel(stimNumsSorted) == nStim && all(stimNumsSorted(:).' == 1:nStim), ...
    'Tall.stimNum must cover 1..%d exactly.', nStim);
TallSorted = Tall(ord);

T0 = TallSorted(1).T;
assert(istable(T0), 'Tall(stim).T must be a table.');
assert(height(T0) >= max(siteRows), ...
    'Tall(stim).T has %d rows; need at least %d.', height(T0), max(siteRows));
vn = string(T0.Properties.VariableNames);

ccIdx = find(vn == "center_color", 1);
if isempty(ccIdx)
    ccIdx = find(contains(lower(vn), "center") & contains(lower(vn), "color"), 1);
end
assert(~isempty(ccIdx), 'Could not find center_color column.');

distIdx = find(vn == "dist_to_nearest_color_px", 1);
assert(~isempty(distIdx), 'Could not find dist_to_nearest_color_px column.');

CC = strings(numel(siteRows), nStim);
Dist = nan(numel(siteRows), nStim);
for stim = 1:nStim
    Ti = TallSorted(stim).T;
    CC(:, stim) = strtrim(string(Ti{siteRows, ccIdx}));
    Dist(:, stim) = double(Ti{siteRows, distIdx});
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

function plot_sem_band_local(ax, t, mu, sem, color, alphaVal)
if ~any(isfinite(mu))
    return;
end
lo = mu(:) - sem(:);
hi = mu(:) + sem(:);
fill(ax, [t(:); flipud(t(:))], [lo; flipud(hi)], color, ...
    'FaceAlpha', alphaVal, 'EdgeColor', 'none');
end
