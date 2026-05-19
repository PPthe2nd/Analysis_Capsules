%% Balanced time course: preferred vs non-preferred color for V4 (10 ms bins)
% Mirrors the current V1 gray-distance PSTH:
% - preferred/nonpreferred curves use balanced complementary pairs
% - the gray/background curve uses only gray RF stimuli that are at least
%   minDistThr pixels from the nearest colored curve

Monkey = 1; % 1 = Nilson, 2 = Figaro
useColorIndexFrom = "early"; % "early" or "late"
ciThr = 0.25; % abs(colorIndex) threshold for calling a site color tuned
useSNRgate = true;
snrThr = 0.7;
minDistThr = 30; % include gray only when the RF center is at least this far from the colored curves
AutoBuildColorTune = true; % compute/save the V4 ColorTune file if it is missing

COL_Y = "yellowArm";
COL_P = "purple";
COL_G = "gray";

cfg = config();

if Monkey == 1
    monkeySuffix = "N";
    tallFile = 'Tall_V4_lines_N.mat';
    colorTuneFile = 'ColorTune_balanced_V4_N.mat';
    resp3binFile = 'SNR_capsules_N_d12.mat';
    respFile = 'Resp_capsules_N_d12.mat';
elseif Monkey == 2
    monkeySuffix = "F";
    tallFile = 'Tall_V4_lines_F.mat';
    colorTuneFile = 'ColorTune_balanced_V4_F.mat';
    resp3binFile = 'SNR_capsules_F_d12.mat';
    respFile = 'Resp_capsules_F_d12.mat';
else
    error('PSTH_colorPref_V4:InvalidMonkey', ...
        'Monkey must be 1 (Nilson) or 2 (Figaro).');
end

tallPath = fullfile(cfg.matDir, tallFile);
colorTunePath = fullfile(cfg.matDir, colorTuneFile);
resp3binPath = fullfile(cfg.matDir, resp3binFile);
respPath = fullfile(cfg.matDir, respFile);

assert(exist(tallPath, 'file') == 2, ...
    'Missing %s. Run Line_Stimuli_V4.m first.', tallPath);
assert(exist(respPath, 'file') == 2, ...
    'Missing %s. Create the high-resolution response summary first.', respPath);

Sgeo = load(tallPath);
assert(isfield(Sgeo, 'Tall_V4') && isstruct(Sgeo.Tall_V4), ...
    '%s must contain struct Tall_V4.', tallFile);
assert(isfield(Sgeo, 'RFrange') && ~isempty(Sgeo.RFrange), ...
    '%s must contain RFrange.', tallFile);

Tall_V4 = Sgeo.Tall_V4;
RFrange = Sgeo.RFrange(:);
nV4 = numel(RFrange);
v4Sites = (1:nV4).';
hasSessionExclusions = ~isempty(site_session_exclusions(monkeySuffix));

R_full = load_capsules_struct_exclusion_aware(respPath, monkeySuffix, 'cfg', cfg);
R = R_full;
R.meanAct = R_full.meanAct(RFrange, :, :);
R.meanSqAct = R_full.meanSqAct(RFrange, :, :);
if ismatrix(R_full.nTrials) && size(R_full.nTrials,1) >= max(RFrange)
    R.nTrials = R_full.nTrials(RFrange, :);
else
    R.nTrials = R_full.nTrials;
end

if exist(colorTunePath, 'file') == 2 && ~hasSessionExclusions
    Sct = load(colorTunePath);
    assert(isfield(Sct, 'ColorTune') && isstruct(Sct.ColorTune), ...
        '%s must contain ColorTune.', colorTuneFile);
    ColorTune = Sct.ColorTune;
else
    assert(AutoBuildColorTune || hasSessionExclusions, ...
        'Missing %s. Run ColorTuning_Capsules_V4.m first or set AutoBuildColorTune = true.', colorTunePath);
    assert(exist(resp3binPath, 'file') == 2, ...
        'Missing %s. Need the 3-bin response summary to build ColorTune.', resp3binPath);

    if hasSessionExclusions
        fprintf(['Session exclusions are active for monkey %s; recomputing V4 ColorTune ' ...
                 'from exclusion-aware 3-bin responses.\n'], char(monkeySuffix));
    else
        fprintf('Color tuning file missing. Computing V4 ColorTune and saving to:\n%s\n', colorTunePath);
    end
    R3_full = load_capsules_struct_exclusion_aware(resp3binPath, monkeySuffix, 'cfg', cfg);
    R3 = R3_full;
    R3.meanAct = R3_full.meanAct(RFrange, :, :);
    R3.meanSqAct = R3_full.meanSqAct(RFrange, :, :);
    if ismatrix(R3_full.nTrials) && size(R3_full.nTrials,1) >= max(RFrange)
        R3.nTrials = R3_full.nTrials(RFrange, :);
    else
        R3.nTrials = R3_full.nTrials;
    end

    SNR3 = compute_snr_per_color_sites(R3, Tall_V4, v4Sites, 'Verbose', true);
    SNRmat3 = [SNR3.yellowEarly(v4Sites), SNR3.yellowLate(v4Sites), ...
               SNR3.purpleEarly(v4Sites), SNR3.purpleLate(v4Sites)];
    [bestSNR3, ~] = max(SNRmat3, [], 2, 'omitnan');

    keepSiteIdx3 = find(isfinite(bestSNR3) & (bestSNR3 > snrThr));
    assert(~isempty(keepSiteIdx3), ...
        'No V4 sites passed the SNR threshold %.2f during ColorTune construction.', snrThr);

    ColorTune = compute_color_tuning_balanced_sites(R3, Tall_V4, v4Sites, keepSiteIdx3, 'Verbose', true);
    ColorTune.thr = snrThr;
    ColorTune.bestSNR = bestSNR3;
    ColorTune.RFrange = RFrange;
    ColorTune.monkeySuffix = monkeySuffix;

    if ~hasSessionExclusions
        save(colorTunePath, 'ColorTune', '-v7.3');
        fprintf('Saved V4 color tuning results to %s\n', colorTunePath);
    end
end

[nCh, nStim2, nBins] = size(R.meanAct);
nStim = nStim2;
assert(nStim == 384, 'Expected 384 stimuli.');
assert(nCh == nV4, 'Localized response struct should have %d V4 rows, got %d.', nV4, nCh);
assert(nStim2 == nStim, 'R.meanAct stimulus dim mismatch.');
assert(size(R.timeWindows,1) == nBins, 'R.timeWindows rows must equal #bins.');

nTrialsRaw = R.nTrials;
if isvector(nTrialsRaw)
    nTrialsByStim = double(nTrialsRaw(:));
    perSiteTrials = false;
    assert(numel(nTrialsByStim) == nStim, 'R.nTrials vector must have %d elements.', nStim);
elseif ismatrix(nTrialsRaw) && size(nTrialsRaw,2) == nStim
    perSiteTrials = true;
    nTrialsByStim = [];
else
    error('R.nTrials must be a vector(%d) or matrix(nSites x %d).', nStim, nStim);
end

tCenters = mean(R.timeWindows, 2);

switch lower(useColorIndexFrom)
    case "early"
        ci = ColorTune.early.colorIndex(:);
    case "late"
        ci = ColorTune.late.colorIndex(:);
    otherwise
        error('useColorIndexFrom must be "early" or "late".');
end

bestSNR = ColorTune.bestSNR(:);
if useSNRgate
    keep = isfinite(ci) & (abs(ci) > ciThr) & isfinite(bestSNR) & (bestSNR > snrThr);
else
    keep = isfinite(ci) & (abs(ci) > ciThr);
end

keepSites = find(keep);
fprintf('Selected %d / %d V4 sites (abs(CI)>%.2f)%s\n', ...
    numel(keepSites), nV4, ciThr, ternary(useSNRgate, sprintf(' & bestSNR>%.2f', snrThr), ''));
assert(~isempty(keepSites), 'No V4 sites passed the preferred-color selection.');

prefIsYellow = ci > 0;

TallStimNums = arrayfun(@(x) x.stimNum, Tall_V4(:));
[sortedStimNums, order] = sort(TallStimNums(:));
assert(all(sortedStimNums(:).' == 1:nStim), 'Tall_V4.stimNum should cover 1..384.');
TallSorted = Tall_V4(order);

T0 = TallSorted(1).T;
assert(istable(T0), 'Tall_V4(stim).T must be a table.');
vn = string(T0.Properties.VariableNames);
ccIdx = find(vn == "center_color", 1);
if isempty(ccIdx)
    ccIdx = find(contains(lower(vn), "center") & contains(lower(vn), "color"), 1);
end
assert(~isempty(ccIdx), 'Could not find center_color column.');

distIdx = find(vn == "dist_to_nearest_color_px", 1);
assert(~isempty(distIdx), 'Could not find dist_to_nearest_color_px column.');

CC = strings(nV4, nStim);
Dist = nan(nV4, nStim);
for stim = 1:nStim
    Ti = TallSorted(stim).T;
    labs = string(Ti{:, ccIdx});
    CC(:, stim) = strtrim(labs);
    Dist(:, stim) = double(Ti{:, distIdx});
end

qualGray = (CC == COL_G) & isfinite(Dist) & (Dist >= minDistThr);
nQualGray = sum(qualGray, 2);
keepGraySite = keep & (nQualGray > 0);

fprintf(['Gray selection (for gray curve only): %d / %d color-selected V4 sites ' ...
         'have >=1 gray stimulus with dist>=%.1f px\n'], ...
    nnz(keepGraySite), numel(keepSites), minDistThr);

pairsA = zeros(nStim/2, 1);
pairsB = zeros(nStim/2, 1);
k = 0;
for a = 1:nStim
    pos = mod(a - 1, 8) + 1;
    if pos <= 4
        k = k + 1;
        pairsA(k) = a;
        pairsB(k) = a + 4;
    end
end
nPairs = k;

muPref = nan(numel(keepSites), nBins);
muNon  = nan(numel(keepSites), nBins);
muGray = nan(numel(keepSites), nBins);

for ii = 1:numel(keepSites)
    site = keepSites(ii);
    ch = v4Sites(site);

    sumY = zeros(nBins,1); NY = 0;
    sumP = zeros(nBins,1); NP = 0;
    sumG = zeros(nBins,1); NG = 0;
    doGray = keepGraySite(site);
    if perSiteTrials
        nTrSite = double(R.nTrials(ch, :)).';
    else
        nTrSite = nTrialsByStim;
    end

    for ip = 1:nPairs
        a = pairsA(ip);
        b = pairsB(ip);

        na = nTrSite(a);
        nb = nTrSite(b);

        ma = squeeze(R.meanAct(ch, a, :));
        mb = squeeze(R.meanAct(ch, b, :));

        ca = CC(site, a);
        cb = CC(site, b);

        if doGray
            if qualGray(site, a) && (na > 0)
                sumG = sumG + na * ma;
                NG = NG + na;
            end
            if qualGray(site, b) && (nb > 0)
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

    if prefIsYellow(site)
        muPref(ii, :) = muY;
        muNon(ii, :) = muP;
    else
        muPref(ii, :) = muP;
        muNon(ii, :) = muY;
    end
end

good = any(isfinite(muPref), 2) & any(isfinite(muNon), 2);
muPref = muPref(good, :);
muNon  = muNon(good, :);
muGray = muGray(good, :);
nPrefSites = size(muPref, 1);
nGraySites = sum(any(isfinite(muGray), 2));
fprintf('After requiring usable balanced Y/P trials: N=%d V4 sites\n', nPrefSites);

mPref = mean(muPref, 1, 'omitnan');
mNon  = mean(muNon,  1, 'omitnan');
mGray = mean(muGray, 1, 'omitnan');

semPref = std(muPref, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(muPref), 1));
semNon  = std(muNon,  0, 1, 'omitnan') ./ sqrt(sum(isfinite(muNon), 1));
semGray = std(muGray, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(muGray), 1));

figure; hold on;
hPref = plot(tCenters, mPref, 'LineWidth', 2);
hNon = plot(tCenters, mNon,  'LineWidth', 2);
hGray = plot(tCenters, mGray, 'LineWidth', 2);

fill([tCenters; flipud(tCenters)], [(mPref-semPref)'; flipud((mPref+semPref)')], ...
    'k', 'FaceAlpha', 0.12, 'EdgeColor', 'none');
fill([tCenters; flipud(tCenters)], [(mNon-semNon)'; flipud((mNon+semNon)')], ...
    'k', 'FaceAlpha', 0.06, 'EdgeColor', 'none');
fill([tCenters; flipud(tCenters)], [(mGray-semGray)'; flipud((mGray+semGray)')], ...
    'k', 'FaceAlpha', 0.04, 'EdgeColor', 'none');

xline(0, 'k-');
xlabel('Time from stimulus onset (ms)');
ylabel('Mean response (a.u.)');
title(sprintf(['V4 Balanced Pref vs Nonpref + Gray ' ...
    '(abs(CI)>%.2f, gray dist>=%.0f px, N pref/non=%d, N gray=%d)'], ...
    ciThr, minDistThr, nPrefSites, nGraySites));
legend([hPref, hNon, hGray], ...
    sprintf('Preferred (N=%d)', nPrefSites), ...
    sprintf('Nonpreferred (N=%d)', nPrefSites), ...
    sprintf('Gray (RF on background, N=%d)', nGraySites), ...
    'Location', 'best');
grid on;

function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end
