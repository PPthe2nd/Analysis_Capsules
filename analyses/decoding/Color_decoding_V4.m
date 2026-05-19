%% Weighted pair-wise color decoding over time for V4 (no fitting, complementary balanced)
% Mirrors the current V1 decoder:
% - each complementary pair is one sample
% - site weights are abs(dprime) from ColorTune.(refWin)
% - gray/unknown RF assignments are ignored for each pair
% - the final population trace is normalized by the smoothed peak of the
%   weighted visual response trace across all object-in-RF conditions

Monkey = 1; % 1 = Nilson, 2 = Figaro
refWin = "early"; % "early" or "late" -> which ColorTune window sets preference + |d'| weights
capWeightPct = 95; % cap |d'| weights at this percentile to limit domination by a few sites
ColorTuneSNRthr = 0.7; % SNR threshold used if the V4 ColorTune file has to be built on the fly
visualPeakWindowMs = [40 250]; % search for the visual-response peak only in this post-stimulus window
visualPeakSmoothMs = 30; % running-mean width for the visual-response trace before peak finding
useSNRsdSpontIfPresent = true; % use spontaneous SD from the 3-bin SNR summary when available
computeSpontSDIfMissing = true; % if no usable SNR SD is available, compute it from pre-0 high-res bins
normalizeBySpontSD = false; % optional legacy normalization: divide both decoder and visual-reference traces by each site's spontaneous SD
AutoBuildColorTune = true; % compute/save ColorTune on the fly if the V4 file is missing

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
    error('Color_decoding_V4:InvalidMonkey', ...
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
assert(exist(resp3binPath, 'file') == 2, ...
    'Missing %s. Need the 3-bin response summary for V4 decoding.', resp3binPath);

Sgeo = load(tallPath);
assert(isfield(Sgeo, 'Tall_V4') && isstruct(Sgeo.Tall_V4), ...
    '%s must contain struct Tall_V4.', tallFile);
assert(isfield(Sgeo, 'RFrange') && ~isempty(Sgeo.RFrange), ...
    '%s must contain RFrange.', tallFile);

Tall_V4 = Sgeo.Tall_V4;
RFrange = Sgeo.RFrange(:);
nV4 = numel(RFrange);
siteRows = (1:nV4).';
hasSessionExclusions = ~isempty(site_session_exclusions(monkeySuffix));

R_full = load_capsules_struct_exclusion_aware(respPath, monkeySuffix, 'cfg', cfg);
Rdec = R_full;
Rdec.meanAct = R_full.meanAct(RFrange, :, :);
Rdec.meanSqAct = R_full.meanSqAct(RFrange, :, :);
if ismatrix(R_full.nTrials) && size(R_full.nTrials,1) >= max(RFrange)
    Rdec.nTrials = R_full.nTrials(RFrange, :);
else
    Rdec.nTrials = R_full.nTrials;
end

R3_full = load_capsules_struct_exclusion_aware(resp3binPath, monkeySuffix, 'cfg', cfg);
R3dec = R3_full;
R3dec.meanAct = R3_full.meanAct(RFrange, :, :);
R3dec.meanSqAct = R3_full.meanSqAct(RFrange, :, :);
if ismatrix(R3_full.nTrials) && size(R3_full.nTrials,1) >= max(RFrange)
    R3dec.nTrials = R3_full.nTrials(RFrange, :);
else
    R3dec.nTrials = R3_full.nTrials;
end

SNR3 = compute_snr_per_color_sites(R3dec, Tall_V4, siteRows, 'Verbose', false);

if exist(colorTunePath, 'file') == 2 && ~hasSessionExclusions
    Sct = load(colorTunePath);
    assert(isfield(Sct, 'ColorTune') && isstruct(Sct.ColorTune), ...
        '%s must contain ColorTune.', colorTuneFile);
    ColorTune = Sct.ColorTune;
else
    assert(AutoBuildColorTune || hasSessionExclusions, ...
        'Missing %s. Run ColorTuning_Capsules_V4.m first or set AutoBuildColorTune = true.', colorTunePath);

    if hasSessionExclusions
        fprintf(['Session exclusions are active for monkey %s; recomputing V4 ColorTune ' ...
                 'from exclusion-aware 3-bin responses.\n'], char(monkeySuffix));
    else
        fprintf('Color tuning file missing. Computing V4 ColorTune and saving to:\n%s\n', colorTunePath);
    end

    SNRmat3 = [SNR3.yellowEarly(siteRows), SNR3.yellowLate(siteRows), ...
               SNR3.purpleEarly(siteRows), SNR3.purpleLate(siteRows)];
    [bestSNR3, ~] = max(SNRmat3, [], 2, 'omitnan');

    keepSiteIdx3 = find(isfinite(bestSNR3) & (bestSNR3 > ColorTuneSNRthr));
    assert(~isempty(keepSiteIdx3), ...
        'No V4 sites passed the SNR threshold %.2f during ColorTune construction.', ColorTuneSNRthr);

    ColorTune = compute_color_tuning_balanced_sites(R3dec, Tall_V4, siteRows, keepSiteIdx3, 'Verbose', true);
    ColorTune.thr = ColorTuneSNRthr;
    ColorTune.bestSNR = bestSNR3;
    ColorTune.RFrange = RFrange;
    ColorTune.monkeySuffix = monkeySuffix;

    if ~hasSessionExclusions
        save(colorTunePath, 'ColorTune', '-v7.3');
        fprintf('Saved V4 color tuning results to %s\n', colorTunePath);
    end
end

[nCh, nStim2, nBins] = size(Rdec.meanAct);
nStim = nStim2;
assert(nStim == 384, 'Expected 384 stimuli.');
assert(nCh == nV4, 'Localized response struct should have %d V4 rows, got %d.', nV4, nCh);
assert(nStim2 == nStim, 'R.meanAct stimulus dim mismatch.');
assert(all(size(Rdec.meanSqAct) == size(Rdec.meanAct)), 'R.meanSqAct must match R.meanAct size.');
assert(size(Rdec.timeWindows,1) == nBins, 'R.timeWindows rows must equal #bins.');

nTrialsRaw = Rdec.nTrials;
if isvector(nTrialsRaw)
    nTrialsByStim = double(nTrialsRaw(:));
    perSiteTrials = false;
    assert(numel(nTrialsByStim) == nStim, 'Rdec.nTrials vector must have %d elements.', nStim);
elseif ismatrix(nTrialsRaw) && size(nTrialsRaw,2) == nStim
    perSiteTrials = true;
    nTrialsByStim = [];
else
    error('Rdec.nTrials must be a vector(%d) or matrix(nSites x %d).', nStim, nStim);
end

tCenters = mean(Rdec.timeWindows, 2);

switch lower(refWin)
    case "early"
        ci = ColorTune.early.colorIndex(:);
        dp = ColorTune.early.dprime(:);
    case "late"
        ci = ColorTune.late.colorIndex(:);
        dp = ColorTune.late.dprime(:);
    otherwise
        error('refWin must be "early" or "late".');
end

prefIsYellow = (ci > 0);
prefIsYellow(~isfinite(ci)) = false;

w = abs(dp);
w(~isfinite(w)) = 0;
wPos = w(w > 0);
if ~isempty(wPos)
    wCap = prctile(wPos, capWeightPct);
    w95 = prctile(wPos, 95);
    w = min(w, wCap);
else
    wCap = NaN;
    w95 = NaN;
end

fprintf('V4 weight summary |d''|: median=%.3f, 95%%=%.3f, max(after cap)=%.3f\n', ...
    median(w(w > 0), 'omitnan'), w95, max(w));
useSites = find(w > 0);
fprintf('Using %d / %d V4 sites with nonzero weight\n', numel(useSites), nV4);
assert(~isempty(useSites), 'No V4 sites have nonzero decoding weight.');
useWeights = w(useSites);

TallStimNums = arrayfun(@(x) x.stimNum, Tall_V4(:));
[sortedStimNums, order] = sort(TallStimNums(:));
assert(all(sortedStimNums(:).' == 1:nStim), 'Tall_V4.stimNum should cover 1..384.');
TallSorted = Tall_V4(order);

T0 = TallSorted(1).T;
assert(istable(T0), 'Tall_V4(stim).T must be a table.');
vn = string(T0.Properties.VariableNames);
colIdx = find(vn == "center_color", 1);
if isempty(colIdx)
    colIdx = find(contains(lower(vn), "center") & contains(lower(vn), "color"), 1);
end
assert(~isempty(colIdx), 'Could not find center_color column in Tall_V4(stim).T.');

CC = strings(nV4, nStim);
for stim = 1:nStim
    CC(:,stim) = strtrim(string(TallSorted(stim).T{:, colIdx}));
end

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
assert(nPairs == 192, 'Expected 192 complementary pairs, got %d.', nPairs);
fprintf('Using %d complementary pairs\n', nPairs);

sdSpont = ones(nV4, 1);
if normalizeBySpontSD
    if useSNRsdSpontIfPresent && isfield(SNR3, 'sdSpont') && numel(SNR3.sdSpont) >= nV4
        sdSpont = SNR3.sdSpont(siteRows);
        sdSpont(~isfinite(sdSpont) | sdSpont <= 0) = 1;
        fprintf('Using sdSpont from 3-bin V4 SNR summary\n');
    elseif computeSpontSDIfMissing
        isSpontBin = (Rdec.timeWindows(:,2) <= 0);
        assert(any(isSpontBin), 'No spontaneous bins found (timeWindows end<=0).');

        bIdx = find(isSpontBin);
        fprintf('Computing sdSpont from %d pre-0 bins using high-resolution V4 responses...\n', numel(bIdx));
        for site = 1:nV4
            if perSiteTrials
                nTrSite = double(Rdec.nTrials(site,:)).';
            else
                nTrSite = nTrialsByStim;
            end
            wStim = nTrSite;
            wStim = wStim / sum(wStim);
            varBins = nan(numel(bIdx),1);
            for ib = 1:numel(bIdx)
                tb = bIdx(ib);

                m = squeeze(Rdec.meanAct(site,:,tb)).';
                m2 = squeeze(Rdec.meanSqAct(site,:,tb)).';

                good = isfinite(m) & isfinite(m2) & (nTrSite > 0);
                if ~any(good)
                    continue;
                end

                ww = wStim(good);
                ww = ww / sum(ww);
                mu = sum(ww .* m(good));
                ex2 = sum(ww .* m2(good));
                varBins(ib) = max(0, ex2 - mu^2);
            end

            sd = sqrt(mean(varBins, 'omitnan'));
            if isfinite(sd) && sd > 0
                sdSpont(site) = sd;
            else
                sdSpont(site) = 1;
            end
        end
        fprintf('Done computing sdSpont\n');
    else
        fprintf('normalizeBySpontSD=true but no usable sdSpont available; using sdSpont=1\n');
    end
end

scorePerPair = nan(nPairs, nBins);
nUsedSitesPerPair = nan(nPairs, 1);

for ip = 1:nPairs
    a = pairsA(ip);
    b = pairsB(ip);

    Ra = squeeze(Rdec.meanAct(siteRows, a, :));
    Rb = squeeze(Rdec.meanAct(siteRows, b, :));

    sumScore = zeros(1, nBins);
    sumW = zeros(1, nBins);
    nUsed = 0;

    for site = useSites(:).'
        ca = CC(site, a);
        cb = CC(site, b);

        if ca == COL_G || cb == COL_G || ca == "" || cb == ""
            continue;
        end

        if ca == COL_Y && cb == COL_P
            rY = Ra(site, :);
            rP = Rb(site, :);
        elseif ca == COL_P && cb == COL_Y
            rY = Rb(site, :);
            rP = Ra(site, :);
        else
            continue;
        end

        wi = w(site);
        if wi <= 0
            continue;
        end

        if prefIsYellow(site)
            d = (rY - rP);
        else
            d = (rP - rY);
        end

        if normalizeBySpontSD
            d = d ./ sdSpont(site);
        end

        sumScore = sumScore + wi * d;
        sumW = sumW + wi;
        nUsed = nUsed + 1;
    end

    nUsedSitesPerPair(ip) = nUsed;

    ok = sumW > 0;
    tmp = nan(1, nBins);
    tmp(ok) = sumScore(ok) ./ sumW(ok);
    scorePerPair(ip, :) = tmp;
end

fprintf('Median #V4 sites contributing per pair: %.1f\n', median(nUsedSitesPerPair, 'omitnan'));

isObject = (CC == COL_Y) | (CC == COL_P);
preMask = (Rdec.timeWindows(:,2) <= 0);
assert(any(preMask), 'No pre-stimulus bins available for baseline subtraction.');

objectTraceSite = nan(numel(useSites), nBins);
for ii = 1:numel(useSites)
    site = useSites(ii);
    if perSiteTrials
        nTrSite = double(Rdec.nTrials(site,:)).';
    else
        nTrSite = nTrialsByStim;
    end
    idxObj = isObject(site, :) & isfinite(nTrSite.');
    if ~any(idxObj)
        continue;
    end

    nObjTrials = nTrSite(idxObj);
    nObjTotal = sum(nObjTrials);
    if ~(isfinite(nObjTotal) && nObjTotal > 0)
        continue;
    end

    for tb = 1:nBins
        rObj = squeeze(Rdec.meanAct(site, idxObj, tb));
        rObj = rObj(:);
        good = isfinite(rObj) & isfinite(nObjTrials) & (nObjTrials > 0);
        if ~any(good)
            continue;
        end

        objectTraceSite(ii, tb) = sum(nObjTrials(good) .* rObj(good)) / sum(nObjTrials(good));
    end
end

objectBase = mean(objectTraceSite(:, preMask), 2, 'omitnan');
objectTraceSiteBS = bsxfun(@minus, objectTraceSite, objectBase);
if normalizeBySpontSD
    objectTraceSiteBS = bsxfun(@rdivide, objectTraceSiteBS, sdSpont(useSites));
end

visualTrace = weighted_nanmean_rows(objectTraceSiteBS, useWeights);
dtMs = median(diff(tCenters));
assert(isfinite(dtMs) && dtMs > 0, 'Could not infer a positive time step from the response windows.');
smoothBins = max(1, round(visualPeakSmoothMs / dtMs));
visualTraceSmooth = running_mean(visualTrace, smoothBins);
peakMask = (tCenters >= visualPeakWindowMs(1)) & (tCenters <= visualPeakWindowMs(2));
assert(any(peakMask), 'No decoder bins fall inside the requested visual peak window.');
peakVals = visualTraceSmooth(peakMask);
peakVals = peakVals(isfinite(peakVals));
assert(~isempty(peakVals), 'Visual peak window does not contain any finite values.');
visualPeak = max(peakVals);
assert(isfinite(visualPeak) && (visualPeak > 0), ...
    'Weighted visual peak must be positive for normalization.');

fprintf(['Using weighted visual peak %.4f for normalization ' ...
         '(smooth=%d ms, window=[%.0f %.0f] ms)\n'], ...
    visualPeak, visualPeakSmoothMs, visualPeakWindowMs(1), visualPeakWindowMs(2));

mScoreRaw = mean(scorePerPair, 1, 'omitnan');
semScoreRaw = std(scorePerPair, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(scorePerPair), 1));
mScore = mScoreRaw / visualPeak;
semScore = semScoreRaw / visualPeak;
scorePerPairNorm = scorePerPair / visualPeak;

figure; hold on;
hScore = plot(tCenters, mScore, 'LineWidth', 2);
fill([tCenters; flipud(tCenters)], [(mScore-semScore)'; flipud((mScore+semScore)')], ...
    'k', 'FaceAlpha', 0.12, 'EdgeColor', 'none');
xline(0, 'k-');
xlabel('Time from stimulus onset (ms)');
ylabel('Weighted pref - nonpref / peak visual response');
title(sprintf(['V4 weighted pair-wise color decoding ' ...
    '(%s weights |d''|, Npairs=%d, Nsites=%d, visual peak = 1)'], ...
    refWin, nPairs, numel(useSites)));
legend(hScore, sprintf('Decoder score (median sites/pair = %.0f)', median(nUsedSitesPerPair, 'omitnan')), ...
    'Location', 'best');
grid on;

if ~isempty(find(tCenters < 0, 1, 'last')) && ~isempty(find(tCenters >= 100, 1, 'first'))
    tbPre = find(tCenters < 0, 1, 'last');
    tbPost = find(tCenters >= 100, 1, 'first');
    fprintf('Mean normalized score pre (%.0f ms): %.4f\n', tCenters(tbPre), mean(scorePerPairNorm(:, tbPre), 'omitnan'));
    fprintf('Mean normalized score post (%.0f ms): %.4f\n', tCenters(tbPost), mean(scorePerPairNorm(:, tbPost), 'omitnan'));
end

ColorDecV4 = struct();
ColorDecV4.Monkey = Monkey;
ColorDecV4.monkeySuffix = monkeySuffix;
ColorDecV4.refWin = refWin;
ColorDecV4.RFrange = RFrange;
ColorDecV4.useSites = useSites;
ColorDecV4.weights = w;
ColorDecV4.weightCap = wCap;
ColorDecV4.sdSpont = sdSpont;
ColorDecV4.scorePerPair = scorePerPair;
ColorDecV4.scorePerPairNorm = scorePerPairNorm;
ColorDecV4.nUsedSitesPerPair = nUsedSitesPerPair;
ColorDecV4.meanScoreRaw = mScoreRaw;
ColorDecV4.semScoreRaw = semScoreRaw;
ColorDecV4.meanScore = mScore;
ColorDecV4.semScore = semScore;
ColorDecV4.visualTrace = visualTrace;
ColorDecV4.visualTraceSmooth = visualTraceSmooth;
ColorDecV4.visualPeak = visualPeak;
ColorDecV4.visualPeakWindowMs = visualPeakWindowMs;
ColorDecV4.visualPeakSmoothMs = visualPeakSmoothMs;
ColorDecV4.normalizeBySpontSD = normalizeBySpontSD;
ColorDecV4.timeCenters = tCenters;

function y = weighted_nanmean_rows(X, w)
X = double(X);
w = double(w(:));
assert(size(X,1) == numel(w), 'weighted_nanmean_rows: size mismatch.');

W = repmat(w, 1, size(X,2));
good = isfinite(X) & isfinite(W) & (W > 0);
X(~good) = 0;
W(~good) = 0;

den = sum(W, 1);
y = sum(X .* W, 1) ./ den;
y(den <= 0) = NaN;
end

function y = running_mean(x, win)
x = double(x(:));
win = max(1, round(win));
kernel = ones(win,1) / win;

valid = isfinite(x);
x0 = x;
x0(~valid) = 0;

num = conv(x0, kernel, 'same');
den = conv(double(valid), kernel, 'same');
y = num ./ den;
y(den <= 0) = NaN;
y = y.';
end
