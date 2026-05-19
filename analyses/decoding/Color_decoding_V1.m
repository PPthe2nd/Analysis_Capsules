%% Weighted pair-wise color decoding over time (no fitting, complementary balanced)
% Each complementary pair is a sample.
% For each pair and each time bin:
%   score = weighted mean over sites of (preferredColorResponse - nonpreferredColorResponse)
% Sites on gray are ignored (missing dimensions).
% Weights w(site) = abs(dprime(site)) from ColorTune.(refWin).dprime.
% The final population trace is normalized by the smoothed peak of the
% weighted visual response trace across all object-in-RF conditions.

% =========================
% REQUIRED INPUTS
% =========================
% R.meanAct    : [1024 x 384 x nBins]  (10 ms bins)
% R.meanSqAct  : [1024 x 384 x nBins]
% R.nTrials    : [384 x 1]
% R.timeWindows: [nBins x 2]  e.g. -200..500 in 10 ms
% Tall_V1      : [384 x 1] struct with fields stimNum and T (table 512 rows) with center_color
% ColorTune.early (and/or .late) with fields colorIndex, dprime
% Optional: SNR.sdSpont (512x1) from your previous SNR script

% =========================
% SETTINGS
% =========================
refWin = "early";      % "early" or "late" -> which ColorTune.* fields define preference + dprime weights
capWeightPct = 95;     % cap |d'| weights at this percentile to avoid a few huge weights dominating
visualPeakWindowMs = [40 250]; % search for the visual-response peak only in this post-stimulus window
visualPeakSmoothMs = 30;       % running-mean width for the visual-response trace before peak finding
useSNRsdSpontIfPresent = true;  % if SNR.sdSpont exists, use it
computeSpontSDIfMissing = true; % if no SNR.sdSpont, compute from R pre-0 bins
normalizeBySpontSD = false;     % optional legacy normalization: divide both decoder and visual-reference traces by spont SD

COL_Y = "yellowArm";
COL_P = "purple";
COL_G = "gray";

cfg = config();
tallPath = fullfile(cfg.matDir, 'Tall_V1_lines_N.mat');
colorTunePath = fullfile(cfg.matDir, 'ColorTune_balanced_V1.mat');
respPath = fullfile(cfg.matDir, 'Resp_capsules_N_d12.mat');
resp3binPath = fullfile(cfg.matDir, 'SNR_capsules_N_d12.mat');
snrNormPath = fullfile(cfg.matDir, 'SNR_V1_byColor_byWindow.mat');
hasSessionExclusions = ~isempty(site_session_exclusions("N"));

assert(exist(tallPath, 'file') == 2, ...
    'Missing %s. Run the V1 line-stimulus geometry setup first.', tallPath);
assert(exist(respPath, 'file') == 2, ...
    'Missing %s. Create the high-resolution V1 response summary first.', respPath);
assert(exist(resp3binPath, 'file') == 2, ...
    'Missing %s. Create the 3-bin V1 response summary first.', resp3binPath);

Sgeo = load(tallPath);
assert(isfield(Sgeo, 'Tall_V1') && isstruct(Sgeo.Tall_V1), ...
    '%s must contain struct Tall_V1.', tallPath);
TallDec = Sgeo.Tall_V1;

Rdec = load_capsules_struct_exclusion_aware(respPath, "N", 'cfg', cfg);
assert(isfield(Rdec, 'meanAct') && size(Rdec.meanAct, 1) >= 512, ...
    '%s must contain a response struct R with at least 512 rows.', respPath);
R3dec = load_capsules_struct_exclusion_aware(resp3binPath, "N", 'cfg', cfg);

if hasSessionExclusions
    fprintf('Session exclusions are active for monkey N; recomputing V1 SNR/color tuning from exclusion-aware responses.\n');
end

SNRnorm = struct();
if hasSessionExclusions || exist(snrNormPath, 'file') ~= 2
    SNRnorm = compute_snr_per_color_sites(R3dec, TallDec, 1:512, 'Verbose', false);
elseif normalizeBySpontSD && useSNRsdSpontIfPresent && exist(snrNormPath, 'file') == 2
    Ssnr = load(snrNormPath);
    if isfield(Ssnr, 'SNR') && isstruct(Ssnr.SNR)
        SNRnorm = Ssnr.SNR;
    end
end

if exist(colorTunePath, 'file') == 2 && ~hasSessionExclusions
    Sct = load(colorTunePath);
    assert(isfield(Sct, 'ColorTune') && isstruct(Sct.ColorTune), ...
        '%s must contain struct ColorTune.', colorTunePath);
    ColorTuneDec = Sct.ColorTune;
else
    SNRmat3 = [SNRnorm.yellowEarly(1:512), SNRnorm.yellowLate(1:512), ...
               SNRnorm.purpleEarly(1:512), SNRnorm.purpleLate(1:512)];
    bestSNR3 = max(SNRmat3, [], 2, 'omitnan');
    keepSiteIdx3 = find(isfinite(bestSNR3) & (bestSNR3 > 0.7));
    assert(~isempty(keepSiteIdx3), ...
        'No V1 sites passed the SNR threshold %.2f during ColorTune construction.', 0.7);
    ColorTuneDec = compute_color_tuning_balanced_sites(R3dec, TallDec, (1:512).', keepSiteIdx3, 'Verbose', false);
    ColorTuneDec.thr = 0.7;
    ColorTuneDec.bestSNR = bestSNR3;
end

% =========================
% BASIC DIMS
% =========================
[nCh, nStim2, nBins] = size(Rdec.meanAct);
nStim = nStim2;
assert(nStim == 384, 'Expected 384 stimuli.');
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

v1Sites = 1:512;      % mapping: site i -> channel i in R
tCenters = mean(Rdec.timeWindows, 2);  % nBins x 1

% =========================
% Choose preference + weights from ColorTune
% =========================
switch lower(refWin)
    case "early"
        ci = ColorTuneDec.early.colorIndex(:);   % 512x1
        dp = ColorTuneDec.early.dprime(:);       % 512x1
    case "late"
        ci = ColorTuneDec.late.colorIndex(:);
        dp = ColorTuneDec.late.dprime(:);
    otherwise
        error('refWin must be "early" or "late".');
end

% Preference direction: CI>0 => prefers yellow, CI<0 => prefers purple
prefIsYellow = (ci > 0);
prefIsYellow(~isfinite(ci)) = false; % arbitrary, but those will get ~0 weight anyway

% Weights: abs(d')
w = abs(dp);
w(~isfinite(w)) = 0;

% Cap extreme weights robustly
wPos = w(w>0);
if ~isempty(wPos)
    wCap = prctile(wPos, capWeightPct);
    w95 = prctile(wPos, 95);
    w = min(w, wCap);
else
    wCap = NaN;
    w95 = NaN;
end

fprintf('Weight summary |d''|: median=%.3f, 95%%=%.3f, max(after cap)=%.3f\n', ...
    median(w(w>0),'omitnan'), w95, max(w));

% If a site has 0 weight, it contributes nothing (no need to threshold by CI)
useSites = find(w > 0);
fprintf('Using %d sites with nonzero weight\n', numel(useSites));
useWeights = w(useSites);

% =========================
% Sort Tall_V1 by stimNum and extract center_color per site x stim
% =========================
TallStimNums = arrayfun(@(x) x.stimNum, TallDec(:));
[sortedStimNums, order] = sort(TallStimNums(:));
assert(all(sortedStimNums(:).' == 1:nStim), 'Tall_V1.stimNum should cover 1..384.');
TallSorted = TallDec(order);

T0 = TallSorted(1).T;
assert(istable(T0), 'Tall_V1(stim).T must be a table.');
vn = string(T0.Properties.VariableNames);
colIdx = find(vn=="center_color", 1);
if isempty(colIdx)
    colIdx = find(contains(lower(vn), "center") & contains(lower(vn), "color"), 1);
end
assert(~isempty(colIdx), 'Could not find center_color column in Tall_V1(stim).T.');

CC = strings(512, nStim);
for stim = 1:nStim
    CC(:,stim) = strtrim(string(TallSorted(stim).T{:, colIdx})); % 512x1
end

% =========================
% Build complementary pairs (1<->5),(2<->6),(3<->7),(4<->8), repeating
% =========================
pairsA = []; pairsB = [];
for a = 1:nStim
    pos = mod(a-1,8) + 1;
    if pos <= 4
        pairsA(end+1,1) = a; %#ok<AGROW>
        pairsB(end+1,1) = a + 4; %#ok<AGROW>
    end
end
nPairs = numel(pairsA);
assert(nPairs == 192, 'Expected 192 complementary pairs, got %d.', nPairs);
fprintf('Using %d complementary pairs\n', nPairs);

% =========================
% Spontaneous SD per site (for normalization)
% =========================
sdSpont = ones(512,1);

if normalizeBySpontSD
    if useSNRsdSpontIfPresent && isfield(SNRnorm,'sdSpont') && numel(SNRnorm.sdSpont) >= 512
        sdSpont = SNRnorm.sdSpont(1:512);
        sdSpont(~isfinite(sdSpont) | sdSpont<=0) = 1;
        fprintf('Using sdSpont from SNR_V1_byColor_byWindow.mat\n');

    elseif computeSpontSDIfMissing
        % Compute from pre-0 bins in the 10 ms representation
        isSpontBin = (Rdec.timeWindows(:,2) <= 0);
        assert(any(isSpontBin), 'No spontaneous bins found (timeWindows end<=0).');

        bIdx = find(isSpontBin);

        fprintf('Computing sdSpont from %d pre-0 bins using R.meanAct/meanSqAct...\n', numel(bIdx));

        for site = 1:512
            ch = v1Sites(site);
            if perSiteTrials
                nTrSite = double(Rdec.nTrials(ch,:)).';
            else
                nTrSite = nTrialsByStim;
            end
            wStim = nTrSite;
            wStim = wStim / sum(wStim);

            varBins = nan(numel(bIdx),1);
            for ib = 1:numel(bIdx)
                tb = bIdx(ib);

                m  = squeeze(Rdec.meanAct(ch,:,tb)).';    % 384x1
                m2 = squeeze(Rdec.meanSqAct(ch,:,tb)).';

                good = isfinite(m) & isfinite(m2) & (nTrSite > 0);
                if ~any(good), continue; end

                ww = wStim(good); ww = ww / sum(ww);
                mu = sum(ww .* m(good));
                ex2 = sum(ww .* m2(good));
                varBins(ib) = max(0, ex2 - mu^2);
            end

            sd = sqrt(mean(varBins, 'omitnan'));
            if isfinite(sd) && sd>0
                sdSpont(site) = sd;
            else
                sdSpont(site) = 1;
            end
        end
        fprintf('Done computing sdSpont\n');
    else
        fprintf('normalizeBySpontSD=true but no sdSpont available; using sdSpont=1\n');
    end
end

% =========================
% MAIN: score per pair per time bin (weighted)
% =========================
scorePerPair = nan(nPairs, nBins);
nUsedSitesPerPair = nan(nPairs,1);

for ip = 1:nPairs
    a = pairsA(ip);
    b = pairsB(ip);

    % Pull responses for all V1 channels for these stimuli
    Ra = squeeze(Rdec.meanAct(v1Sites, a, :)); % 512 x nBins
    Rb = squeeze(Rdec.meanAct(v1Sites, b, :)); % 512 x nBins

    sumScore = zeros(1, nBins);
    sumW     = zeros(1, nBins);

    nUsed = 0;

    for site = useSites(:).'
        ca = CC(site,a);
        cb = CC(site,b);

        % Ignore if either member is gray (missing) or unknown
        if ca==COL_G || cb==COL_G || ca=="" || cb==""
            continue;
        end

        % Determine which stim has Yellow and which has Purple for this site
        if ca==COL_Y && cb==COL_P
            rY = Ra(site,:); rP = Rb(site,:);
        elseif ca==COL_P && cb==COL_Y
            rY = Rb(site,:); rP = Ra(site,:);
        else
            % unexpected label combo (e.g. not a clean Y/P swap), skip
            continue;
        end

        wi = w(site);
        if wi <= 0
            continue;
        end

        % Signed difference according to preference
        if prefIsYellow(site)
            d = (rY - rP);
        else
            d = (rP - rY);
        end

        % Normalize by spont SD (optional legacy mode)
        if normalizeBySpontSD
            d = d ./ sdSpont(site);
        end

        sumScore = sumScore + wi * d;
        sumW     = sumW + wi;

        nUsed = nUsed + 1;
    end

    nUsedSitesPerPair(ip) = nUsed;

    ok = sumW > 0;
    tmp = nan(1, nBins);
    tmp(ok) = sumScore(ok) ./ sumW(ok);
    scorePerPair(ip,:) = tmp;
end

fprintf('Median #sites contributing per pair: %.1f\n', median(nUsedSitesPerPair, 'omitnan'));

% =========================
% Weighted visual-response normalization
% =========================
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

W = repmat(useWeights, 1, nBins);
goodObj = isfinite(objectTraceSiteBS) & isfinite(W) & (W > 0);
Xobj = objectTraceSiteBS;
Xobj(~goodObj) = 0;
W(~goodObj) = 0;
denObj = sum(W, 1);
visualTrace = sum(Xobj .* W, 1) ./ denObj;
visualTrace(denObj <= 0) = NaN;

dtMs = median(diff(tCenters));
assert(isfinite(dtMs) && dtMs > 0, 'Could not infer a positive time step from the response windows.');
smoothBins = max(1, round(visualPeakSmoothMs / dtMs));
kernel = ones(smoothBins,1) / smoothBins;
validVis = isfinite(visualTrace(:));
vis0 = visualTrace(:);
vis0(~validVis) = 0;
numVis = conv(vis0, kernel, 'same');
denVis = conv(double(validVis), kernel, 'same');
visualTraceSmooth = numVis ./ denVis;
visualTraceSmooth(denVis <= 0) = NaN;
visualTraceSmooth = visualTraceSmooth.';

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

% =========================
% Average across pairs
% =========================
mScoreRaw = mean(scorePerPair, 1, 'omitnan');
semScoreRaw = std(scorePerPair, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(scorePerPair),1));
mScore = mScoreRaw / visualPeak;
semScore = semScoreRaw / visualPeak;
scorePerPairNorm = scorePerPair / visualPeak;

figure; hold on;
plot(tCenters, mScore, 'LineWidth', 2);
fill([tCenters; flipud(tCenters)], [(mScore-semScore)'; flipud((mScore+semScore)')], ...
     'k', 'FaceAlpha', 0.12, 'EdgeColor', 'none');
xline(0,'k-');
xlabel('Time from stimulus onset (ms)');
ylabel('Weighted pref - nonpref / peak visual response');
title(sprintf(['Weighted pair-wise color decoding ' ...
    '(%s weights |d''|, Npairs=%d, Nsites=%d, visual peak = 1)'], ...
    refWin, nPairs, numel(useSites)));
grid on;

% Optional: quick pre/post sanity
tbPre  = find(tCenters < 0, 1, 'last');
tbPost = find(tCenters >= 100, 1, 'first');
if ~isempty(tbPre) && ~isempty(tbPost)
    fprintf('Mean normalized score pre (%.0f ms): %.4f\n', tCenters(tbPre), mean(scorePerPairNorm(:,tbPre),'omitnan'));
    fprintf('Mean normalized score post (%.0f ms): %.4f\n', tCenters(tbPost), mean(scorePerPairNorm(:,tbPost),'omitnan'));
end
