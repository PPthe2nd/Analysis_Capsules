function OUT = Attention_ColorPref_Line_Stimuli()
% ATTENTION_COLORPREF_LINE_STIMULI
% Time-resolved responses for 4 conditions on the intersection of:
%   - attention-significant sites (OUT3 pValueTD < pAttentionThresh)
%   - color-significant sites (ColorTune.<window>.p < pColorThresh)
%
% Conditions (per site):
%   1) best color attended      (target + preferred color)
%   2) best color non-attended  (distractor + preferred color)
%   3) worst color attended     (target + non-preferred color)
%   4) worst color non-attended (distractor + non-preferred color)
%
% Preferred color is defined by sign(ColorTune.<window>.colorIndex):
%   > 0 => yellow preferred, < 0 => purple preferred.

%% Settings
P = struct();
P.respCapsulesFile = 'Resp_capsules_N_d20.mat';
P.resultTag = 'd20';

P.attention3binFile = 'SNR_capsules_N_d12.mat';
P.attentionTimeIdx = 3;
P.pAttentionThresh = 0.05;
P.excludeOverlap = true;

P.colorTuneFile = 'ColorTune_balanced_V1.mat';
P.colorTuneWindow = 'early';      % 'early' or 'late'
P.pColorThresh = 0.05;

P.snrNormFile = 'SNR_V1_byColor_byWindow.mat';
P.v1Sites = 1:512;
P.normalizeResponses = true;      % normalize using (x-muSpont)/(muTop-muSpont)
P.minTrialsPerCondition = 1;      % minimum total trial count per condition per site

P.plotSem = true;
P.plotFigure = true;
P.saveResult = true;
P.plotInteractionFigure = true;
P.interactionAlpha = 0.05;
P.interactionWindowMs = [100 400];
P.interactionTest = 'signrank';  % 'signrank' or 'ttest'
P.interactionFdr = true;         % Benjamini-Hochberg over bins

cfg = config();

%% Load required data
S = load(fullfile(cfg.matDir, 'Tall_V1_lines_N.mat'));
assert(isfield(S, 'Tall_V1') && isstruct(S.Tall_V1), ...
    'Tall_V1_lines_N.mat must contain struct Tall_V1.');
Tall_V1 = S.Tall_V1;

S = load(fullfile(cfg.matDir, P.respCapsulesFile));
assert(isfield(S, 'R') && isstruct(S.R), '%s must contain struct R.', P.respCapsulesFile);
R_resp = S.R;
assert(all(isfield(R_resp, {'meanAct','nTrials','timeWindows'})), ...
    'R_resp must contain meanAct, nTrials, timeWindows.');

S = load(fullfile(cfg.matDir, P.colorTuneFile));
assert(isfield(S, 'ColorTune') && isstruct(S.ColorTune), ...
    '%s must contain struct ColorTune.', P.colorTuneFile);
ColorTune = S.ColorTune;
assert(isfield(ColorTune, P.colorTuneWindow), ...
    'ColorTune missing window "%s".', P.colorTuneWindow);

S = load(fullfile(cfg.matDir, P.snrNormFile));
assert(isfield(S, 'SNR') && isstruct(S.SNR), '%s must contain struct SNR.', P.snrNormFile);
SNR = S.SNR;

%% Load or compute attention baseline OUT3 (used only for significance mask)
out3File = fullfile(cfg.resultsDir, 'OUT_attention_modulation_3bin_timeIdx3.mat');
if exist(out3File, 'file') == 2
    S = load(out3File, 'OUT');
    assert(isfield(S, 'OUT') && isstruct(S.OUT), ...
        'File %s must contain struct OUT.', out3File);
    OUT3 = S.OUT;
else
    S = load(fullfile(cfg.matDir, P.attention3binFile));
    assert(isfield(S, 'R') && isstruct(S.R), ...
        '%s must contain struct R.', P.attention3binFile);
    R3 = S.R;
    optsTD = struct('timeIdx', P.attentionTimeIdx, 'excludeOverlap', P.excludeOverlap, 'verbose', true);
    OUT3 = attention_modulation_V1_3bin(R3, Tall_V1, SNR, optsTD);
    OUT = OUT3; %#ok<NASGU>
    meta = struct('created', datestr(now, 30), 'timeIdx', P.attentionTimeIdx); %#ok<NASGU>
    save(out3File, 'OUT', 'meta', '-v7.3');
end

%% Basic dimensions and site indexing
[nCh, nStim, nBins] = size(R_resp.meanAct);
assert(nStim == 384, 'Expected 384 stimuli in %s.', P.respCapsulesFile);
assert(size(R_resp.timeWindows,1) == nBins && size(R_resp.timeWindows,2) == 2, ...
    'R_resp.timeWindows must be [nBins x 2].');

v1Sites = P.v1Sites(:)';
nV1 = numel(v1Sites);
assert(max(v1Sites) <= nCh, 'Requested V1 site index exceeds channels in R_resp.meanAct.');
assert(numel(OUT3.pValueTD) >= max(v1Sites), 'OUT3.pValueTD smaller than max(v1Sites).');
assert(numel(ColorTune.(P.colorTuneWindow).p) >= max(v1Sites), ...
    'ColorTune.%s.p smaller than max(v1Sites).', P.colorTuneWindow);
assert(numel(ColorTune.(P.colorTuneWindow).colorIndex) >= max(v1Sites), ...
    'ColorTune.%s.colorIndex smaller than max(v1Sites).', P.colorTuneWindow);

tCenters = mean(double(R_resp.timeWindows), 2);

%% Site selection: intersection of attention and color significance
pAtt = double(OUT3.pValueTD(v1Sites));
pCol = double(ColorTune.(P.colorTuneWindow).p(v1Sites));
ci = double(ColorTune.(P.colorTuneWindow).colorIndex(v1Sites));

isAttSig = isfinite(pAtt) & (pAtt < P.pAttentionThresh);
isColorSig = isfinite(pCol) & (pCol < P.pColorThresh);
prefIsYellow = ci > 0;

isBoth = isAttSig & isColorSig & isfinite(ci) & (ci ~= 0);
siteLocalAll = find(isBoth);      % indices in 1..nV1
siteGlobalAll = v1Sites(siteLocalAll);

fprintf('Attention-significant sites: %d / %d (p < %.3f)\n', nnz(isAttSig), nV1, P.pAttentionThresh);
fprintf('Color-significant sites: %d / %d (p < %.3f, window=%s)\n', ...
    nnz(isColorSig), nV1, P.pColorThresh, P.colorTuneWindow);
fprintf('Intersection sites (attention & color): %d / %d\n', numel(siteLocalAll), nV1);

%% Sort Tall_V1 by stimulus number and build per-site/stim masks
stimNums = arrayfun(@(x) x.stimNum, Tall_V1(:));
[stimNumsSorted, ord] = sort(stimNums(:));
assert(numel(stimNumsSorted) == nStim && all(stimNumsSorted(:).' == 1:nStim), ...
    'Tall_V1.stimNum must cover 1..%d exactly.', nStim);
TallSorted = Tall_V1(ord);

isTarget = false(nV1, nStim);
isDistr  = false(nV1, nStim);
isBG     = false(nV1, nStim);
isY      = false(nV1, nStim);
isP      = false(nV1, nStim);
isOV     = false(nV1, nStim);

for s = 1:nStim
    assert(isfield(TallSorted(s), 'T') && istable(TallSorted(s).T), ...
        'Tall_V1(%d).T missing or not a table.', s);
    T = TallSorted(s).T;
    assert(height(T) >= max(v1Sites), ...
        'Tall_V1(%d).T has %d rows; need at least %d.', s, height(T), max(v1Sites));
    vn = string(T.Properties.VariableNames);
    assert(ismember("assignment", vn) && ismember("center_color", vn), ...
        'Tall_V1(%d).T must contain assignment and center_color.', s);

    asg = toCellStrLocal(T.assignment(v1Sites));
    col = toCellStrLocal(T.center_color(v1Sites));

    isTarget(:,s) = strcmpi(asg, 'target');
    isDistr(:,s)  = strcmpi(asg, 'distractor');
    isBG(:,s)     = strcmpi(asg, 'background');
    isY(:,s)      = cellfun(@(x) contains(lower(x), 'yellow'), col);
    isP(:,s)      = cellfun(@(x) contains(lower(x), 'purple'), col);

    if ismember('overlap', string(T.Properties.VariableNames))
        ov = double(T.overlap(v1Sites));
        isOV(:,s) = (ov ~= 0);
    end
end

baseMask = ~isBG & (isY | isP);
if P.excludeOverlap
    baseMask = baseMask & ~isOV;
end

%% Prepare response normalization
if P.normalizeResponses
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
    bAll = zeros(nV1,1);
    scaleAll = ones(nV1,1);
end

%% Trial-count handling
nTrialsRaw = R_resp.nTrials;
if isvector(nTrialsRaw)
    nTrialsByStim = double(nTrialsRaw(:)');
    perSiteTrials = false;
    assert(numel(nTrialsByStim) == nStim, 'R_resp.nTrials vector must have %d elements.', nStim);
elseif ismatrix(nTrialsRaw) && (size(nTrialsRaw,2) == nStim)
    perSiteTrials = true;
    nTrialsByStim = [];
    assert(size(nTrialsRaw,1) >= max(v1Sites), ...
        'R_resp.nTrials has %d rows; need at least %d.', size(nTrialsRaw,1), max(v1Sites));
else
    error('R_resp.nTrials must be vector(384) or matrix(nChannels x 384).');
end

%% Build per-site 4-condition time courses
nSiteCand = numel(siteLocalAll);
tcBestAtt = nan(nSiteCand, nBins);
tcBestNon = nan(nSiteCand, nBins);
tcWorstAtt = nan(nSiteCand, nBins);
tcWorstNon = nan(nSiteCand, nBins);

nTrBestAtt = zeros(nSiteCand,1);
nTrBestNon = zeros(nSiteCand,1);
nTrWorstAtt = zeros(nSiteCand,1);
nTrWorstNon = zeros(nSiteCand,1);

nStimBestAtt = zeros(nSiteCand,1);
nStimBestNon = zeros(nSiteCand,1);
nStimWorstAtt = zeros(nSiteCand,1);
nStimWorstNon = zeros(nSiteCand,1);

for k = 1:nSiteCand
    iLocal = siteLocalAll(k);
    ch = v1Sites(iLocal);

    if perSiteTrials
        nTr = double(nTrialsRaw(ch,:));
    else
        nTr = nTrialsByStim;
    end
    nTr(~isfinite(nTr) | nTr < 0) = 0;

    if prefIsYellow(iLocal)
        isPref = isY(iLocal,:);
        isNonpref = isP(iLocal,:);
    else
        isPref = isP(iLocal,:);
        isNonpref = isY(iLocal,:);
    end

    mBase = baseMask(iLocal,:);
    mBestAtt = mBase & isTarget(iLocal,:) & isPref;
    mBestNon = mBase & isDistr(iLocal,:)  & isPref;
    mWorstAtt = mBase & isTarget(iLocal,:) & isNonpref;
    mWorstNon = mBase & isDistr(iLocal,:)  & isNonpref;

    resp = squeeze(double(R_resp.meanAct(ch,:,:))); % [nStim x nBins]
    if size(resp,1) ~= nStim && size(resp,2) == nStim
        resp = resp.';
    end
    assert(size(resp,1) == nStim && size(resp,2) == nBins, ...
        'Unexpected response shape at channel %d.', ch);

    [tcBestAtt(k,:), nTrBestAtt(k), nStimBestAtt(k)] = weightedTimeCourseLocal(resp, nTr, mBestAtt);
    [tcBestNon(k,:), nTrBestNon(k), nStimBestNon(k)] = weightedTimeCourseLocal(resp, nTr, mBestNon);
    [tcWorstAtt(k,:), nTrWorstAtt(k), nStimWorstAtt(k)] = weightedTimeCourseLocal(resp, nTr, mWorstAtt);
    [tcWorstNon(k,:), nTrWorstNon(k), nStimWorstNon(k)] = weightedTimeCourseLocal(resp, nTr, mWorstNon);

    if P.normalizeResponses
        b = bAll(iLocal);
        sc = scaleAll(iLocal);
        if isfinite(sc) && (sc > 0)
            tcBestAtt(k,:) = (tcBestAtt(k,:) - b) ./ sc;
            tcBestNon(k,:) = (tcBestNon(k,:) - b) ./ sc;
            tcWorstAtt(k,:) = (tcWorstAtt(k,:) - b) ./ sc;
            tcWorstNon(k,:) = (tcWorstNon(k,:) - b) ./ sc;
        else
            tcBestAtt(k,:) = nan(1, nBins);
            tcBestNon(k,:) = nan(1, nBins);
            tcWorstAtt(k,:) = nan(1, nBins);
            tcWorstNon(k,:) = nan(1, nBins);
            nTrBestAtt(k) = 0; nTrBestNon(k) = 0;
            nTrWorstAtt(k) = 0; nTrWorstNon(k) = 0;
            nStimBestAtt(k) = 0; nStimBestNon(k) = 0;
            nStimWorstAtt(k) = 0; nStimWorstNon(k) = 0;
        end
    end
end

hasAll4 = (nTrBestAtt >= P.minTrialsPerCondition) & ...
          (nTrBestNon >= P.minTrialsPerCondition) & ...
          (nTrWorstAtt >= P.minTrialsPerCondition) & ...
          (nTrWorstNon >= P.minTrialsPerCondition);

siteLocal = siteLocalAll(hasAll4);
siteGlobal = siteGlobalAll(hasAll4);

tcBestAtt = tcBestAtt(hasAll4,:);
tcBestNon = tcBestNon(hasAll4,:);
tcWorstAtt = tcWorstAtt(hasAll4,:);
tcWorstNon = tcWorstNon(hasAll4,:);

nTrBestAtt = nTrBestAtt(hasAll4);
nTrBestNon = nTrBestNon(hasAll4);
nTrWorstAtt = nTrWorstAtt(hasAll4);
nTrWorstNon = nTrWorstNon(hasAll4);

nStimBestAtt = nStimBestAtt(hasAll4);
nStimBestNon = nStimBestNon(hasAll4);
nStimWorstAtt = nStimWorstAtt(hasAll4);
nStimWorstNon = nStimWorstNon(hasAll4);

fprintf('Sites with all 4 conditions present: %d / %d intersection sites\n', ...
    numel(siteLocal), nSiteCand);
if ~isempty(siteLocal)
    fprintf(['Median trial totals per site: best-att=%g, best-nonatt=%g, ' ...
             'worst-att=%g, worst-nonatt=%g\n'], ...
        median(nTrBestAtt), median(nTrBestNon), median(nTrWorstAtt), median(nTrWorstNon));
end

%% Population means and SEM
mBestAtt = mean(tcBestAtt, 1, 'omitnan');
mBestNon = mean(tcBestNon, 1, 'omitnan');
mWorstAtt = mean(tcWorstAtt, 1, 'omitnan');
mWorstNon = mean(tcWorstNon, 1, 'omitnan');

semBestAtt = std(tcBestAtt, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(tcBestAtt),1));
semBestNon = std(tcBestNon, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(tcBestNon),1));
semWorstAtt = std(tcWorstAtt, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(tcWorstAtt),1));
semWorstNon = std(tcWorstNon, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(tcWorstNon),1));

%% Additive vs multiplicative interaction test
% Interaction (difference-of-differences):
%   I(t) = [bestAtt-bestNon](t) - [worstAtt-worstNon](t)
attEffectBest = tcBestAtt - tcBestNon;
attEffectWorst = tcWorstAtt - tcWorstNon;
interactionBySite = attEffectBest - attEffectWorst;

mAttBest = mean(attEffectBest, 1, 'omitnan');
mAttWorst = mean(attEffectWorst, 1, 'omitnan');
semAttBest = std(attEffectBest, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(attEffectBest),1));
semAttWorst = std(attEffectWorst, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(attEffectWorst),1));

mInteraction = mean(interactionBySite, 1, 'omitnan');
semInteraction = std(interactionBySite, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(interactionBySite),1));

pBin = nan(1, nBins);
nBin = zeros(1, nBins);
testUsedBin = strings(1, nBins);
for b = 1:nBins
    xb = interactionBySite(:, b);
    xb = xb(isfinite(xb));
    nBin(b) = numel(xb);
    if numel(xb) < 3
        continue;
    end
    [pBin(b), testUsedBin(b)] = oneSamplePLocal(xb, P.interactionAlpha, P.interactionTest);
end

if P.interactionFdr
    [qBin, sigBin] = fdrBhLocal(pBin, P.interactionAlpha);
else
    qBin = nan(size(pBin));
    sigBin = isfinite(pBin) & (pBin < P.interactionAlpha);
end

wMask = (tCenters >= P.interactionWindowMs(1)) & (tCenters <= P.interactionWindowMs(2));
IwinSite = mean(interactionBySite(:, wMask), 2, 'omitnan');
IwinSite = IwinSite(isfinite(IwinSite));
nWin = numel(IwinSite);
if nWin >= 3
    [pWin, testUsedWin] = oneSamplePLocal(IwinSite, P.interactionAlpha, P.interactionTest);
else
    pWin = NaN;
    testUsedWin = "n<3";
end
meanIwin = mean(IwinSite, 'omitnan');
medianIwin = median(IwinSite, 'omitnan');

if isfinite(pWin) && (pWin < P.interactionAlpha)
    if meanIwin > 0
        interactionCall = "multiplicative_like";
        interactionNote = "Attention effect is larger for best than worst color.";
    elseif meanIwin < 0
        interactionCall = "inverse_interaction";
        interactionNote = "Attention effect is smaller for best than worst color.";
    else
        interactionCall = "interaction_nonzero_unclear_direction";
        interactionNote = "Interaction is significant but near zero mean.";
    end
else
    interactionCall = "additive_consistent";
    interactionNote = "No significant interaction; consistent with additive attention effect.";
end

fprintf('\nInteraction test (additive vs multiplicative):\n');
fprintf('  Window [%g %g] ms: N=%d, mean(I)=%.6g, median(I)=%.6g, p=%.4g (%s)\n', ...
    P.interactionWindowMs(1), P.interactionWindowMs(2), nWin, meanIwin, medianIwin, pWin, testUsedWin);
if P.interactionFdr
    fprintf('  Bin-wise significant bins (FDR q<%.3f): %d / %d\n', ...
        P.interactionAlpha, nnz(sigBin), nBins);
else
    fprintf('  Bin-wise significant bins (p<%.3f): %d / %d\n', ...
        P.interactionAlpha, nnz(sigBin), nBins);
end
fprintf('  Call: %s (%s)\n\n', interactionCall, interactionNote);

%% Plot
if P.plotFigure
    cBestAtt = [0.10 0.55 0.10];
    cBestNon = [0.35 0.75 0.35];
    cWorstAtt = [0.55 0.15 0.70];
    cWorstNon = [0.80 0.45 0.90];

    hFig = figure('Color', 'w');
    ax = axes('Parent', hFig); %#ok<LAXES>
    hold(ax, 'on');

    if P.plotSem
        addSemPatchLocal(ax, tCenters, mBestAtt, semBestAtt, cBestAtt, 0.16);
        addSemPatchLocal(ax, tCenters, mBestNon, semBestNon, cBestNon, 0.12);
        addSemPatchLocal(ax, tCenters, mWorstAtt, semWorstAtt, cWorstAtt, 0.16);
        addSemPatchLocal(ax, tCenters, mWorstNon, semWorstNon, cWorstNon, 0.12);
    end

    plot(ax, tCenters, mBestAtt, '-',  'Color', cBestAtt, 'LineWidth', 2.2);
    plot(ax, tCenters, mBestNon, '--', 'Color', cBestNon, 'LineWidth', 2.0);
    plot(ax, tCenters, mWorstAtt, '-',  'Color', cWorstAtt, 'LineWidth', 2.2);
    plot(ax, tCenters, mWorstNon, '--', 'Color', cWorstNon, 'LineWidth', 2.0);
    xline(ax, 0, 'k-');

    xlabel(ax, 'Time from stimulus onset (ms)');
    if P.normalizeResponses
        ylabel(ax, 'Normalized response (a.u.)');
    else
        ylabel(ax, 'Response (a.u.)');
    end
    title(ax, sprintf('Best/Worst color x Attention (N=%d sites)', numel(siteLocal)));
    legend(ax, {'Best color attended', 'Best color non-attended', ...
                'Worst color attended', 'Worst color non-attended'}, ...
           'Location', 'best');
    grid(ax, 'on');
end

if P.plotInteractionFigure
    cBest = [0.10 0.55 0.10];
    cWorst = [0.55 0.15 0.70];
    cInter = [0.15 0.15 0.15];

    hInt = figure('Color', 'w');
    figure(hInt);
    ax1 = subplot(2,1,1);
    hold(ax1, 'on');
    addSemPatchLocal(ax1, tCenters, mAttBest, semAttBest, cBest, 0.14);
    addSemPatchLocal(ax1, tCenters, mAttWorst, semAttWorst, cWorst, 0.14);
    plot(ax1, tCenters, mAttBest, '-', 'Color', cBest, 'LineWidth', 2.0);
    plot(ax1, tCenters, mAttWorst, '-', 'Color', cWorst, 'LineWidth', 2.0);
    xline(ax1, 0, 'k-');
    yline(ax1, 0, 'k:');
    ylabel(ax1, '\Delta attention (A-N)');
    title(ax1, 'Attention effect by color level');
    legend(ax1, {'Best color', 'Worst color'}, 'Location', 'best');
    grid(ax1, 'on');

    ax2 = subplot(2,1,2);
    hold(ax2, 'on');
    addSemPatchLocal(ax2, tCenters, mInteraction, semInteraction, cInter, 0.18);
    plot(ax2, tCenters, mInteraction, '-', 'Color', cInter, 'LineWidth', 2.2);
    xline(ax2, 0, 'k-');
    yline(ax2, 0, 'k:');
    if any(sigBin)
        tmpLo = [0; mInteraction(:)-semInteraction(:)];
        tmpLo = tmpLo(isfinite(tmpLo));
        if isempty(tmpLo)
            yMark = 0;
        else
            yMark = min(tmpLo);
        end
        tmpSpan = [mInteraction(:)-semInteraction(:); mInteraction(:)+semInteraction(:)];
        tmpSpan = tmpSpan(isfinite(tmpSpan));
        if numel(tmpSpan) >= 2
            dy = max(tmpSpan) - min(tmpSpan);
        else
            dy = 1;
        end
        yMark = yMark - 0.06 * dy;
        plot(ax2, tCenters(sigBin), yMark * ones(1, nnz(sigBin)), 'o', ...
            'MarkerSize', 4, 'MarkerFaceColor', [0.9 0.2 0.2], 'MarkerEdgeColor', 'none');
    end
    xlabel(ax2, 'Time from stimulus onset (ms)');
    ylabel(ax2, 'Interaction I');
    if P.interactionFdr
        ttl2 = sprintf('I=(BA-BN)-(WA-WN), q<%.2g bins: %d', P.interactionAlpha, nnz(sigBin));
    else
        ttl2 = sprintf('I=(BA-BN)-(WA-WN), p<%.2g bins: %d', P.interactionAlpha, nnz(sigBin));
    end
    title(ax2, ttl2);
    grid(ax2, 'on');
end

%% Package outputs
OUT = struct();
OUT.params = P;
OUT.timeCentersMs = tCenters(:);
OUT.siteIdxLocal = siteLocal(:);
OUT.siteIdxGlobal = siteGlobal(:);
OUT.prefIsYellow = prefIsYellow(siteLocal(:));

OUT.tcBestAtt = tcBestAtt;
OUT.tcBestNon = tcBestNon;
OUT.tcWorstAtt = tcWorstAtt;
OUT.tcWorstNon = tcWorstNon;

OUT.nTrialsBestAtt = nTrBestAtt;
OUT.nTrialsBestNon = nTrBestNon;
OUT.nTrialsWorstAtt = nTrWorstAtt;
OUT.nTrialsWorstNon = nTrWorstNon;

OUT.nStimBestAtt = nStimBestAtt;
OUT.nStimBestNon = nStimBestNon;
OUT.nStimWorstAtt = nStimWorstAtt;
OUT.nStimWorstNon = nStimWorstNon;

OUT.meanBestAtt = mBestAtt;
OUT.meanBestNon = mBestNon;
OUT.meanWorstAtt = mWorstAtt;
OUT.meanWorstNon = mWorstNon;

OUT.semBestAtt = semBestAtt;
OUT.semBestNon = semBestNon;
OUT.semWorstAtt = semWorstAtt;
OUT.semWorstNon = semWorstNon;

OUT.interaction = struct();
OUT.interaction.attEffectBest = attEffectBest;
OUT.interaction.attEffectWorst = attEffectWorst;
OUT.interaction.interactionBySite = interactionBySite;
OUT.interaction.meanAttEffectBest = mAttBest;
OUT.interaction.meanAttEffectWorst = mAttWorst;
OUT.interaction.semAttEffectBest = semAttBest;
OUT.interaction.semAttEffectWorst = semAttWorst;
OUT.interaction.meanInteraction = mInteraction;
OUT.interaction.semInteraction = semInteraction;
OUT.interaction.pBin = pBin;
OUT.interaction.qBin = qBin;
OUT.interaction.sigBin = sigBin;
OUT.interaction.nPerBin = nBin;
OUT.interaction.testUsedPerBin = testUsedBin;
OUT.interaction.windowMs = P.interactionWindowMs;
OUT.interaction.windowMask = wMask(:);
OUT.interaction.windowSiteInteraction = IwinSite(:);
OUT.interaction.nWindowSites = nWin;
OUT.interaction.pWindow = pWin;
OUT.interaction.testUsedWindow = testUsedWin;
OUT.interaction.meanWindow = meanIwin;
OUT.interaction.medianWindow = medianIwin;
OUT.interaction.call = interactionCall;
OUT.interaction.note = interactionNote;
OUT.interaction.alpha = P.interactionAlpha;
OUT.interaction.fdrApplied = logical(P.interactionFdr);

OUT.siteSelection = struct();
OUT.siteSelection.pAttention = pAtt(:);
OUT.siteSelection.pColor = pCol(:);
OUT.siteSelection.colorIndex = ci(:);
OUT.siteSelection.isAttentionSig = isAttSig(:);
OUT.siteSelection.isColorSig = isColorSig(:);
OUT.siteSelection.isBoth = isBoth(:);

%% Save
if isempty(strtrim(P.resultTag))
    [~, baseResp] = fileparts(P.respCapsulesFile);
    resultTag = regexprep(baseResp, '[^A-Za-z0-9_-]', '_');
else
    resultTag = regexprep(P.resultTag, '[^A-Za-z0-9_-]', '_');
end
outFile = fullfile(cfg.resultsDir, sprintf('attention_colorpref_timeseries_%s.mat', resultTag));

if P.saveResult
    save(outFile, 'OUT', '-v7.3');
    fprintf('Saved color x attention time-course result to: %s\n', outFile);
end

end

%% ============================= Local helpers =============================
function C = toCellStrLocal(x)
if iscell(x)
    C = x;
elseif isstring(x)
    C = cellstr(x);
elseif iscategorical(x)
    C = cellstr(x);
elseif ischar(x)
    C = cellstr(x);
else
    C = cellstr(string(x));
end
C = C(:);
for i = 1:numel(C)
    if iscell(C{i}) && numel(C{i}) == 1
        C{i} = C{i}{1};
    end
end
end

function [tc, nTrialsUsed, nStimUsed] = weightedTimeCourseLocal(respStimByTime, nTrialsStim, stimMask)
% respStimByTime: [nStim x nBins], nTrialsStim: [1 x nStim], stimMask: [1 x nStim]
mask = logical(stimMask(:))';
nTrials = double(nTrialsStim(:))';
mask = mask & isfinite(nTrials) & (nTrials > 0);
nStimUsed = nnz(mask);
nTrialsUsed = sum(nTrials(mask));

nBins = size(respStimByTime, 2);
tc = nan(1, nBins);
if nTrialsUsed <= 0
    return;
end

w = nTrials(mask);
X = respStimByTime(mask, :);

W = repmat(w(:), 1, nBins);
ok = isfinite(X);
num = sum((X .* W) .* ok, 1, 'omitnan');
den = sum(W .* ok, 1, 'omitnan');

use = den > 0;
tc(use) = num(use) ./ den(use);
end

function addSemPatchLocal(ax, t, m, s, col, fa)
t = t(:);
m = m(:);
s = s(:);
if numel(t) ~= numel(m) || numel(t) ~= numel(s)
    return;
end
ok = isfinite(t) & isfinite(m) & isfinite(s);
if nnz(ok) < 3
    return;
end
tt = t(ok);
lo = m(ok) - s(ok);
hi = m(ok) + s(ok);
patch(ax, [tt; flipud(tt)], [lo; flipud(hi)], col, ...
    'FaceAlpha', fa, 'EdgeColor', 'none');
end

function [p, testUsed] = oneSamplePLocal(x, alpha, prefTest)
x = x(:);
x = x(isfinite(x));
if numel(x) < 3
    p = NaN;
    testUsed = "n<3";
    return;
end

pref = lower(string(prefTest));
if (pref == "ttest")
    [~, p] = ttest(x, 0, 'Alpha', alpha);
    testUsed = "ttest";
    return;
end

if exist('signrank', 'file') == 2
    p = signrank(x, 0, 'alpha', alpha);
    testUsed = "signrank";
else
    [~, p] = ttest(x, 0, 'Alpha', alpha);
    testUsed = "ttest_fallback";
end
end

function [q, sig] = fdrBhLocal(p, alpha)
origSize = size(p);
p = double(p(:));
q = nan(size(p));
sig = false(size(p));

ok = isfinite(p);
pv = p(ok);
m = numel(pv);
if m == 0
    q = reshape(q, origSize);
    sig = reshape(sig, origSize);
    return;
end

[ps, ord] = sort(pv);
rk = (1:m)';
thr = (rk / m) * alpha;
k = find(ps <= thr, 1, 'last');
if ~isempty(k)
    sigSorted = false(m,1);
    sigSorted(1:k) = true;
    sigTmp = false(m,1);
    sigTmp(ord) = sigSorted;
    sig(ok) = sigTmp;
end

qSorted = nan(m,1);
qRaw = (m ./ rk) .* ps;
qRaw = min(qRaw, 1);
qSorted(end) = qRaw(end);
for i = m-1:-1:1
    qSorted(i) = min(qRaw(i), qSorted(i+1));
end
qTmp = nan(m,1);
qTmp(ord) = qSorted;
q(ok) = qTmp;

q = reshape(q, origSize);
sig = reshape(sig, origSize);
end
