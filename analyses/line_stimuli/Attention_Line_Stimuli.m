% Attention_Line_Stimuli
% Attention modulation summary, still-frame diagnostics, and optional movie wrapper.

%% Run switches
RUN_HISTOGRAM = false;
RUN_QC = false;
RUN_POST_AFFINE_VALUES = false;
RUN_PREP_ANCHOR_KNN = false;
RUN_PREP_NOISE_SIGNAL = false;
RUN_STILLS = false;
RUN_MOVIE = false;
RUN_GROUP_POLYGON_PREP = true;
RUN_GROUP_POLYGON_STILLS = false;
RUN_GROUP_POLYGON_MOVIE = false;
RUN_GROUP_POLYGON_TIMESERIES = false;
RUN_GROUP_POLYGON_SIGMOID_FIT = true;
RUN_GROUP_POLYGON_EXGAUSS_FIT = false;

%% Central parameters
P = struct();
P.timeIdx3 = 3;                      % 3-bin dataset: 300-500 ms
P.groupN = 8;                        % number of along_GC quantile groups
% P.respCapsulesFile = 'Resp_capsules_N_d12.mat'; % high-resolution response windows
% P.resultTag = 'd12';                 % suffix for output files to avoid overwriting runs
% P.targetWin = [300 310];             % target window request [ms] (nearest bin selected)
%P.targetWinPre = [-100 -90];         % pre-stim control window (nearest bin)
P.respCapsulesFile = 'Resp_capsules_N_d20.mat';
P.resultTag = 'd20';
P.targetWin = [300 320];
P.targetWinPre = [-100 -80];
P.targetWin10 = P.targetWin;         % backwards-compatible alias
P.pThresh = 0.05;                    % hard significance rule
P.siteInclusionMode = 'visual'; % 'significant' or 'visual'
P.visualSNRthr = 0.70;               % used when siteInclusionMode='visual'
P.siteWeightMode = 'dprime';    % 'reliability' or 'dprime'
P.siteRange = 1:512;
P.stimIDExample = 38;
P.idxClip = 2;                       % histogram clipping
P.plotAlpha = 0.3;                   % constant across windows
P.plotRobustPct = 95;
P.plotBgColor = [0.5 0.5 0.5];
P.plotCLow = [0.50 0.50 0.50];       % weak effects start at background tone
P.plotCHigh = [0.85 0.05 0.05];
P.neighborN = 1;                      % RF-center smoothing (set 1 to disable)
P.pixelNeighborN = 140;               % image-space KNN averaging after affine projectionNow N
P.pixelSmoothSigma = 0;               % optional Gaussian fallback in image space
P.alphaValueGamma = 4.0;
P.alphaValueMinScale = 0.01;         % keep weak points slightly visible
P.alphaFloorSoft = 0.20;             % soften transition around alpha floor
P.alphaFullAt = 0.20;                % normalized T-D value where alpha reaches 1
P.denMin = 0.03;                     % minimum |OUT3.muT-muD| for site normalization
P.preStimCalibratedAlpha = false;
P.preStimPercentile = 99.9;          % fewer pre-stim points near full opacity
P.plotMarkerSize = 8;
P.alphaByValue = true;
P.hotScale = true;                 % default: gray->red only (no yellow/white hot colors)
P.colorHotMaxFactor = 8.0;           % allow values > cMax to become hotter
P.colorRedAt = 0.20;                 % value where point color reaches red
P.hardAlphaCutoff = true;           % if true: alpha=0 below threshold, alpha=plotAlpha above
P.useSiteScale = true;               % per-site late-phase scaling
P.siteScaleLateWindow = [300 500];   % ms, used to build per-site scale
P.siteScaleStat = 'mean_plus_sd';    % per-site scale from late-phase values
P.siteScaleMin = 1e-6;
P.prepK = 500;                        % K for global post-affine KNN averaging
P.prepNPick = 5;                     % number of anchor combinations for diagnostics
P.prepKList = [10 50 80 200 500 1000];    % K sweep for noise/signal prep
P.preEndMs = 0;                      % pre bins: end <= this time
P.postStartMs = 300;                 % post bins: start >= this time
P.preQuantilePct = 99.98;               % threshold from pre-stim |value| quantile
P.alphaFloorPct = 80;                % suggested alpha floor from pre-stim
P.cMaxPostPct = 95;                  % suggested cMax from post-stim
P.groupPolygonShrink = 0.8;          % boundary() shrink factor for group polygons
P.groupPreQuantilePct = P.preQuantilePct; % grouped threshold from pre-stim quantile
P.groupCMaxPostPct = P.cMaxPostPct;  % grouped color max from post-stim quantile
P.knnUseSiteWeights = true;          % if true: fixed site-weighted KNN averaging
P.groupUseSiteWeights = true;        % if true: fixed site-weighted group mean
P.groupWeightLambda = 1e-6;          % stabilizer in weight denominator
P.groupWeightUseNmatch = true;       % in reliability mode, multiply by sqrt(Nmatch), Nmatch = wY + wP
P.groupWeightClipPct = 95;           % cap extreme site weights at this percentile
P.groupWeightMin = 0;                % optional floor on positive weights (0 = disabled)
P.groupTraceMode = 'overlay';        % 'overlay' (all curves in one axes) or 'subplots'
P.groupTraceSmoothW = 3;             % moving-average width (bins) for grouped traces
P.groupTraceHalfMaxFrac = 0.5;       % first-crossing level as fraction of smoothed max
P.groupTraceHalfSearchStartMs = 0;   % search start for half-max crossing
P.groupTraceShowHalfMarkers = true;  % draw circle markers at half-max crossing
P.groupTraceHalfMarkerSize = 7;      % marker size for half-max points
P.groupSigmoidFitStartMs = 0;        % shared-slope sigmoid fit start (ms)
P.groupSigmoidFitEndMs = 500;        % shared-slope sigmoid fit end (ms)
P.groupSigmoidSmoothW = 3;           % smoothing width (bins) before sigmoid fit
P.groupSigmoidUseAbs = false;        % fit abs(trace) if true
P.groupSigmoidMinTau = 5;            % lower bound on shared slope parameter tau
P.groupSigmoidMaxTau = 400;          % upper bound on shared slope parameter tau
P.groupSigmoidInitTau = 40;          % initial tau (ms)
P.groupSigmoidShowT50 = true;        % show model t50 marker per group
P.groupSigmoidT50MarkerSize = 7;     % marker size for model t50 points
P.groupSigmoidRegressionX = 'alongMid'; % regression x-axis: 'alongMid'|'groupIdx'|'nComb'
P.groupSigmoidExcludeHighestXRegression = false; % exclude highest-x point from t50 regression
P.groupExGaussFitStartMs = -200;     % fit start (ms)
P.groupExGaussFitEndMs = 500;        % fit end (ms)
P.groupExGaussSmoothW = 1;           % moving-average width before fit (bins)
P.groupExGaussMinSigma = 5;          % lower bound for sigma
P.groupExGaussMaxSigma = 400;        % upper bound for sigma
P.groupSigMinRunBins = 5;            % onset criterion: >= this many consecutive significant bins
P.groupSigSearchStartMs = 0;         % onset search starts at this time (ms)
P.stillRunConsistencyCheck = false;  % expensive global KNN check; disable for speed
P.timeLabelRef = 'start';            % 'start' keeps stimulus onset aligned at 0 ms
P.stimOnsetMs = 0;                   % stimulus becomes visible from this bin start (ms)

% Build a safe file suffix for run-specific outputs (e.g., _d12, _d20).
runTag = char(string(P.resultTag));
if isempty(strtrim(runTag))
    [~, runTag] = fileparts(P.respCapsulesFile);
end
runTag = regexprep(runTag, '[^A-Za-z0-9_-]', '_');
tagSuffix = ['_' runTag];

siteInclusionMode = lower(string(P.siteInclusionMode));
siteWeightMode = lower(string(P.siteWeightMode));
assert(any(siteInclusionMode == ["significant","visual"]), ...
    'P.siteInclusionMode must be ''significant'' or ''visual''.');
assert(any(siteWeightMode == ["reliability","dprime"]), ...
    'P.siteWeightMode must be ''reliability'' or ''dprime''.');

selectionTagSuffix = '';
if siteInclusionMode ~= "significant" || siteWeightMode ~= "reliability"
    modeParts = strings(0,1);
    if siteInclusionMode == "visual"
        snrTag = strrep(sprintf('%.2f', P.visualSNRthr), '.', 'p');
        modeParts(end+1) = "allvis"; %#ok<SAGROW>
        modeParts(end+1) = "snr" + string(snrTag); %#ok<SAGROW>
    end
    if siteWeightMode == "dprime"
        modeParts(end+1) = "wdprime"; %#ok<SAGROW>
    elseif siteWeightMode == "reliability" && siteInclusionMode ~= "significant"
        modeParts(end+1) = "wreliability"; %#ok<SAGROW>
    end
    selectionTagSuffix = ['_' char(strjoin(modeParts, '_'))];
end


%% Required context
cfg = config();
TabFile = 'ObjAtt_lines_monkeyN_20220201_B1';

Sgeo = load(fullfile(cfg.logsDir, TabFile));
assert(isfield(Sgeo, 'ALLCOORDS'), ...
    '%s must contain ALLCOORDS.', fullfile(cfg.logsDir, TabFile));
ALLCOORDS = Sgeo.ALLCOORDS;

Srtab = load(fullfile(cfg.logsDir, 'RTAB384.mat'));
assert(isfield(Srtab, 'RTAB384'), ...
    '%s must contain RTAB384.', fullfile(cfg.logsDir, 'RTAB384.mat'));
RTAB384 = Srtab.RTAB384;

S = load(fullfile(cfg.matDir, 'Tall_V1_lines_N.mat'));
assert(isfield(S, 'Tall_V1') && isstruct(S.Tall_V1), ...
    'Tall_V1_lines_N.mat must contain struct Tall_V1.');
Tall_V1 = S.Tall_V1;

% Avoid cross-script contamination when switching between V1 and V4 runs.
clear R_resp OUT_postAffine GGROUP FITS FITG OUTTarget OUTPre

%% Load 3-bin response + normalization and get baseline OUT3
S = load(fullfile(cfg.matDir, 'SNR_capsules_N_d12.mat'));   % loads R
assert(isfield(S, 'R') && isstruct(S.R), 'SNR_capsules_N_d12.mat must contain struct R.');
R3 = S.R;
assert(isfield(R3, 'timeWindows') && size(R3.timeWindows,1) == 3, ...
    'Expected 3 time windows in SNR_capsules_N_d12.mat.');

S = load(fullfile(cfg.matDir, 'SNR_V1_byColor_byWindow.mat')); % loads SNR
assert(isfield(S, 'SNR') && isstruct(S.SNR) && isfield(S.SNR, 'muSpont'), ...
    'SNR_V1_byColor_byWindow.mat must contain normalization struct SNR.');
SNRnorm = S.SNR;

saveTag = 'OUT_attention_modulation_3bin_timeIdx3';
outFile = fullfile(cfg.resultsDir, [saveTag '.mat']);

if exist(outFile, 'file')
    S = load(outFile, 'OUT');
    OUT3 = S.OUT;
    fprintf('Loaded baseline OUT from: %s\n', outFile);
else
    opts3 = struct('timeIdx', P.timeIdx3, 'excludeOverlap', true, 'verbose', true);
    OUT3 = attention_modulation_V1_3bin(R3, Tall_V1, SNRnorm, opts3);
    meta = struct();
    meta.created = datestr(now, 30);
    meta.script = mfilename;
    meta.timeIdx = opts3.timeIdx;
    meta.excludeOverlap = opts3.excludeOverlap;
    meta.note = 'Use OUT.pValueTD from this file as fixed significance mask.';
    OUT = OUT3; %#ok<NASGU>
    save(outFile, 'OUT', 'meta', '-v7.3');
    fprintf('Computed and saved baseline OUT to: %s\n', outFile);
end

fprintf('Median d'': %.3f\n', median(OUT3.dprime, 'omitnan'));
fprintf('Median index: %.3f\n', median(OUT3.index, 'omitnan'));

snrFields = {'yellowEarly','yellowLate','purpleEarly','purpleLate'};
nSitesAll = numel(OUT3.pValueTD);
SNRmatAll = nan(nSitesAll, numel(snrFields));
for k = 1:numel(snrFields)
    assert(isfield(SNRnorm, snrFields{k}), ...
        'SNR normalization struct missing field %s.', snrFields{k});
    v = double(SNRnorm.(snrFields{k})(:));
    assert(numel(v) >= nSitesAll, ...
        'SNR normalization field %s has only %d values; need at least %d.', ...
        snrFields{k}, numel(v), nSitesAll);
    SNRmatAll(:,k) = v(1:nSitesAll);
end
bestSNRAll = max(SNRmatAll, [], 2, 'omitnan');

% Build shared site-selection mask and optional fixed per-site weights.
sigSiteByIndexAll = isfinite(OUT3.pValueTD) & (OUT3.pValueTD < P.pThresh);
visSiteByIndexAll = isfinite(bestSNRAll) & (bestSNRAll > P.visualSNRthr);
selectedSiteByIndexAll = sigSiteByIndexAll;
selectedSiteLabel = sprintf('significant sites (pTD < %.3f)', P.pThresh);
switch char(siteInclusionMode)
    case 'significant'
        selectedSiteByIndexAll = sigSiteByIndexAll;
    case 'visual'
        selectedSiteByIndexAll = visSiteByIndexAll;
        selectedSiteLabel = sprintf('visually responsive sites (bestSNR > %.2f)', P.visualSNRthr);
end
fprintf('Site inclusion mode: %s | kept %d / %d sites\n', ...
    selectedSiteLabel, nnz(selectedSiteByIndexAll), numel(selectedSiteByIndexAll));
if siteWeightMode == "dprime"
    fprintf('Site weight mode: pure |d''| (no sqrt(Nmatch) multiplier)\n');
else
    fprintf('Site weight mode: %s\n', char(siteWeightMode));
end

siteWeightsByIndexAll = [];
if P.knnUseSiteWeights || P.groupUseSiteWeights
    switch char(siteWeightMode)
        case 'reliability'
            assert(all(isfield(OUT3, {'muT','muD','varT','varD'})), ...
                'OUT3 must contain muT/muD/varT/varD for reliability weighting.');
            dSite = double(OUT3.muT(:) - OUT3.muD(:));
            varSite = 0.5 * (double(OUT3.varT(:)) + double(OUT3.varD(:)));
            varSite(~isfinite(varSite) | varSite < 0) = 0;
            den = sqrt(varSite + max(P.groupWeightLambda, 0));
            den(~isfinite(den) | den <= 0) = NaN;
            siteWeightsByIndexAll = abs(dSite) ./ den;
        case 'dprime'
            assert(isfield(OUT3, 'dprime'), ...
                'OUT3 must contain dprime for dprime weighting.');
            siteWeightsByIndexAll = abs(double(OUT3.dprime(:)));
    end

    applyNmatchWeight = P.groupWeightUseNmatch;
    if siteWeightMode == "dprime"
        applyNmatchWeight = false;
    end
    if applyNmatchWeight
        assert(all(isfield(OUT3, {'wY','wP'})), ...
            'OUT3 must contain wY/wP when P.groupWeightUseNmatch=true.');
        nMatch = double(OUT3.wY(:) + OUT3.wP(:));
        nMatch(~isfinite(nMatch) | nMatch < 0) = 0;
        siteWeightsByIndexAll = siteWeightsByIndexAll .* sqrt(nMatch);
    end

    siteWeightsByIndexAll(~isfinite(siteWeightsByIndexAll) | siteWeightsByIndexAll < 0) = 0;
    siteWeightsByIndexAll(~selectedSiteByIndexAll) = 0;

    wPos = siteWeightsByIndexAll(siteWeightsByIndexAll > 0);
    if ~isempty(wPos) && isfinite(P.groupWeightClipPct) && P.groupWeightClipPct > 0 && P.groupWeightClipPct < 100
        wCap = prctile(wPos, P.groupWeightClipPct);
        if isfinite(wCap) && wCap > 0
            siteWeightsByIndexAll = min(siteWeightsByIndexAll, wCap);
        end
    end
    if isfinite(P.groupWeightMin) && P.groupWeightMin > 0
        pos = siteWeightsByIndexAll > 0;
        siteWeightsByIndexAll(pos) = max(siteWeightsByIndexAll(pos), P.groupWeightMin);
    end

    wPos = siteWeightsByIndexAll(siteWeightsByIndexAll > 0);
    if isempty(wPos)
        warning('Site weighting requested but no positive site weights were obtained. Falling back to unweighted averaging.');
        siteWeightsByIndexAll = [];
    else
        fprintf(['Fixed site weights (%s) ready: nPos=%d | min/median/max=%.6g / %.6g / %.6g | clip p%.2f\n'], ...
            char(siteWeightMode), numel(wPos), min(wPos), median(wPos), max(wPos), P.groupWeightClipPct);
    end
end

%% Histogram: attention index distribution
if RUN_HISTOGRAM
    idxAll = OUT3.index;
    isSig = isfinite(OUT3.pValueTD) & (OUT3.pValueTD < P.pThresh);

    vAllRaw = idxAll(isfinite(idxAll));
    vSigRaw = idxAll(isSig & isfinite(idxAll));
    vAll = max(min(vAllRaw, P.idxClip), -P.idxClip);
    vSig = max(min(vSigRaw, P.idxClip), -P.idxClip);

    fprintf('Clipped attention index for display at +/-%.1f (all: %d, sig: %d values clipped)\n', ...
        P.idxClip, nnz(abs(vAllRaw) > P.idxClip), nnz(abs(vSigRaw) > P.idxClip));

    figure('Color','w');
    hold on;
    histogram(vAll, 30, 'FaceColor', [0.80 0.80 0.80], 'EdgeColor', 'none');
    histogram(vSig, 30, 'FaceColor', [0.85 0.20 0.20], 'EdgeColor', 'none');
    xlabel('Attention index');
    ylabel('Number of sites');
    title(sprintf('Attention index (all vs significant, pTD < %.3f)', P.pThresh));
    legend(sprintf('All sites (N=%d)', numel(vAll)), ...
           sprintf('Significant sites (N=%d)', numel(vSig)), ...
           'Location','best');
    grid on;
end

%% Load high-resolution response windows
if RUN_QC || RUN_POST_AFFINE_VALUES || RUN_STILLS || RUN_MOVIE || RUN_GROUP_POLYGON_PREP || RUN_GROUP_POLYGON_STILLS || RUN_GROUP_POLYGON_MOVIE
    S = load(fullfile(cfg.matDir, P.respCapsulesFile));  % loads R
    assert(isfield(S, 'R') && isstruct(S.R), ...
        '%s must contain struct R.', P.respCapsulesFile);
    R_resp = S.R;

    assert(isfield(R_resp, 'timeWindows') && size(R_resp.timeWindows,2) == 2, ...
        '%s must contain R.timeWindows as [nWindows x 2].', P.respCapsulesFile);
    assert(size(R_resp.timeWindows,1) > 3, ...
        ['%s appears to have only %d windows. ' ...
         'Expected many short windows for still/movie plotting.'], P.respCapsulesFile, size(R_resp.timeWindows,1));

    targetWinReq = P.targetWin10;
    if isfield(P, 'targetWin') && isnumeric(P.targetWin) && numel(P.targetWin) == 2
        targetWinReq = P.targetWin;
    end
    [~, timeIdxTarget] = min(sum(abs(R_resp.timeWindows - targetWinReq),2));
    winTarget = R_resp.timeWindows(timeIdxTarget,:);
    assert((winTarget(2)-winTarget(1)) <= 25, ...
        'Selected target bin is [%d %d] ms (duration %.1f ms).', winTarget(1), winTarget(2), winTarget(2)-winTarget(1));

    [~, timeIdxPre] = min(sum(abs(R_resp.timeWindows - P.targetWinPre),2));
    winPre = R_resp.timeWindows(timeIdxPre,:);
    assert((winPre(2)-winPre(1)) <= 25, ...
        'Selected pre-stim bin is [%d %d] ms (duration %.1f ms).', winPre(1), winPre(2), winPre(2)-winPre(1));

    optsTarget = struct('timeIdx',timeIdxTarget,'excludeOverlap',true,'verbose',false);
    optsPre = struct('timeIdx',timeIdxPre,'excludeOverlap',true,'verbose',false);
    OUTTarget = attention_modulation_V1_3bin(R_resp, Tall_V1, SNRnorm, optsTarget);
    OUTPre = attention_modulation_V1_3bin(R_resp, Tall_V1, SNRnorm, optsPre);

end

%% QC diagnostics
if RUN_QC
    d3 = OUT3.muT - OUT3.muD;
    dTarget = OUTTarget.muT - OUTTarget.muD;
    dPre = OUTPre.muT - OUTPre.muD;

    okTarget = isfinite(d3) & isfinite(dTarget);
    if any(okTarget)
        rTarget = corr(d3(okTarget), dTarget(okTarget));
        ddTarget = abs(d3(okTarget) - dTarget(okTarget));
        fprintf('Delta(T-D) 3-bin vs target bin: corr=%.4f, median|diff|=%.6g, max|diff|=%.6g\n', ...
            rTarget, median(ddTarget), max(ddTarget));
    else
        fprintf('Delta(T-D) 3-bin vs target bin unavailable (no finite overlap).\n');
    end

    okPre = isfinite(d3) & isfinite(dPre);
    if any(okPre)
        rPre = corr(d3(okPre), dPre(okPre));
        ddPre = abs(d3(okPre) - dPre(okPre));
        fprintf('Delta(T-D) 3-bin vs pre-stim: corr=%.4f, median|diff|=%.6g, max|diff|=%.6g\n', ...
            rPre, median(ddPre), max(ddPre));
    end

    fprintf('nSig pre-stim (pTD < %.2f): %d / %d sites\n', P.pThresh, ...
        nnz(isfinite(OUTPre.pValueTD) & OUTPre.pValueTD < P.pThresh), numel(OUTPre.pValueTD));

    sigPre = isfinite(OUT3.pValueTD) & (OUT3.pValueTD < P.pThresh) & isfinite(dPre);
    nPrePos = nnz(sigPre & (dPre > 0));   % T > D
    nPreNeg = nnz(sigPre & (dPre < 0));   % D > T
    nPreTot = nPrePos + nPreNeg;
    if nPreTot > 0
        fracPos = nPrePos / nPreTot;
        fracNeg = nPreNeg / nPreTot;
        fprintf(['Pre-stim sign balance (sig sites): T>D: %d, D>T: %d, ' ...
                 'frac(T>D)=%.3f, frac(D>T)=%.3f\n'], nPrePos, nPreNeg, fracPos, fracNeg);
        if abs(fracPos - 0.5) > 0.15
            warning(['Pre-stim sign balance is noticeably asymmetric. ' ...
                     'Check normalization/trial selection if this persists.']);
        end
    else
        fprintf('Pre-stim sign balance (sig sites): no non-zero d(T-D) values.\n');
    end
end

%% Optional export of post-affine signed T-D values for all bins
outValuesFile = fullfile(cfg.resultsDir, sprintf('post_affine_delta_points_allbins_stim%d%s%s.mat', ...
    P.stimIDExample, tagSuffix, selectionTagSuffix));
groupedFile = fullfile(cfg.resultsDir, sprintf('grouped_alonggc_polygons_stim%d_N%d%s.mat', ...
    P.stimIDExample, P.groupN, [tagSuffix selectionTagSuffix]));
if RUN_POST_AFFINE_VALUES
    siteRangeVals = P.siteRange(:)';
    selectedMaskVals = selectedSiteByIndexAll(siteRangeVals);

    OUT_postAffine = compute_projected_delta_points_allbins( ...
        P.stimIDExample, Tall_V1, ALLCOORDS, RTAB384, R_resp, SNRnorm, ...
        'siteRange', P.siteRange, ...
        'excludeOverlap', true, ...
        'stimIdx', 1:384, ...
        'sigSiteMask', selectedMaskVals, ...
        'saveFile', outValuesFile, ...
        'verbose', true);
    fprintf('Post-affine values saved (%s): %s\n', selectedSiteLabel, outValuesFile);
end

%% Optional pre-movie anchor/KNN diagnostics on post-affine values
if RUN_PREP_ANCHOR_KNN
    if ~exist('OUT_postAffine', 'var')
        assert(exist(outValuesFile, 'file') == 2, ...
            'Post-affine values file missing: %s', outValuesFile);
        S = load(outValuesFile);
        if isfield(S, 'OUT')
            OUT_postAffine = S.OUT;
        elseif isfield(S, 'D')
            OUT_postAffine = S.D;
        else
            error('File %s must contain OUT or D.', outValuesFile);
        end
    end

    prepFile = fullfile(cfg.resultsDir, sprintf('anchor_knn_prep_stim%d%s%s.mat', ...
        P.stimIDExample, tagSuffix, selectionTagSuffix));
    PREP = analyze_anchor_knn_timeseries( ...
        OUT_postAffine, Tall_V1, OUT3, P.stimIDExample, ...
        'siteRange', P.siteRange, ...
        'siteMask', selectedSiteByIndexAll, ...
        'pThresh', P.pThresh, ...
        'K', P.prepK, ...
        'nPick', P.prepNPick, ...
        'makePlot', true, ...
        'verbose', true, ...
        'saveFile', prepFile);
    fprintf('Saved pre-movie anchor/KNN diagnostics to: %s\n', prepFile);
end

%% Optional pre-movie noise vs signal calibration across K
if RUN_PREP_NOISE_SIGNAL
    if ~exist('OUT_postAffine', 'var')
        assert(exist(outValuesFile, 'file') == 2, ...
            'Post-affine values file missing: %s', outValuesFile);
        S = load(outValuesFile);
        if isfield(S, 'OUT')
            OUT_postAffine = S.OUT;
        elseif isfield(S, 'D')
            OUT_postAffine = S.D;
        else
            error('File %s must contain OUT or D.', outValuesFile);
        end
    end

    prepNoiseFile = fullfile(cfg.resultsDir, sprintf('knn_noise_signal_prep_stim%d%s%s.mat', ...
        P.stimIDExample, tagSuffix, selectionTagSuffix));
    siteWeightsKNN = [];
    if P.knnUseSiteWeights
        siteWeightsKNN = siteWeightsByIndexAll;
    end
    PREP_NOISE = analyze_knn_noise_signal_thresholds( ...
        OUT_postAffine, Tall_V1, ...
        'KList', P.prepKList, ...
        'preEndMs', P.preEndMs, ...
        'postStartMs', P.postStartMs, ...
        'preQuantilePct', P.preQuantilePct, ...
        'alphaFloorPct', P.alphaFloorPct, ...
        'cMaxPostPct', P.cMaxPostPct, ...
        'kRef', P.prepK, ...
        'enforceK', true, ...
        'siteWeightsByIndex', siteWeightsKNN, ...
        'makePlot', true, ...
        'verbose', true, ...
        'saveFile', prepNoiseFile);
    fprintf('Saved pre-movie noise/signal calibration to: %s\n', prepNoiseFile);
end

%% Optional grouped along_GC polygon prep (equal-count groups over selected combos)
if RUN_GROUP_POLYGON_PREP
    if ~exist('OUT_postAffine', 'var')
        assert(exist(outValuesFile, 'file') == 2, ...
            'Post-affine values file missing: %s', outValuesFile);
        S = load(outValuesFile);
        if isfield(S, 'OUT')
            OUT_postAffine = S.OUT;
        elseif isfield(S, 'D')
            OUT_postAffine = S.D;
        else
            error('File %s must contain OUT or D.', outValuesFile);
        end
    end

    sigSiteByIndex = selectedSiteByIndexAll;
    siteWeightsByIndex = [];
    if P.groupUseSiteWeights
        siteWeightsByIndex = siteWeightsByIndexAll;
    end

    GGROUP = build_grouped_alonggc_polygons_allbins( ...
        OUT_postAffine, Tall_V1, ...
        'nGroups', P.groupN, ...
        'sigSiteByIndex', sigSiteByIndex, ...
        'siteWeightsByIndex', siteWeightsByIndex, ...
        'preEndMs', P.preEndMs, ...
        'postStartMs', P.postStartMs, ...
        'preQuantilePct', P.groupPreQuantilePct, ...
        'cMaxPostPct', P.groupCMaxPostPct, ...
        'polygonShrink', P.groupPolygonShrink, ...
        'saveFile', groupedFile, ...
        'verbose', true);

    fprintf('Saved grouped along_GC polygons to: %s\n', groupedFile);
    if isfield(GGROUP, 'groupSummary')
        disp(GGROUP.groupSummary(:, {'groupIdx','nComb','alongMin','alongMax'}));
    end
end

%% Still plots
if RUN_STILLS

    % Use post-affine points + prep thresholds (same value pipeline as prep).
    if ~exist('OUT_postAffine', 'var')
        assert(exist(outValuesFile, 'file') == 2, ...
            'Post-affine values file missing: %s', outValuesFile);
        S = load(outValuesFile);
        if isfield(S, 'OUT')
            OUT_postAffine = S.OUT;
        elseif isfield(S, 'D')
            OUT_postAffine = S.D;
        else
            error('File %s must contain OUT or D.', outValuesFile);
        end
    end
    assert(isfield(OUT_postAffine, 'bins') && numel(OUT_postAffine.bins) >= max([timeIdxPre timeIdxTarget]), ...
        'OUT_postAffine.bins missing or shorter than requested time indices.');
    assert(isfield(OUT_postAffine.bins, 'stream'), ...
        ['OUT_postAffine.bins.stream missing. Re-run RUN_POST_AFFINE_VALUES with the ' ...
         'updated compute_projected_delta_points_allbins.m']);

    prepNoiseFile = fullfile(cfg.resultsDir, sprintf('knn_noise_signal_prep_stim%d%s%s.mat', ...
        P.stimIDExample, tagSuffix, selectionTagSuffix));
    assert(exist(prepNoiseFile, 'file') == 2, ...
        'Noise/signal prep file missing: %s. Run RUN_PREP_NOISE_SIGNAL first.', prepNoiseFile);
    Sns = load(prepNoiseFile);
    assert(isfield(Sns, 'R') && isfield(Sns.R, 'summary') && istable(Sns.R.summary), ...
        'Prep file %s must contain R.summary table.', prepNoiseFile);
    Tns = Sns.R.summary;

    rowK = find(Tns.K == P.prepK, 1, 'first');
    assert(~isempty(rowK), ...
        'Requested prep K=%d not found in summary. Update P.prepK or rerun prep with this K.', P.prepK);
    Kuse = double(Tns.K(rowK));
    thrUse = double(Tns.thresholdPreQ(rowK));
    assert(isfinite(thrUse) && thrUse > 0, ...
        'Invalid thresholdPreQ at K=%d in prep summary.', Kuse);
    cMaxUse = double(Tns.cMaxSuggest(rowK));
    if ~isfinite(cMaxUse) || cMaxUse <= 0
        cMaxUse = [];
    end
    fprintf('RUN_STILLS using prep settings: K=%d, threshold=%.6g, cMax=%s\n', ...
        Kuse, thrUse, mat2str(cMaxUse));
    if ismember('preExceedFrac', Tns.Properties.VariableNames) && ismember('postExceedFrac', Tns.Properties.VariableNames)
        fprintf('Prep expected exceedance at K=%d: pre=%.2f%%, post=%.2f%%\n', ...
            Kuse, 100*double(Tns.preExceedFrac(rowK)), 100*double(Tns.postExceedFrac(rowK)));
    end
    siteWeightsKNN = [];
    if P.knnUseSiteWeights
        siteWeightsKNN = siteWeightsByIndexAll;
    end

    if P.stillRunConsistencyCheck
        % Consistency check: recompute exceedance directly from OUT_postAffine.
        preMaskBins = R_resp.timeWindows(:,2) <= P.preEndMs;
        postMaskBins = R_resp.timeWindows(:,1) >= P.postStartMs;
        Scheck = summarize_post_affine_threshold_exceedance( ...
            OUT_postAffine, Kuse, thrUse, preMaskBins, postMaskBins, false, siteWeightsKNN);
        fprintf('Consistency check from OUT_postAffine: pre=%.2f%%, post=%.2f%%\n', ...
            100*Scheck.preFrac, 100*Scheck.postFrac);
        if ismember('preExceedFrac', Tns.Properties.VariableNames) && ismember('postExceedFrac', Tns.Properties.VariableNames)
            dPre = abs(Scheck.preFrac - double(Tns.preExceedFrac(rowK)));
            dPost = abs(Scheck.postFrac - double(Tns.postExceedFrac(rowK)));
            fprintf('Prep vs check difference: pre=%.2f%%, post=%.2f%%\n', 100*dPre, 100*dPost);
            if dPre > 0.01 || dPost > 0.01
                warning(['Prep and direct OUT_postAffine check disagree by >1%%. ' ...
                         'Check if prep/out files are stale or generated with different settings.']);
            end
        end
    else
        Scheck = struct('fracByBin', nan(numel(OUT_postAffine.bins),1));
    end

    showStimTarget = (winTarget(1) >= P.stimOnsetMs);
    showStimPre = (winPre(1) >= P.stimOnsetMs);

    % ===== POST-STIM FRAME =====
    hSmall = plot_post_affine_knn_frame( ...
        P.stimIDExample, ALLCOORDS, RTAB384, OUT_postAffine.bins(timeIdxTarget), ...
        'K', Kuse, ...
        'siteWeightsByIndex', siteWeightsKNN, ...
        'alphaFullAt', thrUse, ...
        'colorRedAt', thrUse, ...
        'cMaxFixed', cMaxUse, ...
        'markerSize', P.plotMarkerSize, ...
        'alpha', 1, ...
        'bgColor', P.plotBgColor, ...
        'cLow', P.plotCLow, ...
        'cHigh', P.plotCHigh, ...
        'hotScale', P.hotScale, ...
        'colorHotMaxFactor', P.colorHotMaxFactor, ...
        'hardAlphaCutoff', P.hardAlphaCutoff, ...
        'timeWindow', winTarget, ...
        'timeLabelRef', P.timeLabelRef, ...
        'showStimulus', showStimTarget, ...
        'enforceK', false);
    figure(hSmall.fig);
    set(hSmall.fig, 'Name', 'Attention TD | Small window (~10 ms)', 'NumberTitle', 'off');
    title(hSmall.ax, sprintf('Small-window attention map | [%d %d] ms | stim %d', ...
        winTarget(1), winTarget(2), P.stimIDExample), 'Color','w');
    fprintf('Using target bin %d: [%d %d] ms | plotted >thr: %.2f%%\n', ...
        timeIdxTarget, winTarget(1), winTarget(2), 100*hSmall.fracAboveThreshold);
    if isfinite(Scheck.fracByBin(timeIdxTarget))
        fprintf('  direct bin check >thr: %.2f%%\n', 100*Scheck.fracByBin(timeIdxTarget));
    end

    % ===== PRE-STIM FRAME =====
    hPre = plot_post_affine_knn_frame( ...
        P.stimIDExample, ALLCOORDS, RTAB384, OUT_postAffine.bins(timeIdxPre), ...
        'K', Kuse, ...
        'siteWeightsByIndex', siteWeightsKNN, ...
        'alphaFullAt', thrUse, ...
        'colorRedAt', thrUse, ...
        'cMaxFixed', cMaxUse, ...
        'markerSize', P.plotMarkerSize, ...
        'alpha', 1, ...
        'bgColor', P.plotBgColor, ...
        'cLow', P.plotCLow, ...
        'cHigh', P.plotCHigh, ...
        'hotScale', P.hotScale, ...
        'colorHotMaxFactor', P.colorHotMaxFactor, ...
        'hardAlphaCutoff', P.hardAlphaCutoff, ...
        'timeWindow', winPre, ...
        'timeLabelRef', P.timeLabelRef, ...
        'showStimulus', showStimPre, ...
        'enforceK', false);
    figure(hPre.fig);
    set(hPre.fig, 'Name', 'Attention TD | Pre-stim control', 'NumberTitle', 'off');
    title(hPre.ax, sprintf('Pre-stim control map | [%d %d] ms | stim %d', ...
        winPre(1), winPre(2), P.stimIDExample), 'Color','w');
    fprintf('Using pre-stim bin %d: [%d %d] ms | plotted >thr: %.2f%%\n', ...
        timeIdxPre, winPre(1), winPre(2), 100*hPre.fracAboveThreshold);
    if isfinite(Scheck.fracByBin(timeIdxPre))
        fprintf('  direct bin check >thr: %.2f%%\n', 100*Scheck.fracByBin(timeIdxPre));
    end
end

%% Grouped polygon stills (pre/post)
if RUN_GROUP_POLYGON_STILLS
    if ~exist('GGROUP', 'var')
        assert(exist(groupedFile, 'file') == 2, ...
            'Grouped polygon file missing: %s. Run RUN_GROUP_POLYGON_PREP first.', groupedFile);
        Sg = load(groupedFile);
        assert(isfield(Sg, 'G') && isstruct(Sg.G), ...
            'Grouped file %s must contain struct G.', groupedFile);
        GGROUP = Sg.G;
    end

    assert(isfield(GGROUP, 'groupMeanSigned') && isfield(GGROUP, 'timeWindows') && isfield(GGROUP, 'groupPolygons'), ...
        'Grouped struct missing required fields.');

    twGrp = double(GGROUP.timeWindows);
    [~, timeIdxTargetG] = min(sum(abs(twGrp - winTarget),2));
    [~, timeIdxPreG] = min(sum(abs(twGrp - winPre),2));
    winTargetG = twGrp(timeIdxTargetG,:);
    winPreG = twGrp(timeIdxPreG,:);

    thrGrp = [];
    cMaxGrp = [];
    if isfield(GGROUP, 'calibration') && isstruct(GGROUP.calibration)
        if isfield(GGROUP.calibration, 'thresholdPreQ')
            thrGrp = double(GGROUP.calibration.thresholdPreQ);
        end
        if isfield(GGROUP.calibration, 'cMaxSuggest')
            cMaxGrp = double(GGROUP.calibration.cMaxSuggest);
        end
    end
    if ~isfinite(thrGrp) || thrGrp <= 0
        thrGrp = 0.2;
    end
    if ~isfinite(cMaxGrp) || cMaxGrp <= 0
        cMaxGrp = [];
    end

    showStimTarget = (winTargetG(1) >= P.stimOnsetMs);
    showStimPre = (winPreG(1) >= P.stimOnsetMs);

    hGroupTarget = plot_grouped_alonggc_polygon_frame( ...
        P.stimIDExample, ALLCOORDS, RTAB384, GGROUP, timeIdxTargetG, ...
        'alphaFullAt', thrGrp, ...
        'colorRedAt', thrGrp, ...
        'cMaxFixed', cMaxGrp, ...
        'alpha', 1, ...
        'bgColor', P.plotBgColor, ...
        'cLow', P.plotCLow, ...
        'cHigh', P.plotCHigh, ...
        'hotScale', P.hotScale, ...
        'colorHotMaxFactor', P.colorHotMaxFactor, ...
        'hardAlphaCutoff', P.hardAlphaCutoff, ...
        'timeLabelRef', P.timeLabelRef, ...
        'showStimulus', showStimTarget);
    figure(hGroupTarget.fig);
    set(hGroupTarget.fig, 'Name', sprintf('Grouped polygons N=%d | post-stim', P.groupN), 'NumberTitle', 'off');
    title(hGroupTarget.ax, sprintf('Grouped polygons (N=%d) | [%d %d] ms | stim %d', ...
        P.groupN, winTargetG(1), winTargetG(2), P.stimIDExample), 'Color', 'w');

    hGroupPre = plot_grouped_alonggc_polygon_frame( ...
        P.stimIDExample, ALLCOORDS, RTAB384, GGROUP, timeIdxPreG, ...
        'alphaFullAt', thrGrp, ...
        'colorRedAt', thrGrp, ...
        'cMaxFixed', cMaxGrp, ...
        'alpha', 1, ...
        'bgColor', P.plotBgColor, ...
        'cLow', P.plotCLow, ...
        'cHigh', P.plotCHigh, ...
        'hotScale', P.hotScale, ...
        'colorHotMaxFactor', P.colorHotMaxFactor, ...
        'hardAlphaCutoff', P.hardAlphaCutoff, ...
        'timeLabelRef', P.timeLabelRef, ...
        'showStimulus', showStimPre);
    figure(hGroupPre.fig);
    set(hGroupPre.fig, 'Name', sprintf('Grouped polygons N=%d | pre-stim', P.groupN), 'NumberTitle', 'off');
    title(hGroupPre.ax, sprintf('Grouped polygons (N=%d) | [%d %d] ms | stim %d', ...
        P.groupN, winPreG(1), winPreG(2), P.stimIDExample), 'Color', 'w');

    fprintf(['Grouped stills (N=%d): pre >thr %.2f%%, post >thr %.2f%% | ' ...
             'pre routes T:%d D:%d, post routes T:%d D:%d\n'], ...
        P.groupN, 100*hGroupPre.fracAboveThreshold, 100*hGroupTarget.fracAboveThreshold, ...
        hGroupPre.nTargetGroups, hGroupPre.nDistrGroups, ...
        hGroupTarget.nTargetGroups, hGroupTarget.nDistrGroups);
end

%% Grouped polygon movie rendering
if RUN_GROUP_POLYGON_MOVIE
    if ~exist('GGROUP', 'var')
        assert(exist(groupedFile, 'file') == 2, ...
            'Grouped polygon file missing: %s. Run RUN_GROUP_POLYGON_PREP first.', groupedFile);
        Sg = load(groupedFile);
        assert(isfield(Sg, 'G') && isstruct(Sg.G), ...
            'Grouped file %s must contain struct G.', groupedFile);
        GGROUP = Sg.G;
    end

    thrGrp = [];
    cMaxGrp = [];
    if isfield(GGROUP, 'calibration') && isstruct(GGROUP.calibration)
        if isfield(GGROUP.calibration, 'thresholdPreQ')
            thrGrp = double(GGROUP.calibration.thresholdPreQ);
        end
        if isfield(GGROUP.calibration, 'cMaxSuggest')
            cMaxGrp = double(GGROUP.calibration.cMaxSuggest);
        end
    end
    if ~isfinite(thrGrp) || thrGrp <= 0
        thrGrp = 0.2;
    end
    if ~isfinite(cMaxGrp) || cMaxGrp <= 0
        cMaxGrp = [];
    end

    outMovieGrouped = fullfile(cfg.resultsDir, ...
        sprintf('V1_attentiondiff_groupedpolygons_N%d%s%s.mp4', ...
        P.groupN, tagSuffix, selectionTagSuffix));

    MOV_GROUP = make_grouped_alonggc_polygon_movie( ...
        outMovieGrouped, P.stimIDExample, ALLCOORDS, RTAB384, GGROUP, ...
        'alphaFullAt', thrGrp, ...
        'colorRedAt', thrGrp, ...
        'cMaxFixed', cMaxGrp, ...
        'alpha', P.plotAlpha, ...
        'bgColor', P.plotBgColor, ...
        'cLow', P.plotCLow, ...
        'cHigh', P.plotCHigh, ...
        'hotScale', P.hotScale, ...
        'colorHotMaxFactor', P.colorHotMaxFactor, ...
        'hardAlphaCutoff', P.hardAlphaCutoff, ...
        'timeLabelRef', P.timeLabelRef, ...
        'stimOnsetMs', P.stimOnsetMs, ...
        'frameRate', 10, ...
        'quality', 95, ...
        'verbose', true);
    fprintf('Saved grouped polygon movie to: %s\n', MOV_GROUP.outMovie);
end

%% Grouped activity time-series by along_GC group
if RUN_GROUP_POLYGON_TIMESERIES
    if ~exist('GGROUP', 'var')
        assert(exist(groupedFile, 'file') == 2, ...
            'Grouped polygon file missing: %s. Run RUN_GROUP_POLYGON_PREP first.', groupedFile);
        Sg = load(groupedFile);
        assert(isfield(Sg, 'G') && isstruct(Sg.G), ...
            'Grouped file %s must contain struct G.', groupedFile);
        GGROUP = Sg.G;
    end

    hTrace = plot_grouped_alonggc_timeseries( ...
        GGROUP, ...
        'plotMode', P.groupTraceMode, ...
        'timeRef', 'center', ...
        'onsetMs', P.stimOnsetMs, ...
        'cmapName', 'parula', ...
        'lineWidth', 1.8, ...
        'smoothW', P.groupTraceSmoothW, ...
        'halfMaxFrac', P.groupTraceHalfMaxFrac, ...
        'halfMaxSearchStartMs', P.groupTraceHalfSearchStartMs, ...
        'showHalfMaxMarkers', P.groupTraceShowHalfMarkers, ...
        'halfMaxMarkerSize', P.groupTraceHalfMarkerSize);
    set(hTrace.fig, 'Name', sprintf('Grouped activity traces N=%d', P.groupN), 'NumberTitle', 'off');
    if isfield(hTrace, 'summary') && istable(hTrace.summary)
        disp(hTrace.summary(:, {'groupIdx','tHalf','yHalf','yMaxSmooth'}));
    end
end

%% Shared-slope sigmoid fits for grouped traces
if RUN_GROUP_POLYGON_SIGMOID_FIT
    if ~exist('GGROUP', 'var')
        assert(exist(groupedFile, 'file') == 2, ...
            'Grouped polygon file missing: %s. Run RUN_GROUP_POLYGON_PREP first.', groupedFile);
        Sg = load(groupedFile);
        assert(isfield(Sg, 'G') && isstruct(Sg.G), ...
            'Grouped file %s must contain struct G.', groupedFile);
        GGROUP = Sg.G;
    end

    FITS = fit_grouped_alonggc_sigmoid_sharedslope( ...
        GGROUP, ...
        'timeRef', 'center', ...
        'fitStartMs', P.groupSigmoidFitStartMs, ...
        'fitEndMs', P.groupSigmoidFitEndMs, ...
        'smoothW', P.groupSigmoidSmoothW, ...
        'useAbs', P.groupSigmoidUseAbs, ...
        'minTau', P.groupSigmoidMinTau, ...
        'maxTau', P.groupSigmoidMaxTau, ...
        'initTau', P.groupSigmoidInitTau, ...
        'plotMode', 'overlay', ...
        'cmapName', 'parula', ...
        'lineWidth', 1.8, ...
        'showT50Markers', P.groupSigmoidShowT50, ...
        't50MarkerSize', P.groupSigmoidT50MarkerSize, ...
        'showT50Regression', true, ...
        'regressionX', P.groupSigmoidRegressionX, ...
        'excludeHighestXFromRegression', P.groupSigmoidExcludeHighestXRegression, ...
        'onsetMs', P.stimOnsetMs, ...
        'verbose', true);
    if ~isempty(FITS.fig)
        set(FITS.fig, 'Name', sprintf('Grouped sigmoid fits N=%d', P.groupN), 'NumberTitle', 'off');
    end
    if isfield(FITS, 'figRegression') && ~isempty(FITS.figRegression)
        set(FITS.figRegression, 'Name', sprintf('Grouped t50 regression N=%d', P.groupN), 'NumberTitle', 'off');
    end
    disp(FITS.summary(:, {'groupIdx','A','t50','t50SE','tauShared','rmse','nFitPoints','status'}));
end

%% exGauss_mod fits for grouped traces
if RUN_GROUP_POLYGON_EXGAUSS_FIT
    if ~exist('GGROUP', 'var')
        assert(exist(groupedFile, 'file') == 2, ...
            'Grouped polygon file missing: %s. Run RUN_GROUP_POLYGON_PREP first.', groupedFile);
        Sg = load(groupedFile);
        assert(isfield(Sg, 'G') && isstruct(Sg.G), ...
            'Grouped file %s must contain struct G.', groupedFile);
        GGROUP = Sg.G;
    end

    FITG = fit_grouped_alonggc_exgauss( ...
        GGROUP, ...
        'timeRef', 'center', ...
        'fitStartMs', P.groupExGaussFitStartMs, ...
        'fitEndMs', P.groupExGaussFitEndMs, ...
        'smoothW', P.groupExGaussSmoothW, ...
        'minSigma', P.groupExGaussMinSigma, ...
        'maxSigma', P.groupExGaussMaxSigma, ...
        'sigMinConsecutiveBins', P.groupSigMinRunBins, ...
        'sigSearchStartMs', P.groupSigSearchStartMs, ...
        'sigUseAbs', true, ...
        'cmapName', 'parula', ...
        'plotMode', 'overlay', ...
        'lineWidth', 1.8, ...
        'verbose', true);
    set(FITG.fig, 'Name', sprintf('Grouped exGauss fits N=%d', P.groupN), 'NumberTitle', 'off');
    disp(FITG.summary);
end

%% Optional movie rendering
if RUN_MOVIE
    if ~exist('OUT_postAffine', 'var')
        assert(exist(outValuesFile, 'file') == 2, ...
            'Post-affine values file missing: %s', outValuesFile);
        S = load(outValuesFile);
        if isfield(S, 'OUT')
            OUT_postAffine = S.OUT;
        elseif isfield(S, 'D')
            OUT_postAffine = S.D;
        else
            error('File %s must contain OUT or D.', outValuesFile);
        end
    end
    assert(isfield(OUT_postAffine, 'bins') && ~isempty(OUT_postAffine.bins), ...
        'OUT_postAffine.bins missing/empty.');
    assert(isfield(OUT_postAffine.bins, 'stream'), ...
        ['OUT_postAffine.bins.stream missing. Re-run RUN_POST_AFFINE_VALUES with the ' ...
         'updated compute_projected_delta_points_allbins.m']);

    prepNoiseFile = fullfile(cfg.resultsDir, sprintf('knn_noise_signal_prep_stim%d%s%s.mat', ...
        P.stimIDExample, tagSuffix, selectionTagSuffix));
    assert(exist(prepNoiseFile, 'file') == 2, ...
        'Noise/signal prep file missing: %s. Run RUN_PREP_NOISE_SIGNAL first.', prepNoiseFile);
    Sns = load(prepNoiseFile);
    assert(isfield(Sns, 'R') && isfield(Sns.R, 'summary') && istable(Sns.R.summary), ...
        'Prep file %s must contain R.summary table.', prepNoiseFile);
    Tns = Sns.R.summary;
    rowK = find(Tns.K == P.prepK, 1, 'first');
    assert(~isempty(rowK), ...
        'Requested prep K=%d not found in summary. Update P.prepK or rerun prep.', P.prepK);

    Kuse = double(Tns.K(rowK));
    thrUse = double(Tns.thresholdPreQ(rowK));
    assert(isfinite(thrUse) && thrUse > 0, ...
        'Invalid thresholdPreQ at K=%d in prep summary.', Kuse);
    cMaxUse = double(Tns.cMaxSuggest(rowK));
    if ~isfinite(cMaxUse) || cMaxUse <= 0
        cMaxUse = [];
    end
    siteWeightsKNN = [];
    if P.knnUseSiteWeights
        siteWeightsKNN = siteWeightsByIndexAll;
    end

    outMovie = fullfile(cfg.resultsDir, sprintf('V1_attentiondiff_movie_postaffine_K%d%s%s.mp4', ...
        Kuse, tagSuffix, selectionTagSuffix));
    MOV = make_post_affine_attention_movie( ...
        outMovie, P.stimIDExample, ALLCOORDS, RTAB384, OUT_postAffine, ...
        'K', Kuse, ...
        'siteWeightsByIndex', siteWeightsKNN, ...
        'alphaFullAt', thrUse, ...
        'colorRedAt', thrUse, ...
        'cMaxFixed', cMaxUse, ...
        'markerSize', P.plotMarkerSize, ...
        'alpha', P.plotAlpha, ...
        'bgColor', P.plotBgColor, ...
        'cLow', P.plotCLow, ...
        'cHigh', P.plotCHigh, ...
        'hotScale', P.hotScale, ...
        'colorHotMaxFactor', P.colorHotMaxFactor, ...
        'hardAlphaCutoff', P.hardAlphaCutoff, ...
        'timeLabelRef', P.timeLabelRef, ...
        'stimOnsetMs', P.stimOnsetMs, ...
        'frameRate', 10, ...
        'quality', 95, ...
        'enforceK', false, ...
        'verbose', true);
    fprintf('Saved post-affine movie to: %s\n', MOV.outMovie);
end
