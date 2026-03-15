% Attention_Line_Stimuli_V4
% Grouped V4 attention analysis on the line-task stimulus set.
%
% This mirrors the grouped along_GC path from Attention_Line_Stimuli.m:
% - significance mask from the 3-bin late-window analysis
% - post-affine signed T-D values from 20 ms response bins
% - equal-count grouping of valid (site, quartet) combinations by along_GC
% - grouped pre/post stills and grouped time-series

%% Run switches
RUN_GROUP_POLYGON_PREP = false;
RUN_GROUP_POLYGON_STILLS = false;
RUN_GROUP_POLYGON_MOVIE = false;
RUN_GROUP_POLYGON_TIMESERIES = false;
RUN_GROUP_POLYGON_SIGMOID_FIT = true;

%% Central parameters
Monkey = 1; % 1 = Nilson, 2 = Figaro
ExampleStimulus = 38; % canonical projection frame
TimeIdx3 = 3; % late window in the 3-bin summary
pThresh = 0.05; % hard significance threshold
ExcludeOverlap = true; % match the V1 attention analysis
RespMovieTag = 'd20'; % high-resolution response bins for grouped analysis
TargetWin = [300 320]; % ms, post-stim grouped still
TargetWinPre = [-100 -80]; % ms, pre-stim grouped still
GroupN = 8; % number of equal-count along_GC groups
GroupPolygonShrink = 0.8; % boundary() shrink factor for polygon outlines
PreEndMs = 0; % pre bins satisfy end <= this time
PostStartMs = 300; % post bins satisfy start >= this time
GroupPreQuantilePct = 99.98; % grouped threshold from pre-stim absolute values
GroupCMaxPostPct = 95; % grouped color max suggestion from post-stim absolute values
UseSiteWeights = true; % same reliability-weighted grouped means as V1
GroupWeightLambda = 1e-6; % stabilizer in the reliability denominator
GroupWeightUseNmatch = true; % multiply by sqrt(wY + wP)
GroupWeightClipPct = 95; % cap extreme site weights
PlotAlpha = 0.30; % used only for the grouped movie
PlotBgColor = [0.5 0.5 0.5];
PlotCLow = [0.50 0.50 0.50];
PlotCHigh = [0.85 0.05 0.05];
HotScale = true;
ColorHotMaxFactor = 8.0;
HardAlphaCutoff = true;
GroupTraceMode = 'overlay'; % 'overlay' or 'subplots'
GroupTraceSmoothW = 3; % moving-average width (bins)
GroupTraceHalfMaxFrac = 0.5; % first-crossing level as fraction of smoothed max
GroupTraceHalfSearchStartMs = 0; % search start for half-max crossing
GroupTraceShowHalfMarkers = true;
GroupTraceHalfMarkerSize = 7;
GroupSigmoidFitStartMs = 0; % shared-slope sigmoid fit start (ms)
GroupSigmoidFitEndMs = 500; % shared-slope sigmoid fit end (ms)
GroupSigmoidSmoothW = 3; % smoothing width (bins) before sigmoid fit
GroupSigmoidUseAbs = false; % fit abs(trace) if true
GroupSigmoidMinTau = 5; % lower bound on shared slope parameter tau
GroupSigmoidMaxTau = 400; % upper bound on shared slope parameter tau
GroupSigmoidInitTau = 40; % initial shared tau (ms)
GroupSigmoidShowT50 = true; % show model t50 marker per group
GroupSigmoidT50MarkerSize = 7; % marker size for model t50 points
GroupSigmoidRegressionX = 'alongMid'; % regression x-axis: 'alongMid'|'groupIdx'|'nComb'
GroupSigmoidExcludeHighestXRegression = false; % exclude highest-x point from t50 regression
FrameRate = 10; % grouped movie frame rate
Quality = 95; % grouped movie quality
ReusePostAffine = true; % load saved post-affine export when available
ReuseGrouped = true; % load saved grouped output when available

cfg = config();

if Monkey == 1
    monkeySuffix = "N";
    tallFile = 'Tall_V4_lines_N.mat';
    resp3binFile = 'SNR_capsules_N_d12.mat';
    respMovieFile = sprintf('Resp_capsules_N_%s.mat', RespMovieTag);
elseif Monkey == 2
    monkeySuffix = "F";
    tallFile = 'Tall_V4_lines_F.mat';
    resp3binFile = 'SNR_capsules_F_d12.mat';
    respMovieFile = sprintf('Resp_capsules_F_%s.mat', RespMovieTag);
else
    error('Attention_Line_Stimuli_V4:InvalidMonkey', ...
        'Monkey must be 1 (Nilson) or 2 (Figaro).');
end

tallPath = fullfile(cfg.matDir, tallFile);
resp3binPath = fullfile(cfg.matDir, resp3binFile);
respMoviePath = fullfile(cfg.matDir, respMovieFile);
outValuesFile = fullfile(cfg.resultsDir, ...
    sprintf('post_affine_delta_points_allbins_V4_%s_stim%03d_%s.mat', ...
    char(monkeySuffix), ExampleStimulus, RespMovieTag));
groupedFile = fullfile(cfg.resultsDir, ...
    sprintf('grouped_alonggc_polygons_V4_%s_stim%03d_N%d_%s.mat', ...
    char(monkeySuffix), ExampleStimulus, GroupN, RespMovieTag));
outMovieGrouped = fullfile(cfg.resultsDir, ...
    sprintf('V4_attentiondiff_groupedpolygons_%s_N%d_%s_stim%03d.mp4', ...
    char(monkeySuffix), GroupN, RespMovieTag, ExampleStimulus));

assert(exist(tallPath, 'file') == 2, ...
    'Missing %s. Run Line_Stimuli_V4.m first.', tallPath);
assert(exist(resp3binPath, 'file') == 2, ...
    'Missing %s. Need the 3-bin response summary for the significance mask.', resp3binPath);
assert(exist(respMoviePath, 'file') == 2, ...
    'Missing %s. Need the %s response file for grouped analysis.', respMoviePath, RespMovieTag);

Sgeo = load(tallPath);
assert(isfield(Sgeo, 'Tall_V4') && isstruct(Sgeo.Tall_V4), ...
    '%s must contain struct Tall_V4.', tallFile);
assert(isfield(Sgeo, 'ALLCOORDS') && isfield(Sgeo, 'RTAB384') && isfield(Sgeo, 'RFrange'), ...
    '%s must contain Tall_V4, ALLCOORDS, RTAB384, and RFrange.', tallFile);

Tall_V4 = Sgeo.Tall_V4;
ALLCOORDS = Sgeo.ALLCOORDS;
RTAB384 = Sgeo.RTAB384;
RFrange = Sgeo.RFrange(:);
nV4 = numel(RFrange);
siteRows = (1:nV4).';

%% Load 3-bin response + normalization and get baseline OUT3
S3 = load(resp3binPath);
assert(isfield(S3, 'R') && isstruct(S3.R), ...
    '%s must contain struct R.', resp3binFile);
R3_full = S3.R;
R3dec = R3_full;
R3dec.meanAct = R3_full.meanAct(RFrange, :, :);
R3dec.meanSqAct = R3_full.meanSqAct(RFrange, :, :);
if ismatrix(R3_full.nTrials) && size(R3_full.nTrials,1) >= max(RFrange)
    R3dec.nTrials = R3_full.nTrials(RFrange, :);
else
    R3dec.nTrials = R3_full.nTrials;
end

SNRnorm = compute_snr_per_color_sites(R3dec, Tall_V4, siteRows, 'Verbose', false);
opts3 = struct('v1Sites', siteRows, 'timeIdx', TimeIdx3, ...
    'excludeOverlap', ExcludeOverlap, 'verbose', false, 'epsDen', 1e-6);
OUT3 = attention_modulation_V1_3bin(R3dec, Tall_V4, SNRnorm, opts3);

sigSiteByIndex = isfinite(OUT3.pValueTD(:)) & (OUT3.pValueTD(:) < pThresh);
fprintf('V4 grouped attention mask: %d / %d sites with pTD < %.3f\n', ...
    nnz(sigSiteByIndex), nV4, pThresh);
assert(any(sigSiteByIndex), 'No V4 sites passed pTD < %.3f.', pThresh);

siteWeightsByIndex = [];
if UseSiteWeights
    assert(all(isfield(OUT3, {'muT','muD','varT','varD'})), ...
        'OUT3 must contain muT/muD/varT/varD for site weighting.');
    dSite = double(OUT3.muT(:) - OUT3.muD(:));
    varSite = 0.5 * (double(OUT3.varT(:)) + double(OUT3.varD(:)));
    varSite(~isfinite(varSite) | varSite < 0) = 0;
    den = sqrt(varSite + max(GroupWeightLambda, 0));
    den(~isfinite(den) | den <= 0) = NaN;
    siteWeightsByIndex = abs(dSite) ./ den;

    if GroupWeightUseNmatch
        assert(all(isfield(OUT3, {'wY','wP'})), ...
            'OUT3 must contain wY/wP when GroupWeightUseNmatch=true.');
        nMatch = double(OUT3.wY(:) + OUT3.wP(:));
        nMatch(~isfinite(nMatch) | nMatch < 0) = 0;
        siteWeightsByIndex = siteWeightsByIndex .* sqrt(nMatch);
    end

    siteWeightsByIndex(~isfinite(siteWeightsByIndex) | siteWeightsByIndex < 0) = 0;
    siteWeightsByIndex(~sigSiteByIndex) = 0;

    wPos = siteWeightsByIndex(siteWeightsByIndex > 0);
    if ~isempty(wPos) && isfinite(GroupWeightClipPct) && ...
            GroupWeightClipPct > 0 && GroupWeightClipPct < 100
        wCap = prctile(wPos, GroupWeightClipPct);
        if isfinite(wCap) && wCap > 0
            siteWeightsByIndex = min(siteWeightsByIndex, wCap);
        end
    end

    wPos = siteWeightsByIndex(siteWeightsByIndex > 0);
    if isempty(wPos)
        warning(['V4 grouped attention requested site weighting but no positive ' ...
                 'site weights were obtained. Falling back to unweighted grouped means.']);
        siteWeightsByIndex = [];
    else
        fprintf(['Fixed site weights ready: nPos=%d | min/median/max=%.6g / %.6g / %.6g | ' ...
                 'clip p%.2f\n'], ...
            numel(wPos), min(wPos), median(wPos), max(wPos), GroupWeightClipPct);
    end
end

%% Load high-resolution response bins
Sm = load(respMoviePath);
assert(isfield(Sm, 'R') && isstruct(Sm.R), ...
    '%s must contain struct R.', respMovieFile);
Rmovie_full = Sm.R;
Rmovie = Rmovie_full;
Rmovie.meanAct = Rmovie_full.meanAct(RFrange, :, :);
Rmovie.meanSqAct = Rmovie_full.meanSqAct(RFrange, :, :);
if ismatrix(Rmovie_full.nTrials) && size(Rmovie_full.nTrials,1) >= max(RFrange)
    Rmovie.nTrials = Rmovie_full.nTrials(RFrange, :);
else
    Rmovie.nTrials = Rmovie_full.nTrials;
end

assert(isfield(Rmovie, 'timeWindows') && size(Rmovie.timeWindows,2) == 2, ...
    '%s must contain R.timeWindows as [nBins x 2].', respMovieFile);
assert(size(Rmovie.timeWindows,1) > 3, ...
    '%s must contain many short windows for grouped analysis.', respMovieFile);

[~, timeIdxTarget] = min(sum(abs(Rmovie.timeWindows - TargetWin),2));
winTarget = double(Rmovie.timeWindows(timeIdxTarget,:));
[~, timeIdxPre] = min(sum(abs(Rmovie.timeWindows - TargetWinPre),2));
winPre = double(Rmovie.timeWindows(timeIdxPre,:));

%% Post-affine export for all bins
OUT_postAffine = [];
if ReusePostAffine && exist(outValuesFile, 'file') == 2
    Svals = load(outValuesFile);
    if isfield(Svals, 'OUT') && isstruct(Svals.OUT)
        OUT_postAffine = Svals.OUT;
    elseif isfield(Svals, 'D') && isstruct(Svals.D)
        OUT_postAffine = Svals.D;
    end
end

if isempty(OUT_postAffine)
    OUT_postAffine = compute_projected_delta_points_allbins( ...
        ExampleStimulus, Tall_V4, ALLCOORDS, RTAB384, Rmovie, SNRnorm, ...
        'siteRange', siteRows, ...
        'excludeOverlap', ExcludeOverlap, ...
        'stimIdx', 1:384, ...
        'sigSiteMask', sigSiteByIndex, ...
        'saveFile', outValuesFile, ...
        'verbose', true);
    fprintf('Saved V4 post-affine values to: %s\n', outValuesFile);
else
    fprintf('Loaded V4 post-affine values from: %s\n', outValuesFile);
end

%% Grouped along_GC polygon prep
GGROUP = [];
if ReuseGrouped && exist(groupedFile, 'file') == 2
    Sg = load(groupedFile);
    if isfield(Sg, 'G') && isstruct(Sg.G)
        GGROUP = Sg.G;
    end
end

if RUN_GROUP_POLYGON_PREP || isempty(GGROUP)
    weightsForGrouping = [];
    if UseSiteWeights
        weightsForGrouping = siteWeightsByIndex;
    end

    GGROUP = build_grouped_alonggc_polygons_allbins( ...
        OUT_postAffine, Tall_V4, ...
        'nGroups', GroupN, ...
        'sigSiteByIndex', sigSiteByIndex, ...
        'siteWeightsByIndex', weightsForGrouping, ...
        'preEndMs', PreEndMs, ...
        'postStartMs', PostStartMs, ...
        'preQuantilePct', GroupPreQuantilePct, ...
        'cMaxPostPct', GroupCMaxPostPct, ...
        'polygonShrink', GroupPolygonShrink, ...
        'saveFile', groupedFile, ...
        'verbose', true);

    fprintf('Saved grouped along_GC polygons to: %s\n', groupedFile);
else
    fprintf('Loaded grouped along_GC polygons from: %s\n', groupedFile);
end

if isfield(GGROUP, 'groupSummary') && istable(GGROUP.groupSummary)
    disp(GGROUP.groupSummary(:, {'groupIdx','nComb','alongMin','alongMax'}));
end

%% Grouped polygon stills (pre/post)
if RUN_GROUP_POLYGON_STILLS
    assert(isfield(GGROUP, 'groupMeanSigned') && isfield(GGROUP, 'timeWindows') && ...
        isfield(GGROUP, 'groupPolygons'), ...
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

    showStimTarget = (winTargetG(1) >= 0);
    showStimPre = (winPreG(1) >= 0);

    hGroupTarget = plot_grouped_alonggc_polygon_frame( ...
        ExampleStimulus, ALLCOORDS, RTAB384, GGROUP, timeIdxTargetG, ...
        'alphaFullAt', thrGrp, ...
        'colorRedAt', thrGrp, ...
        'cMaxFixed', cMaxGrp, ...
        'alpha', 1, ...
        'bgColor', PlotBgColor, ...
        'cLow', PlotCLow, ...
        'cHigh', PlotCHigh, ...
        'hotScale', HotScale, ...
        'colorHotMaxFactor', ColorHotMaxFactor, ...
        'hardAlphaCutoff', HardAlphaCutoff, ...
        'timeLabelRef', 'start', ...
        'showStimulus', showStimTarget);
    figure(hGroupTarget.fig);
    set(hGroupTarget.fig, 'Name', sprintf('V4 grouped polygons N=%d | post-stim', GroupN), ...
        'NumberTitle', 'off');
    title(hGroupTarget.ax, sprintf('V4 grouped polygons (N=%d) | [%d %d] ms | stim %d', ...
        GroupN, winTargetG(1), winTargetG(2), ExampleStimulus), 'Color', 'w');

    hGroupPre = plot_grouped_alonggc_polygon_frame( ...
        ExampleStimulus, ALLCOORDS, RTAB384, GGROUP, timeIdxPreG, ...
        'alphaFullAt', thrGrp, ...
        'colorRedAt', thrGrp, ...
        'cMaxFixed', cMaxGrp, ...
        'alpha', 1, ...
        'bgColor', PlotBgColor, ...
        'cLow', PlotCLow, ...
        'cHigh', PlotCHigh, ...
        'hotScale', HotScale, ...
        'colorHotMaxFactor', ColorHotMaxFactor, ...
        'hardAlphaCutoff', HardAlphaCutoff, ...
        'timeLabelRef', 'start', ...
        'showStimulus', showStimPre);
    figure(hGroupPre.fig);
    set(hGroupPre.fig, 'Name', sprintf('V4 grouped polygons N=%d | pre-stim', GroupN), ...
        'NumberTitle', 'off');
    title(hGroupPre.ax, sprintf('V4 grouped polygons (N=%d) | [%d %d] ms | stim %d', ...
        GroupN, winPreG(1), winPreG(2), ExampleStimulus), 'Color', 'w');

    fprintf(['Grouped stills (N=%d): pre >thr %.2f%%, post >thr %.2f%% | ' ...
             'pre routes T:%d D:%d, post routes T:%d D:%d\n'], ...
        GroupN, 100*hGroupPre.fracAboveThreshold, 100*hGroupTarget.fracAboveThreshold, ...
        hGroupPre.nTargetGroups, hGroupPre.nDistrGroups, ...
        hGroupTarget.nTargetGroups, hGroupTarget.nDistrGroups);
end

%% Grouped polygon movie rendering
if RUN_GROUP_POLYGON_MOVIE
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

    MOV_GROUP = make_grouped_alonggc_polygon_movie( ...
        outMovieGrouped, ExampleStimulus, ALLCOORDS, RTAB384, GGROUP, ...
        'alphaFullAt', thrGrp, ...
        'colorRedAt', thrGrp, ...
        'cMaxFixed', cMaxGrp, ...
        'alpha', PlotAlpha, ...
        'bgColor', PlotBgColor, ...
        'cLow', PlotCLow, ...
        'cHigh', PlotCHigh, ...
        'hotScale', HotScale, ...
        'colorHotMaxFactor', ColorHotMaxFactor, ...
        'hardAlphaCutoff', HardAlphaCutoff, ...
        'timeLabelRef', 'start', ...
        'stimOnsetMs', 0, ...
        'frameRate', FrameRate, ...
        'quality', Quality, ...
        'verbose', true);
    fprintf('Saved grouped polygon movie to: %s\n', MOV_GROUP.outMovie);
end

%% Grouped activity time-series by along_GC group
if RUN_GROUP_POLYGON_TIMESERIES
    hTrace = plot_grouped_alonggc_timeseries( ...
        GGROUP, ...
        'plotMode', GroupTraceMode, ...
        'timeRef', 'center', ...
        'onsetMs', 0, ...
        'cmapName', 'parula', ...
        'lineWidth', 1.8, ...
        'smoothW', GroupTraceSmoothW, ...
        'halfMaxFrac', GroupTraceHalfMaxFrac, ...
        'halfMaxSearchStartMs', GroupTraceHalfSearchStartMs, ...
        'showHalfMaxMarkers', GroupTraceShowHalfMarkers, ...
        'halfMaxMarkerSize', GroupTraceHalfMarkerSize);
    set(hTrace.fig, 'Name', sprintf('V4 grouped activity traces N=%d', GroupN), ...
        'NumberTitle', 'off');
    if isfield(hTrace, 'summary') && istable(hTrace.summary)
        disp(hTrace.summary(:, {'groupIdx','tHalf','yHalf','yMaxSmooth'}));
    end
end

%% Shared-slope sigmoid fits for grouped traces
if RUN_GROUP_POLYGON_SIGMOID_FIT
    FITS = fit_grouped_alonggc_sigmoid_sharedslope( ...
        GGROUP, ...
        'timeRef', 'center', ...
        'fitStartMs', GroupSigmoidFitStartMs, ...
        'fitEndMs', GroupSigmoidFitEndMs, ...
        'smoothW', GroupSigmoidSmoothW, ...
        'useAbs', GroupSigmoidUseAbs, ...
        'minTau', GroupSigmoidMinTau, ...
        'maxTau', GroupSigmoidMaxTau, ...
        'initTau', GroupSigmoidInitTau, ...
        'plotMode', 'overlay', ...
        'cmapName', 'parula', ...
        'lineWidth', 1.8, ...
        'showT50Markers', GroupSigmoidShowT50, ...
        't50MarkerSize', GroupSigmoidT50MarkerSize, ...
        'showT50Regression', true, ...
        'regressionX', GroupSigmoidRegressionX, ...
        'excludeHighestXFromRegression', GroupSigmoidExcludeHighestXRegression, ...
        'onsetMs', 0, ...
        'verbose', true);
    if ~isempty(FITS.fig)
        set(FITS.fig, 'Name', sprintf('V4 grouped sigmoid fits N=%d', GroupN), ...
            'NumberTitle', 'off');
    end
    if isfield(FITS, 'figRegression') && ~isempty(FITS.figRegression)
        set(FITS.figRegression, 'Name', sprintf('V4 grouped t50 regression N=%d', GroupN), ...
            'NumberTitle', 'off');
    end
    disp(FITS.summary(:, {'groupIdx','A','t50','t50SE','tauShared','rmse','nFitPoints','status'}));
end

V4Grouped = struct();
V4Grouped.Monkey = Monkey;
V4Grouped.monkeySuffix = monkeySuffix;
V4Grouped.RFrange = RFrange;
V4Grouped.OUT3 = OUT3;
V4Grouped.sigSiteByIndex = sigSiteByIndex;
V4Grouped.siteWeightsByIndex = siteWeightsByIndex;
V4Grouped.OUT_postAffine = OUT_postAffine;
V4Grouped.GGROUP = GGROUP;
V4Grouped.groupedFile = groupedFile;
if exist('FITS', 'var')
    V4Grouped.FITS = FITS;
end
