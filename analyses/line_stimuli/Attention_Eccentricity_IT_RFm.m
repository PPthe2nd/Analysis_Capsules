% Attention_Eccentricity_IT_RFm
% Grouped IT attention analysis on the line-task stimulus set using
% RF.m-center eccentricity bins.
%
% This keeps the same post-affine delta export as the V1/V4 attention
% movies, but groups valid (site, quartet) combinations into 8 equal-count
% bins by IT RF eccentricity. By default it mirrors the active V1 grouped
% path in two ways: all IT responsive sites are included, and grouped means
% use pure |d'| site weights.

%% Run switches
RUN_GROUP_PREP = true;
RUN_GROUP_STILLS = true;
RUN_GROUP_POLYGON_MOVIE = false;
RUN_GROUP_TIMESERIES = true;
RUN_GROUP_POLYGON_SIGMOID_FIT = true;

%% Central parameters
Monkey = 1; % 1 = Nilson, 2 = Figaro
ExampleStimulus = 38; % canonical projection frame
TimeIdx3 = 3; % late window in the 3-bin summary
pThresh = 0.05; % hard significance threshold
ExcludeOverlap = true; % match the V1/V4 attention analysis
RespMovieTag = 'd20'; % high-resolution response bins for grouped analysis
TargetWin = [300 320]; % ms, post-stim grouped still
TargetWinPre = [-100 -80]; % ms, pre-stim grouped still
GroupN = 8; % number of equal-count RF-eccentricity groups
GroupPolygonShrink = 0.8; % boundary() shrink factor for polygon outlines
PreEndMs = 0; % pre bins satisfy end <= this time
PostStartMs = 300; % post bins satisfy start >= this time
GroupPreQuantilePct = 99.98; % grouped threshold from pre-stim absolute values
GroupCMaxPostPct = 95; % grouped color max suggestion from post-stim absolute values
SiteInclusionMode = 'responsive'; % 'responsive' or 'significant'
UseSiteWeights = true;
SiteWeightMode = 'dprime'; % 'dprime' or 'reliability'
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
GroupTraceMode = 'overlay'; % use overlay so labels stay compact
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
GroupSigmoidRegressionX = 'alongMid'; % uses eccentricity midpoint in this IT grouped output
GroupSigmoidExcludeHighestXRegression = false; % exclude highest-x point from t50 regression
FrameRate = 10; % grouped movie frame rate
Quality = 95; % grouped movie quality
MinObjectStim = 1; % IT validity gate, same as Attention_Histogram_IT_RFm
RespThr = 0.7; % IT early quartet-responsiveness threshold
TopKQuartets = 5; % IT responsiveness score uses top-k quartets
ReusePostAffine = true; % load saved post-affine export when available
ReuseGrouped = true; % load saved grouped output when available

cfg = config();

if Monkey == 1
    monkeySuffix = "N";
    tallFile = 'Tall_IT_lines_N.mat';
    resp3binFile = 'SNR_capsules_N_d12.mat';
    respMovieFile = sprintf('Resp_capsules_N_%s.mat', RespMovieTag);
elseif Monkey == 2
    monkeySuffix = "F";
    tallFile = 'Tall_IT_lines_F.mat';
    resp3binFile = 'SNR_capsules_F_d12.mat';
    respMovieFile = sprintf('Resp_capsules_F_%s.mat', RespMovieTag);
else
    error('Attention_Eccentricity_IT_RFm:InvalidMonkey', ...
        'Monkey must be 1 (Nilson) or 2 (Figaro).');
end

gateSuffix = sprintf('pTD%.3f_top%d_%.2f', pThresh, TopKQuartets, RespThr);
siteInclusionMode = lower(string(SiteInclusionMode));
siteWeightMode = lower(string(SiteWeightMode));
assert(any(siteInclusionMode == ["responsive","significant"]), ...
    'SiteInclusionMode must be ''responsive'' or ''significant''.');
assert(any(siteWeightMode == ["dprime","reliability"]), ...
    'SiteWeightMode must be ''dprime'' or ''reliability''.');

switch char(siteInclusionMode)
    case 'responsive'
        inclusionTag = sprintf('allresp_top%d_%.2f', TopKQuartets, RespThr);
        selectionLabel = sprintf('responsive sites (valid + object RF + top-%d > %.2f)', ...
            TopKQuartets, RespThr);
    case 'significant'
        inclusionTag = sprintf('pTD%.3f_top%d_%.2f', pThresh, TopKQuartets, RespThr);
        selectionLabel = sprintf('attention-significant responsive sites (pTD < %.3f)', pThresh);
end
switch char(siteWeightMode)
    case 'dprime'
        weightTag = 'wdprime';
    case 'reliability'
        weightTag = 'wrel';
end
gateSuffix = strrep(sprintf('%s_%s', inclusionTag, weightTag), '.', 'p');

tallPath = fullfile(cfg.matDir, tallFile);
resp3binPath = fullfile(cfg.matDir, resp3binFile);
respMoviePath = fullfile(cfg.matDir, respMovieFile);
outValuesFile = fullfile(cfg.resultsDir, ...
    sprintf('post_affine_delta_points_allbins_IT_%s_stim%03d_%s_%s.mat', ...
    char(monkeySuffix), ExampleStimulus, RespMovieTag, gateSuffix));
groupedFile = fullfile(cfg.resultsDir, ...
    sprintf('grouped_rfecc_polygons_IT_%s_stim%03d_N%d_%s_%s.mat', ...
    char(monkeySuffix), ExampleStimulus, GroupN, RespMovieTag, gateSuffix));
outMovieGrouped = fullfile(cfg.resultsDir, ...
    sprintf('IT_attentiondiff_grouped_rfecc_%s_N%d_%s_stim%03d_%s.mp4', ...
    char(monkeySuffix), GroupN, RespMovieTag, ExampleStimulus, gateSuffix));

assert(exist(tallPath, 'file') == 2, ...
    'Missing %s. Run Line_Stimuli_IT.m first.', tallPath);
assert(exist(resp3binPath, 'file') == 2, ...
    'Missing %s. Need the 3-bin response summary for IT normalization.', resp3binPath);
assert(exist(respMoviePath, 'file') == 2, ...
    'Missing %s. Need the %s response file for grouped analysis.', respMoviePath, RespMovieTag);

Sgeo = load(tallPath);
assert(isfield(Sgeo, 'Tall_IT') && isstruct(Sgeo.Tall_IT), ...
    '%s must contain struct Tall_IT.', tallFile);
assert(isfield(Sgeo, 'ALLCOORDS') && isfield(Sgeo, 'RTAB384') && isfield(Sgeo, 'RFrange'), ...
    '%s must contain Tall_IT, ALLCOORDS, RTAB384, and RFrange.', tallFile);

Tall_IT = Sgeo.Tall_IT;
ALLCOORDS = Sgeo.ALLCOORDS;
RTAB384 = Sgeo.RTAB384;
RFrange = Sgeo.RFrange(:);
nIT = numel(RFrange);
siteRows = (1:nIT).';

Sel = Attention_Histogram_IT_RFm(Monkey, struct( ...
    'PlotFigure', false, ...
    'timeIdx3', TimeIdx3, ...
    'pThresh', pThresh, ...
    'excludeOverlap', ExcludeOverlap, ...
    'MinObjectStim', MinObjectStim, ...
    'RespThr', RespThr, ...
    'TopKQuartets', TopKQuartets));

ATT = Sel.ATT;
keep = Sel.keep(:);
isSig = Sel.isSig(:);
selectedMaskByIndex = keep;
if siteInclusionMode == "significant"
    selectedMaskByIndex = isSig;
end

fprintf('IT eccentricity grouping gate: responsive %d / %d | significant %d / %d responsive\n', ...
    nnz(keep), nIT, nnz(isSig), nnz(keep));
fprintf('Site inclusion mode: %s | kept %d / %d sites\n', ...
    selectionLabel, nnz(selectedMaskByIndex), nIT);
assert(any(selectedMaskByIndex), 'No IT sites passed the selected inclusion gate.');

R3_full = load_capsules_struct_exclusion_aware(resp3binPath, monkeySuffix, 'cfg', cfg);
R3dec = R3_full;
R3dec.meanAct = R3_full.meanAct(RFrange, :, :);
R3dec.meanSqAct = R3_full.meanSqAct(RFrange, :, :);
if ismatrix(R3_full.nTrials) && size(R3_full.nTrials,1) >= max(RFrange)
    R3dec.nTrials = R3_full.nTrials(RFrange, :);
else
    R3dec.nTrials = R3_full.nTrials;
end
SNRnorm = compute_snr_per_color_sites(R3dec, Tall_IT, siteRows, 'Verbose', false);

siteWeightsByIndex = [];
if UseSiteWeights
    switch char(siteWeightMode)
        case 'dprime'
            assert(isfield(ATT, 'dprime'), ...
                'ATT must contain dprime for dprime weighting.');
            siteWeightsByIndex = abs(double(ATT.dprime(:)));
        case 'reliability'
            assert(all(isfield(ATT, {'muT','muD','varT','varD'})), ...
                'ATT must contain muT/muD/varT/varD for site weighting.');
            dSite = double(ATT.muT(:) - ATT.muD(:));
            varSite = 0.5 * (double(ATT.varT(:)) + double(ATT.varD(:)));
            varSite(~isfinite(varSite) | varSite < 0) = 0;
            den = sqrt(varSite + max(GroupWeightLambda, 0));
            den(~isfinite(den) | den <= 0) = NaN;
            siteWeightsByIndex = abs(dSite) ./ den;
    end

    applyNmatchWeight = GroupWeightUseNmatch;
    if siteWeightMode == "dprime"
        applyNmatchWeight = false;
    end
    if applyNmatchWeight
        assert(all(isfield(ATT, {'wY','wP'})), ...
            'ATT must contain wY/wP when GroupWeightUseNmatch=true.');
        nMatch = double(ATT.wY(:) + ATT.wP(:));
        nMatch(~isfinite(nMatch) | nMatch < 0) = 0;
        siteWeightsByIndex = siteWeightsByIndex .* sqrt(nMatch);
    end

    siteWeightsByIndex(~isfinite(siteWeightsByIndex) | siteWeightsByIndex < 0) = 0;
    siteWeightsByIndex(~selectedMaskByIndex) = 0;

    wPos = siteWeightsByIndex(siteWeightsByIndex > 0);
    if ~isempty(wPos) && isfinite(GroupWeightClipPct) && GroupWeightClipPct > 0 && GroupWeightClipPct < 100
        wCap = prctile(wPos, GroupWeightClipPct);
        if isfinite(wCap) && wCap > 0
            siteWeightsByIndex = min(siteWeightsByIndex, wCap);
        end
    end

    wPos = siteWeightsByIndex(siteWeightsByIndex > 0);
    if isempty(wPos)
        warning('IT eccentricity grouping requested site weighting but no positive site weights were obtained. Falling back to unweighted means.');
        siteWeightsByIndex = [];
    else
        fprintf(['Fixed IT site weights ready (%s): nPos=%d | min/median/max=%.6g / %.6g / %.6g | ' ...
                 'clip p%.2f\n'], ...
            char(siteWeightMode), numel(wPos), min(wPos), median(wPos), max(wPos), GroupWeightClipPct);
    end
end

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
        ExampleStimulus, Tall_IT, ALLCOORDS, RTAB384, Rmovie, SNRnorm, ...
        'siteRange', siteRows, ...
        'excludeOverlap', ExcludeOverlap, ...
        'stimIdx', 1:384, ...
        'sigSiteMask', selectedMaskByIndex, ...
        'saveFile', outValuesFile, ...
        'verbose', true);
    fprintf('Saved IT post-affine values to: %s\n', outValuesFile);
else
    fprintf('Loaded IT post-affine values from: %s\n', outValuesFile);
end

rfEcc = load_it_rf_eccentricity_local(cfg, Monkey, RFrange);
assert(numel(rfEcc) == nIT, 'RF eccentricity vector length mismatch.');

G = [];
if ReuseGrouped && exist(groupedFile, 'file') == 2
    Sg = load(groupedFile);
    if isfield(Sg, 'G') && isstruct(Sg.G)
        G = Sg.G;
    end
end

if RUN_GROUP_PREP || isempty(G)
    G = build_grouped_eccentricity_polygons_allbins( ...
        OUT_postAffine, rfEcc, ...
        'nGroups', GroupN, ...
        'sigSiteByIndex', selectedMaskByIndex, ...
        'siteWeightsByIndex', siteWeightsByIndex, ...
        'preEndMs', PreEndMs, ...
        'postStartMs', PostStartMs, ...
        'preQuantilePct', GroupPreQuantilePct, ...
        'cMaxPostPct', GroupCMaxPostPct, ...
        'polygonShrink', GroupPolygonShrink, ...
        'saveFile', groupedFile, ...
        'verbose', true);
    fprintf('Saved grouped RF-eccentricity polygons to: %s\n', groupedFile);
else
    fprintf('Loaded grouped RF-eccentricity polygons from: %s\n', groupedFile);
end

if isfield(G, 'groupSummary') && istable(G.groupSummary)
    disp(G.groupSummary(:, {'groupIdx','nComb','eccMin','eccMax'}));
end

tw = double(G.timeWindows);
[~, timeIdxTargetG] = min(sum(abs(tw - TargetWin),2));
winTargetG = double(tw(timeIdxTargetG,:));
[~, timeIdxPreG] = min(sum(abs(tw - TargetWinPre),2));
winPreG = double(tw(timeIdxPreG,:));

thrGrp = [];
cMaxGrp = [];
if isfield(G, 'calibration') && isstruct(G.calibration)
    if isfield(G.calibration, 'thresholdPreQ')
        thrGrp = double(G.calibration.thresholdPreQ);
    end
    if isfield(G.calibration, 'cMaxSuggest')
        cMaxGrp = double(G.calibration.cMaxSuggest);
    end
end
if ~isfinite(thrGrp) || thrGrp <= 0
    thrGrp = 0.2;
end
if ~isfinite(cMaxGrp) || cMaxGrp <= 0
    cMaxGrp = [];
end

if RUN_GROUP_STILLS
    hGroupTarget = plot_grouped_alonggc_polygon_frame( ...
        ExampleStimulus, ALLCOORDS, RTAB384, G, timeIdxTargetG, ...
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
        'showStimulus', true);
    figure(hGroupTarget.fig);
    set(hGroupTarget.fig, 'Name', sprintf('IT grouped RF-ecc N=%d | post-stim', GroupN), ...
        'NumberTitle', 'off');
    title(hGroupTarget.ax, sprintf('IT grouped RF eccentricity (N=%d) | [%d %d] ms | stim %d', ...
        GroupN, winTargetG(1), winTargetG(2), ExampleStimulus), 'Color', 'w');

    hGroupPre = plot_grouped_alonggc_polygon_frame( ...
        ExampleStimulus, ALLCOORDS, RTAB384, G, timeIdxPreG, ...
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
        'showStimulus', true);
    figure(hGroupPre.fig);
    set(hGroupPre.fig, 'Name', sprintf('IT grouped RF-ecc N=%d | pre-stim', GroupN), ...
        'NumberTitle', 'off');
    title(hGroupPre.ax, sprintf('IT grouped RF eccentricity (N=%d) | [%d %d] ms | stim %d', ...
        GroupN, winPreG(1), winPreG(2), ExampleStimulus), 'Color', 'w');
end

if RUN_GROUP_POLYGON_MOVIE
    MOV_GROUP = make_grouped_alonggc_polygon_movie( ...
        outMovieGrouped, ExampleStimulus, ALLCOORDS, RTAB384, G, ...
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
    fprintf('Saved grouped RF-ecc movie to: %s\n', MOV_GROUP.outMovie);
end

if RUN_GROUP_TIMESERIES
    hTrace = plot_grouped_alonggc_timeseries( ...
        G, ...
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
    set(hTrace.fig, 'Name', sprintf('IT grouped RF eccentricity traces N=%d', GroupN), ...
        'NumberTitle', 'off');
    if isscalar(hTrace.ax)
        title(hTrace.ax, sprintf('IT grouped RF eccentricity traces (N=%d)', GroupN));
    elseif ~isempty(hTrace.tiled)
        title(hTrace.tiled, sprintf('IT grouped RF eccentricity traces (N=%d)', GroupN));
    end

    summaryTable = table( ...
        G.groupSummary.groupIdx, G.groupSummary.nComb, ...
        G.groupSummary.eccMin, G.groupSummary.eccMax, ...
        hTrace.summary.tHalf, hTrace.summary.yHalf, hTrace.summary.yMaxSmooth, ...
        'VariableNames', {'groupIdx','nComb','eccMin','eccMax','tHalf','yHalf','yMaxSmooth'});
    disp(summaryTable);
end

if RUN_GROUP_POLYGON_SIGMOID_FIT
    FITS = fit_grouped_alonggc_sigmoid_sharedslope( ...
        G, ...
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
        set(FITS.fig, 'Name', sprintf('IT grouped RF eccentricity sigmoid fits N=%d', GroupN), ...
            'NumberTitle', 'off');
    end
    if isfield(FITS, 'figRegression') && ~isempty(FITS.figRegression)
        set(FITS.figRegression, 'Name', sprintf('IT grouped RF eccentricity t50 regression N=%d', GroupN), ...
            'NumberTitle', 'off');
    end
    fitCols = {'groupIdx','A','t50','t50SE','tauShared','rmse','nFitPoints','status'};
    if ismember('eccMin', FITS.summary.Properties.VariableNames)
        fitCols = [{'groupIdx','eccMin','eccMax'}, fitCols(2:end)];
    end
    disp(FITS.summary(:, fitCols));
end

function rfEcc = load_it_rf_eccentricity_local(cfg, Monkey, RFrange)
if Monkey == 1
    fileName = 'THINGS_RF1s_N.mat';
    legacySubdir = 'Mr Nilson';
elseif Monkey == 2
    fileName = 'THINGS_RF1s_F.mat';
    legacySubdir = 'Figaro';
else
    error('load_it_rf_eccentricity_local:InvalidMonkey', ...
        'Monkey must be 1 (Nilson) or 2 (Figaro).');
end

candidates = { ...
    fullfile(cfg.matDir, fileName), ...
    fullfile(cfg.dataRoot, legacySubdir, fileName)};

rfPath = '';
for i = 1:numel(candidates)
    if exist(candidates{i}, 'file') == 2
        rfPath = candidates{i};
        break;
    end
end
assert(~isempty(rfPath), 'Could not find RF file %s.', fileName);

Srf = load(rfPath, 'all_centrex', 'all_centrey');
assert(isfield(Srf, 'all_centrex') && isfield(Srf, 'all_centrey'), ...
    '%s must contain all_centrex/all_centrey.', rfPath);

x = double(Srf.all_centrex(:));
y = double(Srf.all_centrey(:));
assert(all(RFrange >= 1 & RFrange <= numel(x)) && all(RFrange <= numel(y)), ...
    'RFrange exceeds RF center arrays in %s.', rfPath);

rfEcc = hypot(x(RFrange), y(RFrange));
end
