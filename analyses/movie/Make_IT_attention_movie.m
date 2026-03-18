% Render an IT attention-effect movie for the line-task stimulus set.
%
% This mirrors the active V4 post-affine movie path, but uses the IT
% RF.m-center geometry and the IT validity/responsiveness gate from
% Attention_Histogram_IT_RFm.

Monkey = 1; % 1 = Nilson, 2 = Figaro
ExampleStimulus = 38; % canonical projection frame for the movie
TimeIdx3 = 3; % late window in the 3-bin summary
pTDthr = 0.05; % hard significance threshold for site inclusion
ExcludeOverlap = true; % match the V1/V4 attention analysis
RespMovieTag = 'd20'; % use 20 ms bins for the time-resolved movie
FrameRate = 10; % movie frame rate
Quality = 95; % MPEG-4 writer quality
PlotMarkerSize = 8; % same size as V1/V4 attention movies
PlotAlpha = 0.30; % constant marker alpha above threshold
PlotBgColor = [0.5 0.5 0.5];
PlotCLow = [0.50 0.50 0.50];
PlotCHigh = [0.85 0.05 0.05];
HotScale = true;
ColorHotMaxFactor = 8.0;
HardAlphaCutoff = true;
PrepK = 80; % K used for the final post-affine movie; larger K visibly oversmooths IT late frames
PrepKList = [10 50 80 200 500 1000]; % calibration sweep
PreEndMs = 0; % pre bins satisfy window end <= this time
PostStartMs = 300; % post bins satisfy window start >= this time
PreQuantilePct = 99.98; % threshold from pre-stim pooled absolute values
AlphaFloorPct = 80; % retained for parity with V1/V4 prep summary
CMaxPostPct = 95; % color max suggestion from post-stim pooled values
UseSiteWeights = true; % reliability-weighted KNN averaging
GroupWeightLambda = 1e-6; % stabilizer in the reliability denominator
GroupWeightUseNmatch = true; % multiply weights by sqrt(wY + wP)
GroupWeightClipPct = 95; % cap extreme positive weights
MinObjectStim = 1; % IT validity gate, same as Attention_Histogram_IT_RFm
RespThr = 0.7; % IT early quartet-responsiveness threshold
TopKQuartets = 5; % IT responsiveness score uses top-k quartets
ReusePostAffine = true; % load saved post-affine export when available
ForceRender = false; % if false, skip when the output movie already exists

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
    error('Make_IT_attention_movie:InvalidMonkey', ...
        'Monkey must be 1 (Nilson) or 2 (Figaro).');
end

gateSuffix = sprintf('pTD%.3f_top%d_%.2f', pTDthr, TopKQuartets, RespThr);
gateSuffix = strrep(gateSuffix, '.', 'p');

tallPath = fullfile(cfg.matDir, tallFile);
resp3binPath = fullfile(cfg.matDir, resp3binFile);
respMoviePath = fullfile(cfg.matDir, respMovieFile);
outValuesFile = fullfile(cfg.resultsDir, ...
    sprintf('post_affine_delta_points_allbins_IT_%s_stim%03d_%s_%s.mat', ...
    char(monkeySuffix), ExampleStimulus, RespMovieTag, gateSuffix));
prepNoiseFile = fullfile(cfg.resultsDir, ...
    sprintf('knn_noise_signal_prep_IT_%s_stim%03d_%s_%s.mat', ...
    char(monkeySuffix), ExampleStimulus, RespMovieTag, gateSuffix));
outMovie = fullfile(cfg.resultsDir, ...
    sprintf('IT_attentiondiff_movie_postaffine_%s_K%d_%s_stim%03d_%s.mp4', ...
    char(monkeySuffix), PrepK, RespMovieTag, ExampleStimulus, gateSuffix));

if exist(outMovie, 'file') == 2 && ~ForceRender
    fprintf('Skipping IT attention movie because output already exists:\n%s\n', outMovie);
    fprintf('Set ForceRender = true to re-render it.\n');
    return;
end

assert(exist(tallPath, 'file') == 2, ...
    'Missing %s. Run Line_Stimuli_IT.m first.', tallPath);
assert(exist(resp3binPath, 'file') == 2, ...
    'Missing %s. Need the 3-bin response summary for IT normalization.', resp3binPath);
assert(exist(respMoviePath, 'file') == 2, ...
    'Missing %s. Need the %s response file for the time-resolved movie.', respMoviePath, RespMovieTag);

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
    'pThresh', pTDthr, ...
    'excludeOverlap', ExcludeOverlap, ...
    'MinObjectStim', MinObjectStim, ...
    'RespThr', RespThr, ...
    'TopKQuartets', TopKQuartets));

ATT = Sel.ATT;
keep = Sel.keep(:);
isSig = Sel.isSig(:);

fprintf('IT attention movie gate: kept %d / %d | significant %d / %d kept\n', ...
    nnz(keep), nIT, nnz(isSig), nnz(keep));
assert(any(isSig), 'No IT sites passed the attention-significant gate.');

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
SNRnorm = compute_snr_per_color_sites(R3dec, Tall_IT, siteRows, 'Verbose', false);

siteWeightsByIndex = [];
if UseSiteWeights
    assert(all(isfield(ATT, {'muT','muD','varT','varD'})), ...
        'ATT must contain muT/muD/varT/varD for site weighting.');
    dSite = double(ATT.muT(:) - ATT.muD(:));
    varSite = 0.5 * (double(ATT.varT(:)) + double(ATT.varD(:)));
    varSite(~isfinite(varSite) | varSite < 0) = 0;
    den = sqrt(varSite + max(GroupWeightLambda, 0));
    den(~isfinite(den) | den <= 0) = NaN;
    siteWeightsByIndex = abs(dSite) ./ den;

    if GroupWeightUseNmatch
        assert(all(isfield(ATT, {'wY','wP'})), ...
            'ATT must contain wY/wP when GroupWeightUseNmatch=true.');
        nMatch = double(ATT.wY(:) + ATT.wP(:));
        nMatch(~isfinite(nMatch) | nMatch < 0) = 0;
        siteWeightsByIndex = siteWeightsByIndex .* sqrt(nMatch);
    end

    siteWeightsByIndex(~isfinite(siteWeightsByIndex) | siteWeightsByIndex < 0) = 0;
    siteWeightsByIndex(~isSig) = 0;

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
        warning(['IT attention movie requested weighted KNN averaging but no positive ' ...
                 'site weights were obtained. Falling back to unweighted averaging.']);
        siteWeightsByIndex = [];
    else
        fprintf(['Fixed IT site weights ready: nPos=%d | min/median/max=%.6g / %.6g / %.6g | ' ...
                 'clip p%.2f\n'], ...
            numel(wPos), min(wPos), median(wPos), max(wPos), GroupWeightClipPct);
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

assert(isfield(Rmovie, 'timeWindows') && size(Rmovie.timeWindows, 2) == 2, ...
    '%s must contain R.timeWindows as [nBins x 2].', respMovieFile);
assert(size(Rmovie.timeWindows, 1) > 3, ...
    '%s must contain many short windows for movie rendering.', respMovieFile);
fprintf('Preparing IT attention movie with %d %s frames for monkey %s\n', ...
    size(Rmovie.timeWindows, 1), RespMovieTag, char(monkeySuffix));

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
        'sigSiteMask', isSig, ...
        'saveFile', outValuesFile, ...
        'verbose', true);
    fprintf('Saved IT post-affine values to: %s\n', outValuesFile);
else
    fprintf('Loaded IT post-affine values from: %s\n', outValuesFile);
end

PREP_NOISE = analyze_knn_noise_signal_thresholds( ...
    OUT_postAffine, Tall_IT, ...
    'KList', PrepKList, ...
    'preEndMs', PreEndMs, ...
    'postStartMs', PostStartMs, ...
    'preQuantilePct', PreQuantilePct, ...
    'alphaFloorPct', AlphaFloorPct, ...
    'cMaxPostPct', CMaxPostPct, ...
    'kRef', PrepK, ...
    'enforceK', true, ...
    'siteWeightsByIndex', siteWeightsByIndex, ...
    'makePlot', true, ...
    'verbose', true, ...
    'saveFile', prepNoiseFile);

Tprep = PREP_NOISE.summary;
rowK = find(Tprep.K == PrepK, 1, 'first');
assert(~isempty(rowK), ...
    'Requested prep K=%d not found in the noise/signal summary.', PrepK);

Kuse = double(Tprep.K(rowK));
thrUse = double(Tprep.thresholdPreQ(rowK));
assert(isfinite(thrUse) && thrUse > 0, ...
    'Invalid thresholdPreQ at K=%d in prep summary.', Kuse);
cMaxUse = double(Tprep.cMaxSuggest(rowK));
if ~isfinite(cMaxUse) || cMaxUse <= 0
    cMaxUse = [];
end

fprintf(['Rendering IT post-affine attention movie | K=%d | threshold=%.6g | ' ...
         'cMax=%s\n'], Kuse, thrUse, mat2str(cMaxUse));

MOV = make_post_affine_attention_movie( ...
    outMovie, ExampleStimulus, ALLCOORDS, RTAB384, OUT_postAffine, ...
    'K', Kuse, ...
    'siteWeightsByIndex', siteWeightsByIndex, ...
    'alphaFullAt', thrUse, ...
    'colorRedAt', thrUse, ...
    'cMaxFixed', cMaxUse, ...
    'markerSize', PlotMarkerSize, ...
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
    'enforceK', false, ...
    'verbose', true);

fprintf('Saved post-affine IT attention movie to: %s\n', MOV.outMovie);
