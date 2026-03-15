% Render a V4 attention-effect movie for the line-task stimulus set.
%
% This follows the active V1 post-affine movie path:
% - significant-site mask from the 3-bin late-window analysis
% - time-resolved bins from the 20 ms response file
% - post-affine KNN pooling with prep-calibrated threshold and color scale

Monkey = 1; % 1 = Nilson, 2 = Figaro
ExampleStimulus = 38; % canonical projection frame for the movie
TimeIdx3 = 3; % late window in the 3-bin summary
pTDthr = 0.05; % hard significance threshold for site inclusion
ExcludeOverlap = true; % match the V1 attention analysis
RespMovieTag = 'd20'; % use 20 ms bins for the time-resolved movie
FrameRate = 10; % movie frame rate
Quality = 95; % MPEG-4 writer quality
PlotMarkerSize = 8; % match V1 movie marker size
PlotAlpha = 0.30; % constant marker alpha above threshold
PlotBgColor = [0.5 0.5 0.5]; % background tone
PlotCLow = [0.50 0.50 0.50]; % low-value color, same as background
PlotCHigh = [0.85 0.05 0.05]; % red at threshold
HotScale = true; % gray->red->hot scale for strong values
ColorHotMaxFactor = 8.0; % allow strong values to move beyond red
HardAlphaCutoff = true; % zero alpha below threshold, constant alpha above
PrepK = 500; % K used for the final post-affine movie
PrepKList = [10 50 80 200 500 1000]; % calibration sweep, same as V1
PreEndMs = 0; % pre bins satisfy window end <= this time
PostStartMs = 300; % post bins satisfy window start >= this time
PreQuantilePct = 99.98; % threshold from pre-stim pooled absolute values
AlphaFloorPct = 80; % retained for parity with the prep summary
CMaxPostPct = 95; % color max suggestion from post-stim pooled values
UseSiteWeights = true; % use the same reliability-weighted KNN averaging as V1
GroupWeightLambda = 1e-6; % stabilizer in the reliability denominator
GroupWeightUseNmatch = true; % multiply weights by sqrt(wY + wP)
GroupWeightClipPct = 95; % cap extreme positive weights
ForceRender = false; % if false, skip when the output movie already exists

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
    error('Make_V4_attention_movie:InvalidMonkey', ...
        'Monkey must be 1 (Nilson) or 2 (Figaro).');
end

tallPath = fullfile(cfg.matDir, tallFile);
resp3binPath = fullfile(cfg.matDir, resp3binFile);
respMoviePath = fullfile(cfg.matDir, respMovieFile);
outValuesFile = fullfile(cfg.resultsDir, ...
    sprintf('post_affine_delta_points_allbins_V4_%s_stim%03d_%s.mat', ...
    char(monkeySuffix), ExampleStimulus, RespMovieTag));
prepNoiseFile = fullfile(cfg.resultsDir, ...
    sprintf('knn_noise_signal_prep_V4_%s_stim%03d_%s.mat', ...
    char(monkeySuffix), ExampleStimulus, RespMovieTag));
outMovie = fullfile(cfg.resultsDir, ...
    sprintf('V4_attentiondiff_movie_postaffine_%s_K%d_%s_stim%03d.mp4', ...
    char(monkeySuffix), PrepK, RespMovieTag, ExampleStimulus));

if exist(outMovie, 'file') == 2 && ~ForceRender
    fprintf('Skipping V4 attention movie because output already exists:\n%s\n', outMovie);
    fprintf('Set ForceRender = true to re-render it.\n');
    return;
end

assert(exist(tallPath, 'file') == 2, ...
    'Missing %s. Run Line_Stimuli_V4.m first.', tallPath);
assert(exist(resp3binPath, 'file') == 2, ...
    'Missing %s. Need the 3-bin response summary for the significance mask.', resp3binPath);
assert(exist(respMoviePath, 'file') == 2, ...
    'Missing %s. Need the %s response file for the time-resolved movie.', respMoviePath, RespMovieTag);

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
optsTD = struct('v1Sites', siteRows, 'timeIdx', TimeIdx3, ...
    'excludeOverlap', ExcludeOverlap, 'verbose', false, 'epsDen', 1e-6);
OUT3 = attention_modulation_V1_3bin(R3dec, Tall_V4, SNRnorm, optsTD);

isSig = isfinite(OUT3.pValueTD) & (OUT3.pValueTD < pTDthr);
fprintf('V4 attention movie mask: %d / %d sites with pTD < %.3f\n', ...
    nnz(isSig), nV4, pTDthr);
assert(any(isSig), 'No V4 sites passed pTD < %.3f.', pTDthr);

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
        warning(['V4 attention movie requested weighted KNN averaging but no positive ' ...
                 'site weights were obtained. Falling back to unweighted averaging.']);
        siteWeightsByIndex = [];
    else
        fprintf(['Fixed site weights ready: nPos=%d | min/median/max=%.6g / %.6g / %.6g | ' ...
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
fprintf('Preparing V4 attention movie with %d %s frames for monkey %s\n', ...
    size(Rmovie.timeWindows, 1), RespMovieTag, char(monkeySuffix));

OUT_postAffine = compute_projected_delta_points_allbins( ...
    ExampleStimulus, Tall_V4, ALLCOORDS, RTAB384, Rmovie, SNRnorm, ...
    'siteRange', siteRows, ...
    'excludeOverlap', ExcludeOverlap, ...
    'stimIdx', 1:384, ...
    'sigSiteMask', isSig, ...
    'saveFile', outValuesFile, ...
    'verbose', true);

PREP_NOISE = analyze_knn_noise_signal_thresholds( ...
    OUT_postAffine, Tall_V4, ...
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

fprintf(['Rendering V4 post-affine attention movie | K=%d | threshold=%.6g | ' ...
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

fprintf('Saved post-affine V4 attention movie to: %s\n', MOV.outMovie);
