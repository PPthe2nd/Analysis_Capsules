function OUT = GaussianOccupancy_Tuning_IT()
% GAUSSIANOCCUPANCY_TUNING_IT
% pRF-style occupancy search for IT using a large circular-Gaussian library.
%
% Step 1:
%   Search a library of circular Gaussians against the pooled non-grey
%   stimulus mask. For each site, select the single Gaussian that best
%   explains the signed stimulus-level response with a weighted
%   baseline+gain fit.
%
% Step 2:
%   Keep the winning Gaussian fixed and fit separate yellow and purple
%   weights:
%       response = baseline + wY * overlapYellow + wP * overlapPurple
%
% Responses are stimulus-level and normalized as:
%   (muStim - muSpont) / sdSpont

%% Settings
P = struct();
P.Monkey = 1;                  % 1 = Nilson, 2 = Figaro
P.imageSize = [1024 768];      % [W H], same convention as render_stim_with_masks2
P.sigmaPx = [15 22 32 47 69 101 149 219 322 475];
P.centerStepFrac = 0.5;        % center spacing = centerStepFrac * sigma
P.padFrac = 0.05;              % 5% padding on each side of occupied bbox
P.kernelCutoffSD = 4;          % Gaussian kernel support truncation
P.minStimuli = 80;             % minimum finite stimulus responses required per site
P.varFloorFrac = 1e-3;         % floor on stimulus-response variance relative to the median
P.sigAlpha = 0.05;             % per-site threshold for approximate p-values
P.makeSummaryFigures = true;
P.makeLibraryConceptFigure = true;
P.makeExampleFigures = true;
P.nExampleSites = 4;           % per window
P.saveResult = true;
P.forceRefit = false;
P.forceRebuildLibrary = false;
P.siteProgressEvery = 10;
P.bboxProgressEvery = 50;
P.overlapProgressEveryStim = 10;
P.vePlotFloorPct = -100;
P.colorDominanceRatio = 2;
P.minDominantPairs = 8;
P.usePrecisionWeightedPairs = true;
P.colorIndexDenomFloor = 0;    % match V1/V4: no extra denominator floor beyond requiring denom > 0
P.colorIndexPlotRange = [-0.6 0.6];
P.colorIndexPlotNBins = 24;

cfg = config();

%% Monkey-specific files
if P.Monkey == 1
    monkeySuffix = "N";
    tallFile = 'Tall_IT_lines_N.mat';
    resp3binFile = 'SNR_capsules_N_d12.mat';
elseif P.Monkey == 2
    monkeySuffix = "F";
    tallFile = 'Tall_IT_lines_F.mat';
    resp3binFile = 'SNR_capsules_F_d12.mat';
else
    error('GaussianOccupancy_Tuning_IT:InvalidMonkey', ...
        'P.Monkey must be 1 (Nilson) or 2 (Figaro).');
end

tallPath = fullfile(cfg.matDir, tallFile);
resp3binPath = fullfile(cfg.matDir, resp3binFile);
cachePath = fullfile(cfg.matDir, sprintf('GaussianOccupancy_Library_IT_%s.mat', char(monkeySuffix)));
outPath = fullfile(cfg.matDir, sprintf('GaussianOccupancy_Tuning_IT_%s.mat', char(monkeySuffix)));

assert(exist(tallPath, 'file') == 2, ...
    'Missing %s. Run Line_Stimuli_IT.m first.', tallPath);
assert(exist(resp3binPath, 'file') == 2, ...
    'Missing %s. Create the 3-bin response summary first.', resp3binPath);

useCachedOut = exist(outPath, 'file') == 2 && ~P.forceRefit;
needsSave = ~useCachedOut;
pairParamsChanged = false;
pairRequiredFields = {'weightedYellowDominant','weightedPurpleDominant', ...
    'weightedPairDiff','colorIndex','sumPairWeight','effPairN', ...
    'muYellowPooled','muPurplePooled','pValuePooled','zStatPooled'};

if useCachedOut
    fprintf('Loading cached IT Gaussian occupancy tuning from %s\n', outPath);
    S = load(outPath);
    assert(isfield(S, 'OUT') && isstruct(S.OUT), ...
        '%s must contain struct OUT.', outPath);
    OUT = S.OUT;
    assert(isfield(OUT, 'Library') && isfield(OUT, 'FitSpatialEarly') && ...
        isfield(OUT, 'FitSpatialLate') && isfield(OUT, 'FitColorEarly') && ...
        isfield(OUT, 'FitColorLate') && isfield(OUT, 'RFrange'), ...
        'Cached OUT is missing required fields.');
    Library = OUT.Library;
    RFrange = OUT.RFrange(:);
    nIT = numel(RFrange);
    signedStimEarly = OUT.signedStimEarly;
    signedStimLate = OUT.signedStimLate;
    varStimEarly = OUT.varStimEarly;
    varStimLate = OUT.varStimLate;
    FitSpatialEarly = OUT.FitSpatialEarly;
    FitSpatialLate = OUT.FitSpatialLate;
    FitColorEarly = OUT.FitColorEarly;
    FitColorLate = OUT.FitColorLate;
    if ~isfield(OUT, 'P') || ~isfield(OUT.P, 'colorDominanceRatio') || ...
            ~isfield(OUT.P, 'minDominantPairs') || ...
            ~isfield(OUT.P, 'usePrecisionWeightedPairs') || ...
            ~isfield(OUT.P, 'colorIndexDenomFloor') || ...
            OUT.P.colorDominanceRatio ~= P.colorDominanceRatio || ...
            OUT.P.minDominantPairs ~= P.minDominantPairs || ...
            OUT.P.usePrecisionWeightedPairs ~= P.usePrecisionWeightedPairs || ...
            OUT.P.colorIndexDenomFloor ~= P.colorIndexDenomFloor
        pairParamsChanged = true;
    end
    if isfield(OUT, 'PairColorEarly') && isfield(OUT, 'PairColorLate')
        PairColorEarly = OUT.PairColorEarly;
        PairColorLate = OUT.PairColorLate;
        if ~all(isfield(PairColorEarly, pairRequiredFields)) || ~all(isfield(PairColorLate, pairRequiredFields))
            pairParamsChanged = true;
        end
    end
else
%% Load geometry and responses
Sgeo = load(tallPath);
assert(isfield(Sgeo, 'Tall_IT') && isfield(Sgeo, 'ALLCOORDS') && isfield(Sgeo, 'RTAB384') && ...
    isfield(Sgeo, 'RFrange'), '%s must contain Tall_IT, ALLCOORDS, RTAB384, and RFrange.', tallFile);

Tall_IT = Sgeo.Tall_IT;
ALLCOORDS = Sgeo.ALLCOORDS;
RTAB384 = Sgeo.RTAB384;
RFrange = Sgeo.RFrange(:);
nIT = numel(RFrange);
siteRows = (1:nIT).';

stimNums = arrayfun(@(x) x.stimNum, Tall_IT(:));
[stimNumsSorted, ordStim] = sort(stimNums(:));
assert(all(stimNumsSorted(:).' == 1:numel(Tall_IT)), ...
    'Tall_IT.stimNum must cover 1..%d exactly.', numel(Tall_IT));
Tall_IT = Tall_IT(ordStim);
nStim = numel(Tall_IT);

R3_full = load_capsules_struct_exclusion_aware(resp3binPath, monkeySuffix, 'cfg', cfg);
R3 = R3_full;
R3.meanAct = R3_full.meanAct(RFrange, :, :);
R3.meanSqAct = R3_full.meanSqAct(RFrange, :, :);
if ismatrix(R3_full.nTrials) && size(R3_full.nTrials, 1) >= max(RFrange)
    R3.nTrials = R3_full.nTrials(RFrange, :);
else
    R3.nTrials = R3_full.nTrials;
end
nTrialsPairMat = double(R3.nTrials);
if isvector(nTrialsPairMat)
    nTrialsPairMat = repmat(nTrialsPairMat(:).', nIT, 1);
end
assert(size(R3.meanAct,1) == nIT, 'Localized IT response rows do not match Tall_IT.');
assert(size(R3.meanAct,2) == nStim, 'Localized IT response stimuli do not match Tall_IT.');

SNR = compute_snr_per_color_sites(R3, Tall_IT, siteRows, 'Verbose', false);
muSpont = SNR.muSpont(siteRows);
sdSpont = SNR.sdSpont(siteRows);
[signedStimEarly, signedStimLate, varStimEarly, varStimLate] = ...
    compute_signed_stim_responses_with_var(R3, muSpont, sdSpont);

%% Build or load overlap cache
[Library, overlapYellow, overlapPurple] = load_or_build_library_cache( ...
    cachePath, ALLCOORDS, RTAB384, nStim, P);

fprintf('Gaussian library size: %d templates\n', numel(Library.sigmaPx));

Xyellow = double(overlapYellow);
Xpurple = double(overlapPurple);
Xall = Xyellow + Xpurple;

%% Fit spatial search and color split
spatialFields = {'status','nObs','bestGaussianIdx','baseline','gain','r2TrainPct', ...
    'effectRangeObs','effectAbsPeakObs','fStatApprox','pValueApprox','dfModel', ...
    'dfError','rssNullWeighted','rssFullWeighted','centerX','centerY','sigmaPx'};
colorFields = {'status','nObs','bestGaussianIdx','baseline','wYellow','wPurple', ...
    'r2TrainPct','effectRangeObs','effectAbsPeakObs','deltaR2OverSpatial', ...
    'fStatApprox','pValueApprox','fStatAddColor','pValueAddColor','dfModel', ...
    'dfError','rssNullWeighted','rssSpatialWeighted','rssFullWeighted','centerX','centerY','sigmaPx'};

FitSpatialEarly = init_fit_struct(nIT, spatialFields);
FitSpatialLate = init_fit_struct(nIT, spatialFields);
FitColorEarly = init_fit_struct(nIT, colorFields);
FitColorLate = init_fit_struct(nIT, colorFields);

fprintf('Scanning Gaussian library for IT occupancy tuning (%d sites, %d stimuli, %d templates)\n', ...
    nIT, nStim, numel(Library.sigmaPx));
tFit = tic;

for iSite = 1:nIT
    yEarly = signedStimEarly(iSite, :).';
    yLate = signedStimLate(iSite, :).';
    vEarly = varStimEarly(iSite, :).';
    vLate = varStimLate(iSite, :).';

    FitSpatialEarly(iSite) = fit_best_gaussian_occupancy(Xall, yEarly, vEarly, Library, P);
    FitSpatialLate(iSite) = fit_best_gaussian_occupancy(Xall, yLate, vLate, Library, P);

    FitColorEarly(iSite) = fit_color_split_for_best_gaussian(Xyellow, Xpurple, ...
        yEarly, vEarly, FitSpatialEarly(iSite), Library, P);
    FitColorLate(iSite) = fit_color_split_for_best_gaussian(Xyellow, Xpurple, ...
        yLate, vLate, FitSpatialLate(iSite), Library, P);

    if P.siteProgressEvery > 0 && (iSite == 1 || mod(iSite, P.siteProgressEvery) == 0 || iSite == nIT)
        elapsedSec = toc(tFit);
        rate = iSite / max(elapsedSec, eps);
        etaSec = (nIT - iSite) / max(rate, eps);
        fprintf('  Site %d / %d | elapsed %.1fs | ETA %.1fs\n', ...
            iSite, nIT, elapsedSec, etaSec);
    end
end
end

pairFields = {'status','nPairs','meanYellowDominant','meanPurpleDominant','meanPairDiff', ...
    'weightedYellowDominant','weightedPurpleDominant','weightedPairDiff', ...
    'muYellowPooled','muPurplePooled','varYellowPooled','varPurplePooled','NYPooled','NPPooled', ...
    'zStatPooled','pValuePooled','colorIndex','colorIndexDenom','meanPairVar','sumPairWeight','effPairN','tStat','pValue','df', ...
    'bestGaussianIdx','centerX','centerY','sigmaPx'};
nStim = size(signedStimEarly, 2);
[pairsA, pairsB] = build_complementary_pairs(nStim);

if pairParamsChanged || ~exist('PairColorEarly', 'var') || ~exist('PairColorLate', 'var')
    if ~exist('Xyellow', 'var') || ~exist('Xpurple', 'var')
        if ~exist('ALLCOORDS', 'var') || ~exist('RTAB384', 'var')
            SgeoPair = load(tallPath, 'ALLCOORDS', 'RTAB384');
            assert(isfield(SgeoPair, 'ALLCOORDS') && isfield(SgeoPair, 'RTAB384'), ...
                '%s must contain ALLCOORDS and RTAB384.', tallPath);
            ALLCOORDS = SgeoPair.ALLCOORDS;
            RTAB384 = SgeoPair.RTAB384;
        end
        [LibraryPair, overlapYellow, overlapPurple] = load_or_build_library_cache( ...
            cachePath, ALLCOORDS, RTAB384, nStim, P);
        Xyellow = double(overlapYellow);
        Xpurple = double(overlapPurple);
        if ~exist('Library', 'var') || isempty(Library)
            Library = LibraryPair;
        end
    end
    if ~exist('nTrialsPairMat', 'var')
        Rpair = load_capsules_struct_exclusion_aware(resp3binPath, monkeySuffix, 'cfg', cfg);
        assert(isfield(Rpair, 'nTrials'), '%s must contain R.nTrials.', resp3binFile);
        if ismatrix(Rpair.nTrials) && size(Rpair.nTrials, 1) >= max(RFrange)
            nTrialsPairMat = double(Rpair.nTrials(RFrange, :));
        else
            nTrialsPairMat = repmat(double(Rpair.nTrials(:)).', nIT, 1);
        end
    end

    PairColorEarly = init_fit_struct(nIT, pairFields);
    PairColorLate = init_fit_struct(nIT, pairFields);
    fprintf('Computing paired complementary-color test within selected Gaussian (%d sites)\n', nIT);
    tPair = tic;
    for iSite = 1:nIT
        PairColorEarly(iSite) = paired_color_dominance_test(Xyellow, Xpurple, ...
            signedStimEarly(iSite, :).', varStimEarly(iSite, :).', FitSpatialEarly(iSite), ...
            Library, pairsA, pairsB, nTrialsPairMat(iSite, :).', P);
        PairColorLate(iSite) = paired_color_dominance_test(Xyellow, Xpurple, ...
            signedStimLate(iSite, :).', varStimLate(iSite, :).', FitSpatialLate(iSite), ...
            Library, pairsA, pairsB, nTrialsPairMat(iSite, :).', P);

        if P.siteProgressEvery > 0 && (iSite == 1 || mod(iSite, P.siteProgressEvery) == 0 || iSite == nIT)
            elapsedSec = toc(tPair);
            rate = iSite / max(elapsedSec, eps);
            etaSec = (nIT - iSite) / max(rate, eps);
            fprintf('  Paired test site %d / %d | elapsed %.1fs | ETA %.1fs\n', ...
                iSite, nIT, elapsedSec, etaSec);
        end
    end
    needsSave = true;
end

%% Summaries
veSpatialEarly = [FitSpatialEarly.r2TrainPct].';
veSpatialLate = [FitSpatialLate.r2TrainPct].';
veColorEarly = [FitColorEarly.r2TrainPct].';
veColorLate = [FitColorLate.r2TrainPct].';
effSpatialEarly = [FitSpatialEarly.effectRangeObs].';
effSpatialLate = [FitSpatialLate.effectRangeObs].';
deltaColorEarly = [FitColorEarly.deltaR2OverSpatial].';
deltaColorLate = [FitColorLate.deltaR2OverSpatial].';
wYellowEarly = [FitColorEarly.wYellow].';
wPurpleEarly = [FitColorEarly.wPurple].';
wYellowLate = [FitColorLate.wYellow].';
wPurpleLate = [FitColorLate.wPurple].';
wOccEarly = 0.5 * (wYellowEarly + wPurpleEarly);
wOccLate = 0.5 * (wYellowLate + wPurpleLate);
wColorEarly = 0.5 * (wYellowEarly - wPurpleEarly);
wColorLate = 0.5 * (wYellowLate - wPurpleLate);
wColorDiffEarly = [FitColorEarly.wYellow].' - [FitColorEarly.wPurple].';
wColorDiffLate = [FitColorLate.wYellow].' - [FitColorLate.wPurple].';
pSpatialEarly = [FitSpatialEarly.pValueApprox].';
pSpatialLate = [FitSpatialLate.pValueApprox].';
pColorAddEarly = [FitColorEarly.pValueAddColor].';
pColorAddLate = [FitColorLate.pValueAddColor].';
pairDiffEarly = [PairColorEarly.weightedPairDiff].';
pairDiffLate = [PairColorLate.weightedPairDiff].';
pairColorIdxEarly = [PairColorEarly.colorIndex].';
pairColorIdxLate = [PairColorLate.colorIndex].';
pairPEarly = [PairColorEarly.pValue].';
pairPLate = [PairColorLate.pValue].';
pooledPEarly = [PairColorEarly.pValuePooled].';
pooledPLate = [PairColorLate.pValuePooled].';

isSpatialTunedEarly = isfinite(pSpatialEarly) & (pSpatialEarly < P.sigAlpha);
isSpatialTunedLate = isfinite(pSpatialLate) & (pSpatialLate < P.sigAlpha);
isColorSplitEarly = isfinite(pColorAddEarly) & (pColorAddEarly < P.sigAlpha);
isColorSplitLate = isfinite(pColorAddLate) & (pColorAddLate < P.sigAlpha);
isPairColorEarly = isfinite(pairPEarly) & (pairPEarly < P.sigAlpha);
isPairColorLate = isfinite(pairPLate) & (pairPLate < P.sigAlpha);
isRFColorEarly = isfinite(pooledPEarly) & (pooledPEarly < P.sigAlpha);
isRFColorLate = isfinite(pooledPLate) & (pooledPLate < P.sigAlpha);

fprintf('Usable IT sites for occupancy search (early): %d / %d\n', nnz(isfinite(veSpatialEarly)), nIT);
fprintf('Usable IT sites for occupancy search (late):  %d / %d\n', nnz(isfinite(veSpatialLate)), nIT);
if any(isfinite(veSpatialEarly))
    fprintf('Median occupancy VE early: %.2f%%\n', median(veSpatialEarly(isfinite(veSpatialEarly))));
end
if any(isfinite(veSpatialLate))
    fprintf('Median occupancy VE late:  %.2f%%\n', median(veSpatialLate(isfinite(veSpatialLate))));
end
fprintf('Spatial occupancy tuned sites at p < %.3f | early=%d late=%d\n', ...
    P.sigAlpha, nnz(isSpatialTunedEarly), nnz(isSpatialTunedLate));
fprintf('Separate yellow/purple weights improve fit at p < %.3f | early=%d late=%d\n', ...
    P.sigAlpha, nnz(isColorSplitEarly), nnz(isColorSplitLate));
fprintf('Precision-weighted paired complementary-color effect at p < %.3f with RF dominance >= %.1fx | early=%d late=%d\n', ...
    P.sigAlpha, P.colorDominanceRatio, nnz(isPairColorEarly), nnz(isPairColorLate));
fprintf('Pooled RF color effect at p < %.3f with RF dominance >= %.1fx | early=%d late=%d\n', ...
    P.sigAlpha, P.colorDominanceRatio, nnz(isRFColorEarly), nnz(isRFColorLate));
fprintf('Finite RF color index with nonzero denominator | early=%d late=%d\n', ...
    nnz(isfinite(pairColorIdxEarly)), nnz(isfinite(pairColorIdxLate)));

%% Figures
if P.makeLibraryConceptFigure || P.makeExampleFigures
    if ~exist('ALLCOORDS', 'var') || ~exist('RTAB384', 'var')
        Sref = load(tallPath, 'ALLCOORDS', 'RTAB384');
        assert(isfield(Sref, 'ALLCOORDS') && isfield(Sref, 'RTAB384'), ...
            '%s must contain ALLCOORDS and RTAB384.', tallPath);
        ALLCOORDS = Sref.ALLCOORDS;
        RTAB384 = Sref.RTAB384;
    end
    if ~exist('nStim', 'var')
        nStim = size(signedStimEarly, 2);
    end
    unionMask = build_union_occupancy_mask(ALLCOORDS, RTAB384, nStim, Library.imageSize(1), Library.imageSize(2));
end

if P.makeSummaryFigures
    if P.makeLibraryConceptFigure
        make_library_concept_figure(Library, unionMask, nStim, ...
            sprintf('IT Gaussian occupancy library concept (%s)', char(monkeySuffix)));
    end
    plot_ve_histogram(veSpatialEarly, veSpatialLate, P.vePlotFloorPct, ...
        sprintf('IT occupancy search variance explained (%s)', char(monkeySuffix)));
    plot_effect_histogram(effSpatialEarly, effSpatialLate, ...
        sprintf('IT occupancy search effect size (%s)', char(monkeySuffix)));
    plot_best_sigma_histogram(FitSpatialEarly, FitSpatialLate, P.sigmaPx, ...
        sprintf('IT best Gaussian sigma (%s)', char(monkeySuffix)));
    plot_delta_color_histogram(deltaColorEarly, deltaColorLate, ...
        sprintf('IT color split improvement (%s)', char(monkeySuffix)));
    plot_color_weight_components(wOccEarly, wColorEarly, isColorSplitEarly, ...
        wOccLate, wColorLate, isColorSplitLate, P.sigAlpha, ...
        sprintf('IT occupancy and color-bias weights within selected Gaussian (%s)', char(monkeySuffix)));
    plot_paired_color_difference(pairDiffEarly, isPairColorEarly, pairDiffLate, isPairColorLate, ...
        P.colorDominanceRatio, sprintf('IT precision-weighted paired complementary-color effect (%s)', char(monkeySuffix)));
    plot_rf_color_index_histogram(pairColorIdxEarly, isRFColorEarly, pairColorIdxLate, isRFColorLate, ...
        P.colorIndexPlotRange, P.colorIndexPlotNBins, sprintf('IT RF color index (%s)', char(monkeySuffix)));
end

if P.makeExampleFigures
    make_spatial_example_figure(Library, FitSpatialEarly, veSpatialEarly, RFrange, unionMask, ...
        sprintf('IT occupancy search examples early (%s)', char(monkeySuffix)), P.nExampleSites);
    make_spatial_example_figure(Library, FitSpatialLate, veSpatialLate, RFrange, unionMask, ...
        sprintf('IT occupancy search examples late (%s)', char(monkeySuffix)), P.nExampleSites);
end

TspEarly = build_spatial_table(FitSpatialEarly, RFrange);
TspLate = build_spatial_table(FitSpatialLate, RFrange);
TcolEarly = build_color_table(FitColorEarly, FitSpatialEarly, RFrange);
TcolLate = build_color_table(FitColorLate, FitSpatialLate, RFrange);
TpairEarly = build_paired_color_table(PairColorEarly, RFrange);
TpairLate = build_paired_color_table(PairColorLate, RFrange);

disp('Top IT occupancy-search sites by variance explained (early):');
disp(TspEarly(1:min(10,height(TspEarly)), :));
disp('Top IT occupancy-search sites by variance explained (late):');
disp(TspLate(1:min(10,height(TspLate)), :));
disp('Top IT sites with the largest color split improvement (late):');
disp(TcolLate(1:min(10,height(TcolLate)), :));
disp('Top IT sites in precision-weighted paired complementary-color test (late):');
disp(TpairLate(1:min(10,height(TpairLate)), :));

%% Pack output
OUT.P = P;
OUT.monkeySuffix = monkeySuffix;
OUT.RFrange = RFrange;
OUT.Library = Library;
OUT.signedStimEarly = signedStimEarly;
OUT.signedStimLate = signedStimLate;
OUT.varStimEarly = varStimEarly;
OUT.varStimLate = varStimLate;
OUT.FitSpatialEarly = FitSpatialEarly;
OUT.FitSpatialLate = FitSpatialLate;
OUT.FitColorEarly = FitColorEarly;
OUT.FitColorLate = FitColorLate;
OUT.PairColorEarly = PairColorEarly;
OUT.PairColorLate = PairColorLate;
OUT.isSpatialTunedEarly = isSpatialTunedEarly;
OUT.isSpatialTunedLate = isSpatialTunedLate;
OUT.isColorSplitEarly = isColorSplitEarly;
OUT.isColorSplitLate = isColorSplitLate;
OUT.isPairColorEarly = isPairColorEarly;
OUT.isPairColorLate = isPairColorLate;
OUT.isRFColorEarly = isRFColorEarly;
OUT.isRFColorLate = isRFColorLate;
OUT.TableSpatialEarly = TspEarly;
OUT.TableSpatialLate = TspLate;
OUT.TableColorEarly = TcolEarly;
OUT.TableColorLate = TcolLate;
OUT.TablePairColorEarly = TpairEarly;
OUT.TablePairColorLate = TpairLate;
OUT.cachePath = cachePath;

if needsSave && P.saveResult
    save(outPath, 'OUT', '-v7.3');
    fprintf('Saved IT Gaussian occupancy tuning to %s\n', outPath);
end
end

function [Library, overlapYellow, overlapPurple] = load_or_build_library_cache(cachePath, ALLCOORDS, RTAB384, nStim, P)
needBuild = P.forceRebuildLibrary || exist(cachePath, 'file') ~= 2;
if ~needBuild
    S = load(cachePath);
    if isfield(S, 'Library') && isfield(S, 'overlapYellow') && isfield(S, 'overlapPurple') && ...
            isfield(S, 'CacheParams') && cache_compatible(S.CacheParams, P, nStim)
        fprintf('Loading cached Gaussian occupancy library from %s\n', cachePath);
        Library = S.Library;
        overlapYellow = S.overlapYellow;
        overlapPurple = S.overlapPurple;
        return;
    end
    needBuild = true;
end

[bbox, W, H] = compute_non_grey_bbox(ALLCOORDS, RTAB384, nStim, P);
Library = build_gaussian_library(bbox, W, H, P);
[overlapYellow, overlapPurple] = build_overlap_matrices(ALLCOORDS, RTAB384, nStim, Library, W, H, P);

CacheParams = struct();
CacheParams.imageSize = P.imageSize;
CacheParams.sigmaPx = P.sigmaPx;
CacheParams.centerStepFrac = P.centerStepFrac;
CacheParams.padFrac = P.padFrac;
CacheParams.kernelCutoffSD = P.kernelCutoffSD;
CacheParams.nStim = nStim;

save(cachePath, 'Library', 'overlapYellow', 'overlapPurple', 'CacheParams', '-v7.3');
fprintf('Saved Gaussian occupancy library cache to %s\n', cachePath);
end

function tf = cache_compatible(CacheParams, P, nStim)
tf = isequal(CacheParams.imageSize, P.imageSize) && ...
    isequal(CacheParams.sigmaPx(:), P.sigmaPx(:)) && ...
    isequal(CacheParams.centerStepFrac, P.centerStepFrac) && ...
    isequal(CacheParams.padFrac, P.padFrac) && ...
    isequal(CacheParams.kernelCutoffSD, P.kernelCutoffSD) && ...
    isequal(CacheParams.nStim, nStim);
end

function [bbox, W, H] = compute_non_grey_bbox(ALLCOORDS, RTAB384, nStim, P)
W = P.imageSize(1);
H = P.imageSize(2);
colMin = W;
colMax = 1;
rowMin = H;
rowMax = 1;

fprintf('Scanning stimulus masks for occupied bbox (%d stimuli)\n', nStim);
tScan = tic;

for stimNum = 1:nStim
    [~, masks] = render_stim_with_masks2(ALLCOORDS, RTAB384, stimNum, ...
        'ImageSize', [W H], 'DrawDots', false);
    maskAll = masks.yellowArm | masks.purple;
    cols = find(any(maskAll, 1));
    rows = find(any(maskAll, 2));
    if ~isempty(cols)
        colMin = min(colMin, cols(1));
        colMax = max(colMax, cols(end));
    end
    if ~isempty(rows)
        rowMin = min(rowMin, rows(1));
        rowMax = max(rowMax, rows(end));
    end

    if P.bboxProgressEvery > 0 && (stimNum == 1 || mod(stimNum, P.bboxProgressEvery) == 0 || stimNum == nStim)
        elapsedSec = toc(tScan);
        rate = stimNum / max(elapsedSec, eps);
        etaSec = (nStim - stimNum) / max(rate, eps);
        fprintf('  BBox stim %d / %d | elapsed %.1fs | ETA %.1fs\n', ...
            stimNum, nStim, elapsedSec, etaSec);
    end
end

if colMax < colMin || rowMax < rowMin
    error('GaussianOccupancy_Tuning_IT:EmptyBBox', 'Could not determine occupied bbox from stimulus masks.');
end

padCols = max(1, round(P.padFrac * (colMax - colMin + 1)));
padRows = max(1, round(P.padFrac * (rowMax - rowMin + 1)));
bbox = struct();
bbox.colMin = max(1, colMin - padCols);
bbox.colMax = min(W, colMax + padCols);
bbox.rowMin = max(1, rowMin - padRows);
bbox.rowMax = min(H, rowMax + padRows);
fprintf('Occupied bbox with padding: cols [%d %d], rows [%d %d]\n', ...
    bbox.colMin, bbox.colMax, bbox.rowMin, bbox.rowMax);
end

function Library = build_gaussian_library(bbox, W, H, P)
nSigma = numel(P.sigmaPx);
nTotal = 0;
countPerSigma = zeros(nSigma, 1);
colAll = [];
rowAll = [];
sigmaAll = [];
sigmaIdxAll = [];

fprintf('Building circular-Gaussian library\n');
for iSigma = 1:nSigma
    sigma = P.sigmaPx(iSigma);
    step = P.centerStepFrac * sigma;
    cols = unique(round(bbox.colMin:step:bbox.colMax));
    rows = unique(round(bbox.rowMin:step:bbox.rowMax));
    if isempty(cols)
        cols = round((bbox.colMin + bbox.colMax) / 2);
    end
    if isempty(rows)
        rows = round((bbox.rowMin + bbox.rowMax) / 2);
    end
    [cc, rr] = meshgrid(cols, rows);
    countPerSigma(iSigma) = numel(cc);
    nTotal = nTotal + countPerSigma(iSigma);
    fprintf('  sigma %.1f px | step %.1f | centers %d x %d = %d\n', ...
        sigma, step, numel(cols), numel(rows), countPerSigma(iSigma));
    colAll = [colAll; cc(:)]; %#ok<AGROW>
    rowAll = [rowAll; rr(:)]; %#ok<AGROW>
    sigmaAll = [sigmaAll; repmat(sigma, countPerSigma(iSigma), 1)]; %#ok<AGROW>
    sigmaIdxAll = [sigmaIdxAll; repmat(iSigma, countPerSigma(iSigma), 1)]; %#ok<AGROW>
end
fprintf('Total Gaussian templates: %d\n', nTotal);

Library = struct();
Library.col = int32(colAll);
Library.row = int32(rowAll);
Library.centerX = double(colAll) - W/2;
Library.centerY = H/2 - double(rowAll);
Library.sigmaPx = double(sigmaAll);
Library.sigmaIdx = int32(sigmaIdxAll);
Library.countPerSigma = countPerSigma;
Library.bbox = bbox;
Library.imageSize = [W H];
Library.sigmaListPx = P.sigmaPx(:);
end

function [overlapYellow, overlapPurple] = build_overlap_matrices(ALLCOORDS, RTAB384, nStim, Library, W, H, P)
nGauss = numel(Library.sigmaPx);
overlapYellow = zeros(nStim, nGauss, 'single');
overlapPurple = zeros(nStim, nGauss, 'single');
nSigma = numel(P.sigmaPx);

kernelList = cell(nSigma, 1);
gaussCols = cell(nSigma, 1);
gaussRows = cell(nSigma, 1);
gaussIdx = cell(nSigma, 1);
for iSigma = 1:nSigma
    kernelList{iSigma} = make_gaussian_kernel_1d(P.sigmaPx(iSigma), P.kernelCutoffSD, max(W, H));
    idx = find(Library.sigmaIdx == iSigma);
    gaussIdx{iSigma} = idx;
    gaussCols{iSigma} = double(Library.col(idx));
    gaussRows{iSigma} = double(Library.row(idx));
end

totalEval = nStim * nGauss;
fprintf('Building yellow/purple overlap matrices (%d stimuli x %d templates = %d overlap samples)\n', ...
    nStim, nGauss, totalEval);
tBuild = tic;

for stimNum = 1:nStim
    [~, masks] = render_stim_with_masks2(ALLCOORDS, RTAB384, stimNum, ...
        'ImageSize', [W H], 'DrawDots', false);
    maskY = double(masks.yellowArm);
    maskP = double(masks.purple);

    for iSigma = 1:nSigma
        g = kernelList{iSigma};
        blurY = conv2(g(:), g(:)', maskY, 'same');
        blurP = conv2(g(:), g(:)', maskP, 'same');
        idx = gaussIdx{iSigma};
        overlapYellow(stimNum, idx) = single(sample_image_at(blurY, gaussRows{iSigma}, gaussCols{iSigma}));
        overlapPurple(stimNum, idx) = single(sample_image_at(blurP, gaussRows{iSigma}, gaussCols{iSigma}));
    end

    if P.overlapProgressEveryStim > 0 && ...
            (stimNum == 1 || mod(stimNum, P.overlapProgressEveryStim) == 0 || stimNum == nStim)
        elapsedSec = toc(tBuild);
        processed = stimNum * nGauss;
        rate = processed / max(elapsedSec, eps);
        etaSec = (totalEval - processed) / max(rate, eps);
        fprintf('  Overlap stim %d / %d | gaussian evals %d / %d | elapsed %.1fs | ETA %.1fs\n', ...
            stimNum, nStim, processed, totalEval, elapsedSec, etaSec);
    end
end
end

function g = make_gaussian_kernel_1d(sigma, cutoffSD, maxDim)
radius = min(ceil(cutoffSD * sigma), maxDim);
x = -radius:radius;
g = exp(-0.5 * (x ./ sigma).^2);
g = g ./ sum(g);
end

function vals = sample_image_at(img, rows, cols)
rows = min(max(round(rows), 1), size(img,1));
cols = min(max(round(cols), 1), size(img,2));
vals = img(sub2ind(size(img), rows(:), cols(:))).';
end

function [signedEarly, signedLate, varEarly, varLate] = compute_signed_stim_responses_with_var(R3, muSpont, sdSpont)
[nSites, nStim, ~] = size(R3.meanAct);

if isvector(R3.nTrials)
    nTrialsMat = repmat(double(R3.nTrials(:))', nSites, 1);
else
    nTrialsMat = double(R3.nTrials);
end

rEarly = squeeze(R3.meanAct(:,:,2));
rLate = squeeze(R3.meanAct(:,:,3));
rSqEarly = squeeze(R3.meanSqAct(:,:,2));
rSqLate = squeeze(R3.meanSqAct(:,:,3));

varStimEarly = nan(nSites, nStim);
varStimLate = nan(nSites, nStim);

validEarly = isfinite(rEarly) & isfinite(rSqEarly) & isfinite(nTrialsMat) & (nTrialsMat > 1);
validLate = isfinite(rLate) & isfinite(rSqLate) & isfinite(nTrialsMat) & (nTrialsMat > 1);

sampleVarEarly = max(0, rSqEarly - rEarly.^2);
sampleVarLate = max(0, rSqLate - rLate.^2);

sampleVarEarly(validEarly) = sampleVarEarly(validEarly) .* ...
    (nTrialsMat(validEarly) ./ max(nTrialsMat(validEarly) - 1, 1));
sampleVarLate(validLate) = sampleVarLate(validLate) .* ...
    (nTrialsMat(validLate) ./ max(nTrialsMat(validLate) - 1, 1));

varStimEarly(validEarly) = sampleVarEarly(validEarly) ./ nTrialsMat(validEarly);
varStimLate(validLate) = sampleVarLate(validLate) ./ nTrialsMat(validLate);

signedEarly = bsxfun(@rdivide, bsxfun(@minus, rEarly, muSpont), sdSpont);
signedLate = bsxfun(@rdivide, bsxfun(@minus, rLate, muSpont), sdSpont);
varEarly = bsxfun(@rdivide, varStimEarly, sdSpont.^2);
varLate = bsxfun(@rdivide, varStimLate, sdSpont.^2);

badNoise = ~isfinite(sdSpont) | (sdSpont <= 0);
signedEarly(badNoise, :) = NaN;
signedLate(badNoise, :) = NaN;
varEarly(badNoise, :) = NaN;
varLate(badNoise, :) = NaN;
end

function F = init_fit_struct(nSites, fields)
tmp = struct();
for i = 1:numel(fields)
    tmp.(fields{i}) = NaN;
end
tmp.status = "not_fit";
F = repmat(tmp, nSites, 1);
end

function F = fit_best_gaussian_occupancy(Xall, y, varY, Library, P)
fields = {'status','nObs','bestGaussianIdx','baseline','gain','r2TrainPct', ...
    'effectRangeObs','effectAbsPeakObs','fStatApprox','pValueApprox','dfModel', ...
    'dfError','rssNullWeighted','rssFullWeighted','centerX','centerY','sigmaPx'};
F = struct();
for i = 1:numel(fields)
    F.(fields{i}) = NaN;
end
F.status = "not_fit";

valid = isfinite(y) & isfinite(varY) & (varY > 0);
if nnz(valid) < P.minStimuli
    F.status = "too_few_stimuli";
    F.nObs = nnz(valid);
    return;
end

yv = y(valid);
vv = varY(valid);
Xv = Xall(valid, :);
F.nObs = numel(yv);

varFloor = max(median(vv(isfinite(vv) & vv > 0)) * P.varFloorFrac, 1e-6);
w = 1 ./ max(vv, varFloor);
S0 = sum(w);
Sy = sum(w .* yv);
yBar = Sy / S0;
rssNull = sum(w .* (yv - yBar).^2);

Sx = w.' * Xv;
Sxx = w.' * (Xv .* Xv);
Sxy = (w .* yv).' * Xv;

denom = Sxx - (Sx.^2) ./ S0;
numer = Sxy - (Sx .* Sy) ./ S0;
ok = denom > 0 & isfinite(denom) & isfinite(numer);
if ~any(ok)
    F.status = "no_valid_template";
    return;
end

gain = nan(size(denom));
gain(ok) = numer(ok) ./ denom(ok);
rssFull = inf(size(denom));
rssFull(ok) = rssNull - (numer(ok).^2) ./ denom(ok);
rssFull = max(rssFull, 0);

[bestRss, bestIdx] = min(rssFull); %#ok<ASGLU>
if ~isfinite(bestIdx) || ~isfinite(rssFull(bestIdx))
    F.status = "template_search_failed";
    return;
end

xBest = Xv(:, bestIdx);
gainBest = gain(bestIdx);
baselineBest = yBar - gainBest * (Sx(bestIdx) / S0);
yHat = baselineBest + gainBest .* xBest;

F.status = "ok";
F.bestGaussianIdx = bestIdx;
F.baseline = baselineBest;
F.gain = gainBest;
if rssNull > 0
    F.r2TrainPct = 100 * (1 - rssFull(bestIdx) / rssNull);
end
F.effectRangeObs = max(yHat) - min(yHat);
F.effectAbsPeakObs = max(abs(yHat));
F.dfModel = 1;
F.dfError = numel(yv) - 2;
F.rssNullWeighted = rssNull;
F.rssFullWeighted = rssFull(bestIdx);
F.centerX = Library.centerX(bestIdx);
F.centerY = Library.centerY(bestIdx);
F.sigmaPx = Library.sigmaPx(bestIdx);

if F.dfError > 0 && rssNull > rssFull(bestIdx)
    F.fStatApprox = ((rssNull - rssFull(bestIdx)) / F.dfModel) / (rssFull(bestIdx) / F.dfError);
    F.pValueApprox = 1 - fcdf_local(F.fStatApprox, F.dfModel, F.dfError);
else
    F.fStatApprox = 0;
    F.pValueApprox = 1;
end
end

function F = fit_color_split_for_best_gaussian(Xyellow, Xpurple, y, varY, Fspatial, Library, P)
fields = {'status','nObs','bestGaussianIdx','baseline','wYellow','wPurple', ...
    'r2TrainPct','effectRangeObs','effectAbsPeakObs','deltaR2OverSpatial', ...
    'fStatApprox','pValueApprox','fStatAddColor','pValueAddColor','dfModel', ...
    'dfError','rssNullWeighted','rssSpatialWeighted','rssFullWeighted','centerX','centerY','sigmaPx'};
F = struct();
for i = 1:numel(fields)
    F.(fields{i}) = NaN;
end
F.status = "not_fit";

if string(Fspatial.status) ~= "ok" || ~isfinite(Fspatial.bestGaussianIdx)
    F.status = "no_spatial_model";
    return;
end

bestIdx = Fspatial.bestGaussianIdx;
valid = isfinite(y) & isfinite(varY) & (varY > 0);
if nnz(valid) < P.minStimuli
    F.status = "too_few_stimuli";
    F.nObs = nnz(valid);
    return;
end

yv = y(valid);
vv = varY(valid);
xY = Xyellow(valid, bestIdx);
xP = Xpurple(valid, bestIdx);
F.nObs = numel(yv);

varFloor = max(median(vv(isfinite(vv) & vv > 0)) * P.varFloorFrac, 1e-6);
w = 1 ./ max(vv, varFloor);
sqrtW = sqrt(w);
X = [ones(numel(yv),1), xY, xP];
Xw = bsxfun(@times, X, sqrtW);
yw = sqrtW .* yv;

beta = Xw \ yw;
if any(~isfinite(beta))
    F.status = "fit_failed";
    return;
end

yHat = X * beta;
wMean = sum(w .* yv) / sum(w);
rssNull = sum(w .* (yv - wMean).^2);
rssFull = sum(w .* (yv - yHat).^2);
rssSpatial = Fspatial.rssFullWeighted;

F.status = "ok";
F.bestGaussianIdx = bestIdx;
F.baseline = beta(1);
F.wYellow = beta(2);
F.wPurple = beta(3);
if rssNull > 0
    F.r2TrainPct = 100 * (1 - rssFull / rssNull);
    F.deltaR2OverSpatial = 100 * (rssSpatial - rssFull) / rssNull;
end
F.effectRangeObs = max(yHat) - min(yHat);
F.effectAbsPeakObs = max(abs(yHat));
F.dfModel = 2;
F.dfError = numel(yv) - 3;
F.rssNullWeighted = rssNull;
F.rssSpatialWeighted = rssSpatial;
F.rssFullWeighted = rssFull;
F.centerX = Library.centerX(bestIdx);
F.centerY = Library.centerY(bestIdx);
F.sigmaPx = Library.sigmaPx(bestIdx);

if F.dfError > 0 && rssNull > rssFull
    F.fStatApprox = ((rssNull - rssFull) / F.dfModel) / (rssFull / F.dfError);
    F.pValueApprox = 1 - fcdf_local(F.fStatApprox, F.dfModel, F.dfError);
else
    F.fStatApprox = 0;
    F.pValueApprox = 1;
end

if F.dfError > 0 && rssSpatial > rssFull
    F.fStatAddColor = ((rssSpatial - rssFull) / 1) / (rssFull / F.dfError);
    F.pValueAddColor = 1 - fcdf_local(F.fStatAddColor, 1, F.dfError);
else
    F.fStatAddColor = 0;
    F.pValueAddColor = 1;
end
end

function F = paired_color_dominance_test(Xyellow, Xpurple, y, varY, Fspatial, Library, pairsA, pairsB, nTrialsStim, P)
fields = {'status','nPairs','meanYellowDominant','meanPurpleDominant','meanPairDiff', ...
    'weightedYellowDominant','weightedPurpleDominant','weightedPairDiff', ...
    'muYellowPooled','muPurplePooled','varYellowPooled','varPurplePooled','NYPooled','NPPooled', ...
    'zStatPooled','pValuePooled','colorIndex','colorIndexDenom','meanPairVar','sumPairWeight','effPairN','tStat','pValue','df', ...
    'bestGaussianIdx','centerX','centerY','sigmaPx'};
F = struct();
for i = 1:numel(fields)
    F.(fields{i}) = NaN;
end
F.status = "not_fit";

if string(Fspatial.status) ~= "ok" || ~isfinite(Fspatial.bestGaussianIdx)
    F.status = "no_spatial_model";
    return;
end

bestIdx = Fspatial.bestGaussianIdx;
xY = Xyellow(:, bestIdx);
xP = Xpurple(:, bestIdx);

pairDiff = nan(numel(pairsA), 1);
pairVar = nan(numel(pairsA), 1);
respYellow = nan(numel(pairsA), 1);
respPurple = nan(numel(pairsA), 1);
nKeepWeighted = 0;
nKeepPooled = 0;
sumY = 0; sumY2 = 0; NY = 0;
sumP = 0; sumP2 = 0; NP = 0;

for k = 1:numel(pairsA)
    a = pairsA(k);
    b = pairsB(k);
    if any(~isfinite([y(a), y(b), xY(a), xY(b), xP(a), xP(b)]))
        continue;
    end
    na = nTrialsStim(a);
    nb = nTrialsStim(b);
    nEff = min(na, nb);
    if ~(isfinite(nEff) && nEff > 0)
        continue;
    end

    isYellowA = is_color_dominant(xY(a), xP(a), P.colorDominanceRatio);
    isPurpleA = is_color_dominant(xP(a), xY(a), P.colorDominanceRatio);
    isYellowB = is_color_dominant(xY(b), xP(b), P.colorDominanceRatio);
    isPurpleB = is_color_dominant(xP(b), xY(b), P.colorDominanceRatio);

    if isYellowA && isPurpleB
        yStim = a;
        pStim = b;
    elseif isYellowB && isPurpleA
        yStim = b;
        pStim = a;
    else
        continue;
    end

    meanSqY = recover_normalized_mean_square(y(yStim), varY(yStim), nTrialsStim(yStim));
    meanSqP = recover_normalized_mean_square(y(pStim), varY(pStim), nTrialsStim(pStim));
    if isfinite(meanSqY)
        sumY = sumY + nEff * y(yStim);
        sumY2 = sumY2 + nEff * meanSqY;
        NY = NY + nEff;
    end
    if isfinite(meanSqP)
        sumP = sumP + nEff * y(pStim);
        sumP2 = sumP2 + nEff * meanSqP;
        NP = NP + nEff;
    end
    if isfinite(meanSqY) && isfinite(meanSqP)
        nKeepPooled = nKeepPooled + 1;
    end

    if ~(isfinite(varY(yStim)) && isfinite(varY(pStim)))
        continue;
    end
    nKeepWeighted = nKeepWeighted + 1;
    respYellow(nKeepWeighted) = y(yStim);
    respPurple(nKeepWeighted) = y(pStim);
    pairDiff(nKeepWeighted) = y(yStim) - y(pStim);
    pairVar(nKeepWeighted) = max(varY(yStim), 0) + max(varY(pStim), 0);
end

F.nPairs = nKeepPooled;
F.bestGaussianIdx = bestIdx;
F.centerX = Library.centerX(bestIdx);
F.centerY = Library.centerY(bestIdx);
F.sigmaPx = Library.sigmaPx(bestIdx);

if nKeepWeighted < P.minDominantPairs
    F.status = "too_few_pairs";
else
    F.status = "ok";
end

[muYpool, varYpool, muPpool, varPpool, NYpool, NPpool, cind, ~, zPool, pPool] = ...
    compute_color_metrics(sumY, sumY2, NY, sumP, sumP2, NP);
F.muYellowPooled = muYpool;
F.muPurplePooled = muPpool;
F.varYellowPooled = varYpool;
F.varPurplePooled = varPpool;
F.NYPooled = NYpool;
F.NPPooled = NPpool;
F.zStatPooled = zPool;
F.pValuePooled = pPool;
denomCI = abs((muYpool + muPpool) / 2);
F.colorIndexDenom = denomCI;
if denomCI >= P.colorIndexDenomFloor && isfinite(cind)
    F.colorIndex = cind;
end
if nKeepWeighted < P.minDominantPairs
    return;
end

pairDiff = pairDiff(1:nKeepWeighted);
pairVar = pairVar(1:nKeepWeighted);
respYellow = respYellow(1:nKeepWeighted);
respPurple = respPurple(1:nKeepWeighted);
posPairVar = pairVar(isfinite(pairVar) & pairVar > 0);
if isempty(posPairVar)
    pairVarFloor = 1e-6;
else
    pairVarFloor = max(median(posPairVar) * P.varFloorFrac, 1e-6);
end
wPair = 1 ./ max(pairVar, pairVarFloor);
sumW = sum(wPair);
meanDiff = mean(pairDiff);
weightedDiff = sum(wPair .* pairDiff) / sumW;
F.meanYellowDominant = mean(respYellow);
F.meanPurpleDominant = mean(respPurple);
F.meanPairDiff = meanDiff;
F.weightedYellowDominant = sum(wPair .* respYellow) / sumW;
F.weightedPurpleDominant = sum(wPair .* respPurple) / sumW;
F.weightedPairDiff = weightedDiff;
F.meanPairVar = mean(pairVar);
F.sumPairWeight = sumW;
F.effPairN = (sumW^2) / sum(wPair.^2);
F.df = nKeepWeighted - 1;

if nKeepWeighted >= 2
    resid = pairDiff - weightedDiff;
    rssWeighted = sum(wPair .* (resid .^ 2));
    if isfinite(rssWeighted) && F.df > 0
        sigma2Hat = rssWeighted / F.df;
        seWeighted = sqrt(sigma2Hat / sumW);
    else
        sigma2Hat = NaN; %#ok<NASGU>
        seWeighted = NaN;
    end
    if isfinite(seWeighted) && seWeighted > 0
        F.tStat = weightedDiff / seWeighted;
        F.pValue = 1 - fcdf_local(F.tStat.^2, 1, F.df);
    elseif weightedDiff ~= 0
        F.tStat = sign(weightedDiff) * Inf;
        F.pValue = 0;
    else
        F.tStat = 0;
        F.pValue = 1;
    end
else
    F.tStat = 0;
    F.pValue = 1;
end

function meanSq = recover_normalized_mean_square(muNorm, varMeanNorm, nOrig)
meanSq = NaN;
if ~(isfinite(muNorm) && isfinite(nOrig) && nOrig > 0)
    return;
end
if nOrig <= 1
    meanSq = muNorm.^2;
    return;
end
if ~isfinite(varMeanNorm) || varMeanNorm < 0
    return;
end
sampleVarBessel = varMeanNorm .* nOrig;
sampleVarBiased = sampleVarBessel .* ((nOrig - 1) ./ nOrig);
meanSq = sampleVarBiased + muNorm.^2;
end

function [muY,varY,muP,varP,NY,NP,colorIndex,dprime,z,p] = compute_color_metrics(sumY,sumY2,NY,sumP,sumP2,NP)
muY = NaN; varY = NaN; muP = NaN; varP = NaN;
colorIndex = NaN; dprime = NaN; z = NaN; p = NaN;

if NY <= 1 || NP <= 1
    return;
end

muY = sumY / NY;
muP = sumP / NP;

Ex2Y = sumY2 / NY;
Ex2P = sumP2 / NP;

varY = max(0, Ex2Y - muY^2);
varP = max(0, Ex2P - muP^2);

varY = varY * (NY / max(NY - 1, 1));
varP = varP * (NP / max(NP - 1, 1));

denom = abs((muY + muP) / 2);
if denom > 0
    colorIndex = (muY - muP) / denom;
end

denomDP = sqrt(0.5 * (varY + varP));
if denomDP > 0
    dprime = (muY - muP) / denomDP;
end

denomZ = sqrt(varY / NY + varP / NP);
if denomZ > 0
    z = (muY - muP) / denomZ;
    p = 2 * normcdf(-abs(z), 0, 1);
end
end
end

function plot_ve_histogram(veEarly, veLate, floorPct, ttl)
vePlotEarly = veEarly;
vePlotLate = veLate;
vePlotEarly(isfinite(vePlotEarly) & (vePlotEarly < floorPct)) = floorPct;
vePlotLate(isfinite(vePlotLate) & (vePlotLate < floorPct)) = floorPct;

figure('Color', 'w');
histogram(vePlotEarly(isfinite(vePlotEarly)), 30, 'FaceColor', [0.25 0.45 0.85], 'EdgeColor', 'none');
hold on;
histogram(vePlotLate(isfinite(vePlotLate)), 30, 'FaceColor', [0.85 0.35 0.25], 'EdgeColor', 'none');
xlabel('Variance explained (%)');
ylabel('N sites');
title(ttl);
legend('Early', 'Late');
grid on;
end

function plot_effect_histogram(effEarly, effLate, ttl)
figure('Color', 'w');
histogram(effEarly(isfinite(effEarly)), 30, 'FaceColor', [0.25 0.45 0.85], 'EdgeColor', 'none');
hold on;
histogram(effLate(isfinite(effLate)), 30, 'FaceColor', [0.85 0.35 0.25], 'EdgeColor', 'none');
xlabel('Fitted modulation range (spont. SD units)');
ylabel('N sites');
title(ttl);
legend('Early', 'Late');
grid on;
end

function plot_best_sigma_histogram(FitEarly, FitLate, sigmaList, ttl)
sigEarly = [FitEarly.sigmaPx].';
sigLate = [FitLate.sigmaPx].';
cntEarly = zeros(numel(sigmaList), 1);
cntLate = zeros(numel(sigmaList), 1);
for i = 1:numel(sigmaList)
    cntEarly(i) = nnz(isfinite(sigEarly) & sigEarly == sigmaList(i));
    cntLate(i) = nnz(isfinite(sigLate) & sigLate == sigmaList(i));
end

figure('Color', 'w');
barX = 1:numel(sigmaList);
hBar = bar(barX, [cntEarly cntLate], 'grouped', 'BarWidth', 0.9);
hBar(1).FaceColor = [0.25 0.45 0.85];
hBar(1).EdgeColor = 'none';
hBar(2).FaceColor = [0.85 0.35 0.25];
hBar(2).EdgeColor = 'none';
xlabel('Best sigma (px)');
ylabel('N sites');
title(ttl);
set(gca, 'XTick', barX, 'XTickLabel', arrayfun(@(s) sprintf('%.0f', s), sigmaList, 'UniformOutput', false));
xlim([0.5, numel(sigmaList) + 0.5]);
legend('Early', 'Late');
grid on;
set(gca, 'XGrid', 'off', 'Box', 'on');
end

function plot_delta_color_histogram(deltaEarly, deltaLate, ttl)
figure('Color', 'w');
histogram(deltaEarly(isfinite(deltaEarly)), 30, 'FaceColor', [0.25 0.45 0.85], 'EdgeColor', 'none');
hold on;
histogram(deltaLate(isfinite(deltaLate)), 30, 'FaceColor', [0.85 0.35 0.25], 'EdgeColor', 'none');
xline(0, 'k-');
xlabel('Color-split improvement in variance explained (%)');
ylabel('N sites');
title(ttl);
legend('Early', 'Late');
grid on;
end

function plot_color_weight_components(wOccEarly, wColorEarly, sigEarly, wOccLate, wColorLate, sigLate, alphaLevel, ttl)
figure('Color', 'w');
useTiled = exist('tiledlayout', 'file') == 2;
if useTiled
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
end

wOccRange = choose_symmetric_plot_range([wOccEarly(:); wOccLate(:)], 0.98);
wColorRange = choose_symmetric_plot_range([wColorEarly(:); wColorLate(:)], 0.98);
nBins = 24;

if useTiled, nexttile; else, subplot(2, 2, 1); end
plot_overflow_histogram(wOccEarly, sigEarly, wOccRange, nBins);
xline(0, 'k-');
xlabel('Overall occupancy weight  wOcc = (wYellow + wPurple)/2');
ylabel('N sites');
title(sprintf('Early | significant split p < %.2f', alphaLevel));
legend('All sites', 'Significant color split');
grid on;

if useTiled, nexttile; else, subplot(2, 2, 2); end
plot_overflow_histogram(wOccLate, sigLate, wOccRange, nBins);
xline(0, 'k-');
xlabel('Overall occupancy weight  wOcc = (wYellow + wPurple)/2');
ylabel('N sites');
title(sprintf('Late | significant split p < %.2f', alphaLevel));
legend('All sites', 'Significant color split');
grid on;

if useTiled, nexttile; else, subplot(2, 2, 3); end
plot_overflow_histogram(wColorEarly, sigEarly, wColorRange, nBins);
xline(0, 'k-');
xlabel('Color-bias weight  wColor = (wYellow - wPurple)/2');
ylabel('N sites');
title('Early | positive = yellow > purple');
legend('All sites', 'Significant color split');
grid on;

if useTiled, nexttile; else, subplot(2, 2, 4); end
plot_overflow_histogram(wColorLate, sigLate, wColorRange, nBins);
xline(0, 'k-');
xlabel('Color-bias weight  wColor = (wYellow - wPurple)/2');
ylabel('N sites');
title('Late | positive = yellow > purple');
legend('All sites', 'Significant color split');
grid on;

if exist('sgtitle', 'file') == 2
    sgtitle(ttl);
end
end

function plot_overflow_histogram(vals, sigMask, plotRange, nCentralBins)
vals = vals(:);
sigMask = sigMask(:);
valid = isfinite(vals);
vals = vals(valid);
sigMask = sigMask(valid) & isfinite(sigMask(valid));

lo = plotRange(1);
hi = plotRange(2);
if ~(isfinite(lo) && isfinite(hi) && lo < hi)
    error('plot_overflow_histogram:InvalidRange', 'plotRange must be [lo hi] with lo < hi.');
end
nCentralBins = max(round(nCentralBins), 4);
step = (hi - lo) / nCentralBins;
edges = [-Inf, lo:step:hi, Inf];
centers = [lo - step/2, lo + step*(0.5:(nCentralBins-0.5)), hi + step/2];

countAll = histcounts(vals, edges);
countSig = histcounts(vals(sigMask), edges);

bar(centers, countAll, 1.0, 'FaceColor', [0.82 0.82 0.82], 'EdgeColor', 'none');
hold on;
bar(centers, countSig, 1.0, 'FaceColor', [0.80 0.10 0.10], 'EdgeColor', 'none');

xlim([centers(1) - step/2, centers(end) + step/2]);
xticks([centers(1), 0, centers(end)]);
xticklabels({sprintf('<%.2g', lo), '0', sprintf('>%.2g', hi)});
end

function plotRange = choose_symmetric_plot_range(vals, pctKeep)
vals = vals(isfinite(vals));
if isempty(vals)
    plotRange = [-1 1];
    return;
end
maxAbs = prctile(abs(vals), pctKeep * 100);
if ~isfinite(maxAbs) || maxAbs <= 0
    maxAbs = max(abs(vals));
end
if ~isfinite(maxAbs) || maxAbs <= 0
    maxAbs = 1;
end
plotRange = 1.05 * [-maxAbs maxAbs];
end

function plot_paired_color_difference(diffEarly, sigEarly, diffLate, sigLate, ratioThr, ttl)
figure('Color', 'w');
useTiled = exist('tiledlayout', 'file') == 2;
if useTiled
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
end

if useTiled, nexttile; else, subplot(1, 2, 1); end
histogram(diffEarly(isfinite(diffEarly)), 30, 'FaceColor', [0.82 0.82 0.82], 'EdgeColor', 'none');
hold on;
histogram(diffEarly(isfinite(diffEarly) & sigEarly), 30, 'FaceColor', [0.80 0.10 0.10], 'EdgeColor', 'none');
xline(0, 'k-');
xlabel('Precision-weighted paired difference (yellow-dom - purple-dom)');
ylabel('N sites');
title(sprintf('Early | RF dominance >= %.1fx', ratioThr));
legend('All sites', 'Significant paired effect');
grid on;

if useTiled, nexttile; else, subplot(1, 2, 2); end
histogram(diffLate(isfinite(diffLate)), 30, 'FaceColor', [0.82 0.82 0.82], 'EdgeColor', 'none');
hold on;
histogram(diffLate(isfinite(diffLate) & sigLate), 30, 'FaceColor', [0.80 0.10 0.10], 'EdgeColor', 'none');
xline(0, 'k-');
xlabel('Precision-weighted paired difference (yellow-dom - purple-dom)');
ylabel('N sites');
title(sprintf('Late | RF dominance >= %.1fx', ratioThr));
legend('All sites', 'Significant paired effect');
grid on;

if exist('sgtitle', 'file') == 2
    sgtitle(ttl);
end
end

function plot_rf_color_index_histogram(ciEarly, sigEarly, ciLate, sigLate, plotRange, nCentralBins, ttl)
figure('Color', 'w');
useTiled = exist('tiledlayout', 'file') == 2;
if useTiled
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
end

if useTiled, nexttile; else, subplot(1, 2, 1); end
plot_overflow_histogram(ciEarly, sigEarly, plotRange, nCentralBins);
xline(0, 'k-');
xlabel('RF color index (yellow - purple)');
ylabel('N sites');
title('Early');
legend('All sites', 'Significant pooled effect');
grid on;

if useTiled, nexttile; else, subplot(1, 2, 2); end
plot_overflow_histogram(ciLate, sigLate, plotRange, nCentralBins);
xline(0, 'k-');
xlabel('RF color index (yellow - purple)');
ylabel('N sites');
title('Late');
legend('All sites', 'Significant pooled effect');
grid on;

if exist('sgtitle', 'file') == 2
    sgtitle(ttl);
end
end

function make_library_concept_figure(Library, unionMask, nStim, figName)
W = Library.imageSize(1);
H = Library.imageSize(2);
sigmaList = double(Library.sigmaListPx(:));
[xOcc, yOcc] = mask_to_xy(unionMask, W, H, 30000);
bboxPos = bbox_rect_position(Library.bbox, W, H);
xSpan = bboxPos(3);
ySpan = bboxPos(4);
xMargin = max(20, 0.08 * xSpan);
yMargin = max(20, 0.08 * ySpan);
xLim = [bboxPos(1) - xMargin, bboxPos(1) + bboxPos(3) + xMargin];
yLim = [bboxPos(2) - yMargin, bboxPos(2) + bboxPos(4) + yMargin];

figure('Color', 'w');
useTiled = exist('tiledlayout', 'file') == 2;
if useTiled
    tl = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    ax1 = nexttile(tl);
    ax2 = nexttile(tl);
    ax3 = nexttile(tl);
else
    ax1 = subplot(1, 3, 1);
    ax2 = subplot(1, 3, 2);
    ax3 = subplot(1, 3, 3);
end

plot(ax1, xOcc, yOcc, '.', 'Color', [0.84 0.84 0.84], 'MarkerSize', 3);
hold(ax1, 'on');
draw_library_bbox(ax1, bboxPos);
title(ax1, sprintf('Union occupancy\n%d stimuli', nStim));
xlabel(ax1, 'x');
ylabel(ax1, 'y');
style_library_axes(ax1, xLim, yLim);

plot(ax2, xOcc, yOcc, '.', 'Color', [0.9 0.9 0.9], 'MarkerSize', 3);
hold(ax2, 'on');
scatter(ax2, Library.centerX, Library.centerY, 10, Library.sigmaPx, ...
    'filled', 'MarkerEdgeColor', 'none');
draw_library_bbox(ax2, bboxPos);
title(ax2, sprintf('Candidate centers\n%d templates', numel(Library.sigmaPx)));
xlabel(ax2, 'x');
style_library_axes(ax2, xLim, yLim);
cb = colorbar(ax2);
ylabel(cb, 'sigma (px)');
cb.Ticks = sigmaList(:).';
cb.TickLabels = arrayfun(@(s) sprintf('%.0f', s), sigmaList, 'UniformOutput', false);

plot(ax3, xOcc, yOcc, '.', 'Color', [0.9 0.9 0.9], 'MarkerSize', 3);
hold(ax3, 'on');
draw_library_bbox(ax3, bboxPos);
[x0, y0] = representative_library_center(Library, bboxPos);
sigmaShowIdx = unique(round(linspace(1, numel(sigmaList), min(4, numel(sigmaList)))));
sigmaShow = sigmaList(sigmaShowIdx);
cmap = parula(max(numel(sigmaList), 2));
hCirc = zeros(numel(sigmaShow), 1);
for i = 1:numel(sigmaShow)
    idxSig = find(abs(sigmaList - sigmaShow(i)) < 1e-9, 1, 'first');
    if isempty(idxSig)
        idxSig = sigmaShowIdx(i);
    end
    hCirc(i) = plot_sigma_circle(ax3, x0, y0, sigmaShow(i), cmap(idxSig, :), 1.5);
end
plot(ax3, x0, y0, 'k+', 'MarkerSize', 8, 'LineWidth', 1.2);
title(ax3, sprintf('Representative footprints\ncenter at (%.0f, %.0f)', x0, y0));
xlabel(ax3, 'x');
style_library_axes(ax3, xLim, yLim);
legend(ax3, hCirc, arrayfun(@(s) sprintf('\\sigma = %.0f px', s), sigmaShow, 'UniformOutput', false), ...
    'Location', 'southoutside');

colormap(parula(max(numel(sigmaList), 2)));
if exist('sgtitle', 'file') == 2
    sgtitle(sprintf('%s | same library for every site', figName));
end
end

function unionMask = build_union_occupancy_mask(ALLCOORDS, RTAB384, nStim, W, H)
unionMask = false(H, W);
for stimNum = 1:nStim
    [~, masks] = render_stim_with_masks2(ALLCOORDS, RTAB384, stimNum, ...
        'ImageSize', [W H], 'DrawDots', false);
    unionMask = unionMask | masks.yellowArm | masks.purple;
end
end

function [xOcc, yOcc] = mask_to_xy(unionMask, W, H, maxPoints)
[rowOcc, colOcc] = find(unionMask);
if nargin >= 4 && isfinite(maxPoints) && numel(rowOcc) > maxPoints
    idx = round(linspace(1, numel(rowOcc), maxPoints));
    rowOcc = rowOcc(idx);
    colOcc = colOcc(idx);
end
xOcc = double(colOcc) - W/2;
yOcc = H/2 - double(rowOcc);
end

function bboxPos = bbox_rect_position(bbox, W, H)
xMin = double(bbox.colMin) - W/2;
xMax = double(bbox.colMax) - W/2;
yMin = H/2 - double(bbox.rowMax);
yMax = H/2 - double(bbox.rowMin);
bboxPos = [xMin, yMin, xMax - xMin, yMax - yMin];
end

function draw_library_bbox(ax, bboxPos)
rectangle(ax, 'Position', bboxPos, 'EdgeColor', [0.85 0.2 0.2], 'LineWidth', 1.5, ...
    'LineStyle', '--');
end

function style_library_axes(ax, xLim, yLim)
axis(ax, 'equal');
xlim(ax, xLim);
ylim(ax, yLim);
grid(ax, 'on');
box(ax, 'on');
end

function [x0, y0] = representative_library_center(Library, bboxPos)
targetX = bboxPos(1) + bboxPos(3) / 2;
targetY = bboxPos(2) + bboxPos(4) / 2;
[~, idx] = min((Library.centerX - targetX).^2 + (Library.centerY - targetY).^2);
x0 = Library.centerX(idx);
y0 = Library.centerY(idx);
end

function h = plot_sigma_circle(ax, x0, y0, sigma, clr, lineWidth)
t = linspace(0, 2*pi, 181);
h = plot(ax, x0 + sigma * cos(t), y0 + sigma * sin(t), ...
    'Color', clr, 'LineWidth', lineWidth);
end

function make_spatial_example_figure(Library, FitStruct, ve, globalSites, unionMask, figName, nExampleSites)
siteSel = select_diverse_spatial_examples(FitStruct, ve, nExampleSites);
if isempty(siteSel)
    return;
end

W = Library.imageSize(1);
H = Library.imageSize(2);
[xOcc, yOcc] = mask_to_xy(unionMask, W, H, 30000);
bboxPos = bbox_rect_position(Library.bbox, W, H);
xSpan = bboxPos(3);
ySpan = bboxPos(4);
xMargin = max(20, 0.08 * xSpan);
yMargin = max(20, 0.08 * ySpan);
xLim = [bboxPos(1) - xMargin, bboxPos(1) + bboxPos(3) + xMargin];
yLim = [bboxPos(2) - yMargin, bboxPos(2) + bboxPos(4) + yMargin];
nPanels = numel(siteSel);
nCols = min(2, nPanels);
nRows = ceil(nPanels / nCols);

figure('Color', 'w');
useTiled = exist('tiledlayout', 'file') == 2;
if useTiled
    tl = tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
end
for i = 1:numel(siteSel)
    s = siteSel(i);
    if useTiled
        ax = nexttile(tl);
    else
        ax = subplot(nRows, nCols, i);
    end
    plot(ax, xOcc, yOcc, '.', 'Color', [0.88 0.88 0.88], 'MarkerSize', 3);
    hold(ax, 'on');
    draw_library_bbox(ax, bboxPos);
    fitClr = fit_gain_color(FitStruct(s).gain);
    draw_gaussian_footprint(ax, FitStruct(s).centerX, FitStruct(s).centerY, FitStruct(s).sigmaPx, fitClr);
    plot(ax, FitStruct(s).centerX, FitStruct(s).centerY, 'k+', 'MarkerSize', 8, 'LineWidth', 1.2);
    xlabel(ax, 'x');
    if i == 1
        ylabel(ax, 'y');
    end
    title(ax, sprintf('site %d | VE %.1f%%\n\\sigma %.0f px | gain %.2f', ...
        globalSites(s), FitStruct(s).r2TrainPct, FitStruct(s).sigmaPx, FitStruct(s).gain), ...
        'FontSize', 10);
    style_library_axes(ax, xLim, yLim);
end
if exist('sgtitle', 'file') == 2
    sgtitle(figName);
end
end

function siteSel = select_diverse_spatial_examples(FitStruct, ve, nExampleSites)
ok = strcmp(string({FitStruct.status}).', "ok") & isfinite(ve);
idx = find(ok);
if isempty(idx)
    siteSel = [];
    return;
end

sigmaAll = [FitStruct.sigmaPx].';
sigmaIdx = sigmaAll(idx);
isGoodSigma = isfinite(sigmaIdx);
idx = idx(isGoodSigma);
sigmaIdx = sigmaIdx(isGoodSigma);
if isempty(idx)
    siteSel = [];
    return;
end

uSigma = unique(sigmaIdx, 'sorted');
repIdx = nan(numel(uSigma), 1);
repVE = nan(numel(uSigma), 1);
for iSig = 1:numel(uSigma)
    cand = idx(sigmaIdx == uSigma(iSig));
    [repVE(iSig), bestLocal] = max(ve(cand)); %#ok<ASGLU>
    repIdx(iSig) = cand(bestLocal);
end

nTake = min(nExampleSites, numel(repIdx));
if numel(repIdx) > nTake
    pick = round(linspace(1, numel(repIdx), nTake));
    pick = unique(pick, 'stable');
    if numel(pick) < nTake
        fill = setdiff(1:numel(repIdx), pick, 'stable');
        pick = [pick(:); fill(1:(nTake - numel(pick))).']; %#ok<AGROW>
    end
    siteSel = repIdx(pick);
else
    siteSel = repIdx;
end

if numel(siteSel) < nExampleSites
    [~, ordRest] = sort(ve(idx), 'descend');
    rest = idx(ordRest);
    rest = rest(~ismember(rest, siteSel));
    need = min(nExampleSites - numel(siteSel), numel(rest));
    siteSel = [siteSel(:); rest(1:need)];
end

[~, ordFinal] = sortrows([[FitStruct(siteSel).sigmaPx].', ve(siteSel)], [1 -2]);
siteSel = siteSel(ordFinal);
end

function [pairsA, pairsB] = build_complementary_pairs(nStim)
pairsA = [];
pairsB = [];
for stimA = 1:nStim
    pos = mod(stimA - 1, 8) + 1;
    if pos <= 4 && stimA + 4 <= nStim
        pairsA(end+1,1) = stimA; %#ok<AGROW>
        pairsB(end+1,1) = stimA + 4; %#ok<AGROW>
    end
end
end

function tf = is_color_dominant(mainOcc, otherOcc, ratioThr)
occFloor = 1e-9;
tf = isfinite(mainOcc) && isfinite(otherOcc) && (mainOcc > 0) && ...
    (mainOcc >= ratioThr * max(otherOcc, occFloor));
end

function clr = fit_gain_color(gainVal)
if ~isfinite(gainVal)
    clr = [0.35 0.35 0.35];
elseif gainVal >= 0
    clr = [0.88 0.42 0.16];
else
    clr = [0.19 0.45 0.78];
end
end

function draw_gaussian_footprint(ax, x0, y0, sigma, clr)
[x1, y1] = circle_points(x0, y0, sigma);
[x2, y2] = circle_points(x0, y0, 2 * sigma);
patch(ax, x1, y1, clr, 'FaceAlpha', 0.18, 'EdgeColor', 'none');
plot(ax, x1, y1, '-', 'Color', clr, 'LineWidth', 1.8);
plot(ax, x2, y2, '--', 'Color', clr, 'LineWidth', 1.1);
end

function [xCirc, yCirc] = circle_points(x0, y0, radius)
t = linspace(0, 2*pi, 181);
xCirc = x0 + radius * cos(t);
yCirc = y0 + radius * sin(t);
end

function T = build_spatial_table(FitStruct, RFrange)
n = numel(FitStruct);
statusText = string({FitStruct.status}).';
T = table((1:n).', RFrange(:), [FitStruct.nObs].', [FitStruct.centerX].', ...
    [FitStruct.centerY].', [FitStruct.sigmaPx].', [FitStruct.gain].', ...
    [FitStruct.effectRangeObs].', [FitStruct.r2TrainPct].', [FitStruct.pValueApprox].', ...
    statusText, 'VariableNames', ...
    {'localSite','globalSiteInR','nObs','centerX','centerY','sigmaPx','gain', ...
    'effectRangeObs','r2TrainPct','pValueApprox','status'});
T = sortrows(T(isfinite(T.r2TrainPct), :), 'r2TrainPct', 'descend');
end

function T = build_color_table(FitColor, FitSpatial, RFrange)
n = numel(FitColor);
statusText = string({FitColor.status}).';
wColorDiff = [FitColor.wYellow].' - [FitColor.wPurple].';
wOccMean = 0.5 * ([FitColor.wYellow].' + [FitColor.wPurple].');
T = table((1:n).', RFrange(:), [FitColor.nObs].', [FitColor.centerX].', [FitColor.centerY].', ...
    [FitColor.sigmaPx].', [FitColor.wYellow].', [FitColor.wPurple].', wColorDiff, wOccMean, ...
    [FitColor.deltaR2OverSpatial].', [FitColor.pValueAddColor].', [FitSpatial.r2TrainPct].', ...
    statusText, 'VariableNames', ...
    {'localSite','globalSiteInR','nObs','centerX','centerY','sigmaPx','wYellow','wPurple', ...
    'wYellowMinusPurple','wOccupancyMean', ...
    'deltaR2OverSpatial','pValueAddColor','spatialR2TrainPct','status'});
T = sortrows(T(isfinite(T.deltaR2OverSpatial), :), 'deltaR2OverSpatial', 'descend');
end

function T = build_paired_color_table(PairColor, RFrange)
n = numel(PairColor);
statusText = string({PairColor.status}).';
T = table((1:n).', RFrange(:), [PairColor.nPairs].', [PairColor.centerX].', [PairColor.centerY].', ...
    [PairColor.sigmaPx].', [PairColor.meanYellowDominant].', [PairColor.meanPurpleDominant].', ...
    [PairColor.weightedYellowDominant].', [PairColor.weightedPurpleDominant].', ...
    [PairColor.muYellowPooled].', [PairColor.muPurplePooled].', ...
    [PairColor.varYellowPooled].', [PairColor.varPurplePooled].', ...
    [PairColor.NYPooled].', [PairColor.NPPooled].', [PairColor.zStatPooled].', [PairColor.pValuePooled].', ...
    [PairColor.meanPairDiff].', [PairColor.weightedPairDiff].', [PairColor.colorIndex].', [PairColor.colorIndexDenom].', [PairColor.meanPairVar].', ...
    [PairColor.sumPairWeight].', [PairColor.effPairN].', [PairColor.tStat].', [PairColor.pValue].', ...
    statusText, 'VariableNames', ...
    {'localSite','globalSiteInR','nPairs','centerX','centerY','sigmaPx', ...
    'meanYellowDominant','meanPurpleDominant','weightedYellowDominant','weightedPurpleDominant', ...
    'muYellowPooled','muPurplePooled','varYellowPooled','varPurplePooled','NYPooled','NPPooled','zStatPooled','pValuePooled', ...
    'meanPairDiff','weightedPairDiff','colorIndex','colorIndexDenom','meanPairVar','sumPairWeight','effPairN','tStat','pValue','status'});
T = sortrows(T(isfinite(T.pValuePooled), :), 'pValuePooled', 'ascend');
end

function p = fcdf_local(x, df1, df2)
z = (df1 .* x) ./ (df1 .* x + df2);
z = min(max(z, 0), 1);
p = betainc(z, df1/2, df2/2);
end
