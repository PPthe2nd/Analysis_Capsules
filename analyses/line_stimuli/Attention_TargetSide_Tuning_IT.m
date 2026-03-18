function OUT = Attention_TargetSide_Tuning_IT()
% ATTENTION_TARGETSIDE_TUNING_IT
% Fit IT tuning to the screen direction of the attended capsule relative to
% the unattended capsule.
%
% For each 4-stimulus quartet, the response is the color-balanced
% attention-direction contrast between the two opposite attended
% directions:
%   deltaQuartet = muDirA - muDirB
% where muDirA is the trial-weighted mean across the two color-swapped
% stimuli with one attended direction, and muDirB is the corresponding mean
% across the two color-swapped stimuli with the opposite attended
% direction. The normalized effect is:
%   (muDirA - muDirB) / sdSpont
%
% The corresponding screen angle is the target-side direction for DirA,
% defined as the vector distractor -> target:
%   0 deg = target to the right of distractor
%   90 deg = target above distractor
%   180 deg = target to the left
%   270 deg = target below
%
% A weighted first-harmonic model is fit separately in the early and late
% 3-bin windows:
%   direction: y = b + c*cos(phi) + s*sin(phi)
%
% Because the contrast is formed within each quartet, fixed stimulus
% rotation/orientation effects are canceled out. The fit therefore asks
% whether the attention effect itself has a preferred screen direction.

%% Settings
P = struct();
P.Monkey = 1;           % 1 = Nilson, 2 = Figaro
P.minQuartets = 20;     % minimum number of usable direction-contrast quartets required per site
P.makeSummaryFigures = true;
P.makeExampleFigures = true;
P.makeAngleReferenceFigure = true;
P.nExampleSites = 4;
P.saveResult = true;
P.forceRefit = false;
P.progressEvery = 10;
P.vePlotFloorPct = -100;
P.sigAlpha = 0.05;
P.varFloorFrac = 1e-3;

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
    error('Attention_TargetSide_Tuning_IT:InvalidMonkey', ...
        'P.Monkey must be 1 (Nilson) or 2 (Figaro).');
end

tallPath = fullfile(cfg.matDir, tallFile);
resp3binPath = fullfile(cfg.matDir, resp3binFile);
outPath = fullfile(cfg.matDir, sprintf('Attention_TargetSide_Tuning_IT_directiondelta_%s.mat', char(monkeySuffix)));

assert(exist(tallPath, 'file') == 2, ...
    'Missing %s. Run Line_Stimuli_IT.m first.', tallPath);
assert(exist(resp3binPath, 'file') == 2, ...
    'Missing %s. Create the 3-bin response summary first.', resp3binPath);

useCached = exist(outPath, 'file') == 2 && ~P.forceRefit;
needsSave = ~useCached;

if useCached
    fprintf('Loading cached IT target-side tuning from %s\n', outPath);
    S = load(outPath);
    assert(isfield(S, 'OUT') && isstruct(S.OUT), ...
        '%s must contain struct OUT.', outPath);
    OUT = S.OUT;
    assert(isfield(OUT, 'QuartetTable') && isfield(OUT, 'FitEarly') && ...
        isfield(OUT, 'FitLate') && isfield(OUT, 'deltaQuartetEarly') && ...
        isfield(OUT, 'deltaQuartetLate') && isfield(OUT, 'RFrange'), ...
        'Cached OUT is missing required fields.');

    QuartetTable = OUT.QuartetTable;
    thetaDeg = QuartetTable.targetDirDeg;
    FitEarly = OUT.FitEarly;
    FitLate = OUT.FitLate;
    deltaQuartetEarly = OUT.deltaQuartetEarly;
    deltaQuartetLate = OUT.deltaQuartetLate;
    varQuartetEarly = OUT.varQuartetEarly;
    varQuartetLate = OUT.varQuartetLate;
    RFrange = OUT.RFrange(:);
    nIT = numel(RFrange);
else
    %% Load geometry
    Sgeo = load(tallPath);
    assert(isfield(Sgeo, 'Tall_IT') && isstruct(Sgeo.Tall_IT), ...
        '%s must contain struct Tall_IT.', tallFile);
    assert(isfield(Sgeo, 'ALLCOORDS') && isfield(Sgeo, 'RTAB384') && isfield(Sgeo, 'RFrange'), ...
        '%s must contain ALLCOORDS, RTAB384 and RFrange.', tallFile);

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

    [quartetPairsA, quartetPairsB, thetaDeg, stimRef, pairMidXY] = ...
        build_quartets_and_target_direction(ALLCOORDS, nStim);
    nQuartets = size(quartetPairsA, 1);

    QuartetTable = table((1:nQuartets).', stimRef(:), pairMidXY(:,1), pairMidXY(:,2), thetaDeg(:), ...
        'VariableNames', {'quartetIdx','stimRef','pairMidX','pairMidY','targetDirDeg'});

    %% Load 3-bin responses and localize to IT rows
    Sresp = load(resp3binPath);
    assert(isfield(Sresp, 'R') && isstruct(Sresp.R), ...
        '%s must contain struct R.', resp3binFile);

    R3_full = Sresp.R;
    R3 = R3_full;
    R3.meanAct = R3_full.meanAct(RFrange, :, :);
    R3.meanSqAct = R3_full.meanSqAct(RFrange, :, :);
    if ismatrix(R3_full.nTrials) && size(R3_full.nTrials,1) >= max(RFrange)
        R3.nTrials = R3_full.nTrials(RFrange, :);
    else
        R3.nTrials = R3_full.nTrials;
    end

    assert(size(R3.meanAct,1) == nIT, 'Localized IT response rows do not match Tall_IT.');
    assert(size(R3.meanAct,2) == nStim, 'Localized IT response stimuli do not match Tall_IT.');

    SNR = compute_snr_per_color_sites(R3, Tall_IT, siteRows, 'Verbose', false);
    muSpont = SNR.muSpont(siteRows);
    sdSpont = SNR.sdSpont(siteRows);

    [deltaQuartetEarly, deltaQuartetLate, varQuartetEarly, varQuartetLate] = ...
        compute_direction_delta_quartets_with_var(R3, quartetPairsA, quartetPairsB, muSpont, sdSpont);

    %% Fit direction model
    fields = {'status','nObs','baseline','coefCos','coefSin','amplitude','prefDeg', ...
        'rmse','r2TrainPct','effectRangeObs','effectAbsPeakObs','fStatApprox', ...
        'pValueApprox','dfModel','dfError','rssNullWeighted','rssFullWeighted','modelKind'};
    FitEarly = init_fit_struct(nIT, fields, "direction");
    FitLate = init_fit_struct(nIT, fields, "direction");

    fprintf('Fitting IT target-side attention-direction tuning for monkey %s (%d sites, %d quartets)\n', ...
        char(monkeySuffix), nIT, nQuartets);
    tFit = tic;

    for iSite = 1:nIT
        yEarly = deltaQuartetEarly(iSite, :).';
        yLate = deltaQuartetLate(iSite, :).';
        vEarly = varQuartetEarly(iSite, :).';
        vLate = varQuartetLate(iSite, :).';

        FitEarly(iSite) = fit_weighted_harmonic(thetaDeg, yEarly, vEarly, 1, "direction", P);
        FitLate(iSite) = fit_weighted_harmonic(thetaDeg, yLate, vLate, 1, "direction", P);

        if P.progressEvery > 0 && (iSite == 1 || mod(iSite, P.progressEvery) == 0 || iSite == nIT)
            elapsedSec = toc(tFit);
            rate = iSite / max(elapsedSec, eps);
            etaSec = (nIT - iSite) / max(rate, eps);
            fprintf('  Site %d / %d | elapsed %.1fs | ETA %.1fs\n', ...
                iSite, nIT, elapsedSec, etaSec);
        end
    end
end

if (P.makeSummaryFigures || P.makeExampleFigures || P.makeAngleReferenceFigure) && ...
        (~exist('ALLCOORDS', 'var') || ~exist('RTAB384', 'var'))
    Sref = load(tallPath, 'ALLCOORDS', 'RTAB384');
    assert(isfield(Sref, 'ALLCOORDS') && isfield(Sref, 'RTAB384'), ...
        '%s must contain ALLCOORDS and RTAB384.', tallPath);
    ALLCOORDS = Sref.ALLCOORDS;
    RTAB384 = Sref.RTAB384;
end

veEarly = [FitEarly.r2TrainPct].';
veLate = [FitLate.r2TrainPct].';
effEarly = [FitEarly.effectRangeObs].';
effLate = [FitLate.effectRangeObs].';
pEarly = [FitEarly.pValueApprox].';
pLate = [FitLate.pValueApprox].';
isDirEarly = isfinite(pEarly) & (pEarly < P.sigAlpha);
isDirLate = isfinite(pLate) & (pLate < P.sigAlpha);

fprintf('Usable IT sites for target-direction fitting (early): %d / %d\n', ...
    nnz(isfinite(veEarly)), nIT);
fprintf('Usable IT sites for target-direction fitting (late):  %d / %d\n', ...
    nnz(isfinite(veLate)), nIT);
fprintf('Direction-tuned sites at p < %.3f | early=%d late=%d\n', ...
    P.sigAlpha, nnz(isDirEarly), nnz(isDirLate));

%% Summary figures
if P.makeSummaryFigures
    figSum = figure('Color', 'w');
    useTiled = exist('tiledlayout', 'file') == 2;
    if useTiled
        tiledlayout(1,3, 'TileSpacing', 'compact', 'Padding', 'compact');
    end

    if useTiled, nexttile; else, subplot(1,3,1); end
    plot_ve_histogram(veEarly, veLate, P.vePlotFloorPct, ...
        sprintf('Target-direction VE (%s)', char(monkeySuffix)));

    if useTiled
        axTmp = nexttile;
    else
        axTmp = subplot(1,3,2);
    end
    axPos = axTmp.Position;
    delete(axTmp);
    pax = polaraxes('Parent', figSum, 'Position', axPos);
    plot_pref_direction_polar_histogram(pax, mod([FitLate(isDirLate).prefDeg], 360), ...
        sprintf('Late significant direction sites (%s)', char(monkeySuffix)));

    if useTiled, nexttile; else, subplot(1,3,3); end
    plot_effect_histogram(effEarly, effLate, ...
        sprintf('Target-direction effect size (%s)', char(monkeySuffix)));
end

if P.makeAngleReferenceFigure
    make_target_side_reference_figure(ALLCOORDS, RTAB384, QuartetTable, ...
        sprintf('IT target-side angle reference (%s)', char(monkeySuffix)));
end

%% Example figures
if P.makeExampleFigures
    siteSelEarly = select_top_examples(veEarly, P.nExampleSites);
    siteSelLate = select_top_examples(veLate, P.nExampleSites);

    make_example_angle_fit_figure(thetaDeg, deltaQuartetEarly, FitEarly, ...
        sprintf('IT target-direction fits early (%s)', char(monkeySuffix)), siteSelEarly, 1);
    make_example_angle_fit_figure(thetaDeg, deltaQuartetLate, FitLate, ...
        sprintf('IT target-direction fits late (%s)', char(monkeySuffix)), siteSelLate, 1);
end

%% Tables
Tearly = build_fit_table(FitEarly, RFrange, 'r2TrainPct');
Tlate = build_fit_table(FitLate, RFrange, 'r2TrainPct');

disp('Top IT target-side direction sites by variance explained (late):');
disp(Tlate(1:min(10,height(Tlate)), :));

%% Pack output
OUT.P = P;
OUT.monkeySuffix = monkeySuffix;
OUT.RFrange = RFrange;
OUT.QuartetTable = QuartetTable;
OUT.deltaQuartetEarly = deltaQuartetEarly;
OUT.deltaQuartetLate = deltaQuartetLate;
OUT.varQuartetEarly = varQuartetEarly;
OUT.varQuartetLate = varQuartetLate;
OUT.FitEarly = FitEarly;
OUT.FitLate = FitLate;
OUT.isDirectionTunedEarly = isDirEarly;
OUT.isDirectionTunedLate = isDirLate;
OUT.TableEarly = Tearly;
OUT.TableLate = Tlate;

if needsSave && P.saveResult
    save(outPath, 'OUT', '-v7.3');
    fprintf('Saved IT target-side tuning to %s\n', outPath);
end
end

function [quartetPairsA, quartetPairsB, thetaDeg, stimRef, pairMidXY] = build_quartets_and_target_direction(ALLCOORDS, nStim)
nQuartets = nStim / 4;
quartetPairsA = zeros(nQuartets, 2);
quartetPairsB = zeros(nQuartets, 2);
thetaDeg = nan(nQuartets, 1);
stimRef = nan(nQuartets, 1);
pairMidXY = nan(nQuartets, 2);
q = 0;

for base = 0:8:(nStim - 8)
    q = q + 1;
    quartetPairsA(q,:) = base + [1 5];
    quartetPairsB(q,:) = base + [2 6];
    [thetaDeg(q), stimRef(q), pairMidXY(q,:)] = target_direction_from_quartet(ALLCOORDS, quartetPairsA(q,:), quartetPairsB(q,:));

    q = q + 1;
    quartetPairsA(q,:) = base + [3 7];
    quartetPairsB(q,:) = base + [4 8];
    [thetaDeg(q), stimRef(q), pairMidXY(q,:)] = target_direction_from_quartet(ALLCOORDS, quartetPairsA(q,:), quartetPairsB(q,:));
end

quartetPairsA = quartetPairsA(1:q, :);
quartetPairsB = quartetPairsB(1:q, :);
thetaDeg = thetaDeg(1:q);
stimRef = stimRef(1:q);
pairMidXY = pairMidXY(1:q, :);
end

function [thetaDeg, stimRef, pairMid] = target_direction_from_quartet(ALLCOORDS, pairA, pairB)
stimRef = pairA(1);
[thetaA, pairMidA] = target_side_from_pair_local(ALLCOORDS, pairA);
[thetaB, pairMidB] = target_side_from_pair_local(ALLCOORDS, pairB);

assert(norm(pairMidB - pairMidA) < 1e-6, ...
    'Pair midpoint differs within direction quartet [%s] vs [%s].', num2str(pairA), num2str(pairB));
assert(abs(abs(angdiff_deg_local(thetaB, thetaA)) - 180) < 1e-6, ...
    'Target directions are not opposite within direction quartet [%s] vs [%s].', num2str(pairA), num2str(pairB));

thetaDeg = mod(thetaA, 360);
pairMid = pairMidA;
end

function [thetaDeg, pairMid] = target_side_from_pair_local(ALLCOORDS, stimPair)
[theta0, pairMid0] = target_side_angle_from_stim(ALLCOORDS, stimPair(1));
for k = 2:numel(stimPair)
    [thetaK, pairMidK] = target_side_angle_from_stim(ALLCOORDS, stimPair(k));
    assert(norm(pairMidK - pairMid0) < 1e-6, ...
        'Pair midpoint differs within pair [%s].', num2str(stimPair));
    assert(abs(angdiff_deg_local(thetaK, theta0)) < 1e-6, ...
        'Target-side angle differs within pair [%s].', num2str(stimPair));
end
thetaDeg = mod(theta0, 360);
pairMid = pairMid0;
end

function [thetaDeg, pairMid] = target_side_angle_from_stim(ALLCOORDS, stimNum)
f = sprintf('stim_%d', stimNum);
tFig = double(ALLCOORDS.(f).t_fig(:))';
tBack = double(ALLCOORDS.(f).t_back(:))';

v = tFig - tBack;
if norm(v) <= 0
    error('Attention_TargetSide_Tuning_IT:DegenerateTargetSide', ...
        'Degenerate target-side vector for stimulus %d.', stimNum);
end

thetaDeg = mod(atan2d(v(2), v(1)), 360);
pairMid = 0.5 * (tFig + tBack);
end

function [deltaEarly, deltaLate, varEarly, varLate] = compute_direction_delta_quartets_with_var( ...
        R3, quartetPairsA, quartetPairsB, muSpont, sdSpont)
[nSites, ~, ~] = size(R3.meanAct);
nQuartets = size(quartetPairsA, 1);

if isvector(R3.nTrials)
    nTrialsAll = double(R3.nTrials(:));
    perSiteTrials = false;
else
    perSiteTrials = true;
end

deltaEarly = nan(nSites, nQuartets);
deltaLate = nan(nSites, nQuartets);
varEarly = nan(nSites, nQuartets);
varLate = nan(nSites, nQuartets);

for iSite = 1:nSites
    if perSiteTrials
        nTrSite = double(R3.nTrials(iSite, :)).';
    else
        nTrSite = nTrialsAll;
    end

    rEarly = squeeze(R3.meanAct(iSite,:,2)).';
    rLate = squeeze(R3.meanAct(iSite,:,3)).';
    rSqEarly = squeeze(R3.meanSqAct(iSite,:,2)).';
    rSqLate = squeeze(R3.meanSqAct(iSite,:,3)).';

    for qIdx = 1:nQuartets
        [muAEarly, varAEarly] = compute_pairmean_for_window_local( ...
            quartetPairsA(qIdx,:), nTrSite, rEarly, rSqEarly, muSpont(iSite), sdSpont(iSite));
        [muBEarly, varBEarly] = compute_pairmean_for_window_local( ...
            quartetPairsB(qIdx,:), nTrSite, rEarly, rSqEarly, muSpont(iSite), sdSpont(iSite));
        if isfinite(muAEarly) && isfinite(muBEarly)
            deltaEarly(iSite, qIdx) = muAEarly - muBEarly;
            varEarly(iSite, qIdx) = varAEarly + varBEarly;
        end

        [muALate, varALate] = compute_pairmean_for_window_local( ...
            quartetPairsA(qIdx,:), nTrSite, rLate, rSqLate, muSpont(iSite), sdSpont(iSite));
        [muBLate, varBLate] = compute_pairmean_for_window_local( ...
            quartetPairsB(qIdx,:), nTrSite, rLate, rSqLate, muSpont(iSite), sdSpont(iSite));
        if isfinite(muALate) && isfinite(muBLate)
            deltaLate(iSite, qIdx) = muALate - muBLate;
            varLate(iSite, qIdx) = varALate + varBLate;
        end
    end
end

badNoise = ~isfinite(sdSpont) | (sdSpont <= 0);
deltaEarly(badNoise, :) = NaN;
deltaLate(badNoise, :) = NaN;
varEarly(badNoise, :) = NaN;
varLate(badNoise, :) = NaN;
end

function [pairMeanNorm, varNorm] = compute_pairmean_for_window_local(stimSet, nTrSite, rWin, rSqWin, muSpont, sdSpont)
pairMeanNorm = NaN;
varNorm = NaN;

nUse = nTrSite(stimSet);
muUse = rWin(stimSet);
msqUse = rSqWin(stimSet);
good = isfinite(muUse) & isfinite(msqUse) & isfinite(nUse) & (nUse > 1);

if ~(any(good) && isfinite(sdSpont) && (sdSpont > 0))
    return;
end

nUse = nUse(good);
muUse = muUse(good);
msqUse = msqUse(good);
varStim = max(0, msqUse - muUse.^2);
varStim = varStim .* (nUse ./ max(nUse - 1, 1));
muPair = sum(nUse .* muUse) / sum(nUse);
varPair = sum(nUse .* varStim) / (sum(nUse)^2);

pairMeanNorm = (muPair - muSpont) / sdSpont;
varNorm = varPair / (sdSpont^2);
end

function F = init_fit_struct(nSites, fields, modelKind)
tmp = struct();
for i = 1:numel(fields)
    tmp.(fields{i}) = NaN;
end
tmp.status = "not_fit";
tmp.modelKind = modelKind;
F = repmat(tmp, nSites, 1);
end

function F = fit_weighted_harmonic(thetaDeg, y, varY, harmonic, modelKind, P)
fields = {'status','nObs','baseline','coefCos','coefSin','amplitude','prefDeg', ...
    'rmse','r2TrainPct','effectRangeObs','effectAbsPeakObs','fStatApprox', ...
    'pValueApprox','dfModel','dfError','rssNullWeighted','rssFullWeighted','modelKind'};
F = struct();
for i = 1:numel(fields)
    F.(fields{i}) = NaN;
end
F.status = "not_fit";
F.modelKind = modelKind;

valid = isfinite(thetaDeg) & isfinite(y) & isfinite(varY) & (varY > 0);
theta = thetaDeg(valid);
yv = y(valid);
vv = varY(valid);
F.nObs = numel(yv);

if numel(yv) < P.minQuartets
    F.status = "too_few_quartets";
    return;
end

varFloor = max(median(vv(isfinite(vv) & vv > 0)) * P.varFloorFrac, 1e-6);
vv = max(vv, varFloor);
sqrtW = sqrt(1 ./ vv);

X = [ones(numel(theta),1), cosd(harmonic * theta), sind(harmonic * theta)];
Xw = bsxfun(@times, X, sqrtW);
yw = sqrtW .* yv;

beta = Xw \ yw;
yHat = X * beta;
wMean = sum((1 ./ vv) .* yv) / sum(1 ./ vv);
rssNull = sum((1 ./ vv) .* (yv - wMean).^2);
rssFull = sum((1 ./ vv) .* (yv - yHat).^2);
dfModel = 2;
dfError = numel(yv) - size(X,2);

F.status = "ok";
F.baseline = beta(1);
F.coefCos = beta(2);
F.coefSin = beta(3);
F.amplitude = hypot(beta(2), beta(3));
if harmonic == 1
    F.prefDeg = mod(atan2d(beta(3), beta(2)), 360);
else
    F.prefDeg = mod(0.5 * atan2d(beta(3), beta(2)), 180);
end
F.rmse = sqrt(mean((yv - yHat).^2));
if rssNull > 0
    F.r2TrainPct = 100 * (1 - rssFull / rssNull);
end
F.effectRangeObs = max(yHat) - min(yHat);
F.effectAbsPeakObs = max(abs(yHat));
F.rssNullWeighted = rssNull;
F.rssFullWeighted = rssFull;
F.dfModel = dfModel;
F.dfError = dfError;

if dfError > 0 && rssNull > rssFull
    F.fStatApprox = ((rssNull - rssFull) / dfModel) / (rssFull / dfError);
    F.pValueApprox = 1 - fcdf_local(F.fStatApprox, dfModel, dfError);
else
    F.fStatApprox = 0;
    F.pValueApprox = 1;
end
end

function siteSel = select_top_examples(ve, nSelect)
idx = find(isfinite(ve));
if isempty(idx)
    siteSel = [];
    return;
end
[~, ord] = sort(ve(idx), 'descend');
siteSel = idx(ord(1:min(nSelect, numel(ord))));
end

function plot_ve_histogram(veEarly, veLate, floorPct, ttl)
vePlotEarly = veEarly;
vePlotLate = veLate;
vePlotEarly(isfinite(vePlotEarly) & (vePlotEarly < floorPct)) = floorPct;
vePlotLate(isfinite(vePlotLate) & (vePlotLate < floorPct)) = floorPct;

allVals = [vePlotEarly(isfinite(vePlotEarly)); vePlotLate(isfinite(vePlotLate))];
if isempty(allVals)
    edges = linspace(floorPct, floorPct + 1, 31);
else
    lo = min(allVals);
    hi = max(allVals);
    if ~(isfinite(lo) && isfinite(hi)) || hi <= lo
        hi = lo + 1;
    end
    edges = linspace(lo, hi, 31);
end

histogram(vePlotEarly(isfinite(vePlotEarly)), edges, 'FaceColor', [0.25 0.45 0.85], 'EdgeColor', 'none');
hold on;
histogram(vePlotLate(isfinite(vePlotLate)), edges, 'FaceColor', [0.85 0.35 0.25], 'EdgeColor', 'none');
xlabel('Variance explained (%)');
ylabel('N sites');
title(ttl);
legend('Early','Late');
grid on;
end

function plot_effect_histogram(effEarly, effLate, ttl)
allVals = [effEarly(isfinite(effEarly)); effLate(isfinite(effLate))];
if isempty(allVals)
    edges = linspace(0, 1, 31);
else
    lo = min(allVals);
    hi = max(allVals);
    if ~(isfinite(lo) && isfinite(hi)) || hi <= lo
        hi = lo + 1;
    end
    edges = linspace(lo, hi, 31);
end

histogram(effEarly(isfinite(effEarly)), edges, 'FaceColor', [0.25 0.45 0.85], 'EdgeColor', 'none');
hold on;
histogram(effLate(isfinite(effLate)), edges, 'FaceColor', [0.85 0.35 0.25], 'EdgeColor', 'none');
xlabel('Fitted modulation range (spont. SD units)');
ylabel('N sites');
title(ttl);
legend('Early','Late');
grid on;
end

function make_target_side_reference_figure(ALLCOORDS, RTAB384, QuartetTable, figName)
dirTargets = [0 90 180 270];
theta = double(QuartetTable.targetDirDeg(:));
stimRef = double(QuartetTable.stimRef(:));

dirStim = nan(size(dirTargets));
dirTheta = nan(size(dirTargets));
for i = 1:numel(dirTargets)
    d = abs(angdiff_deg_local(theta, dirTargets(i)));
    [~, idx] = min(d);
    dirStim(i) = stimRef(idx);
    dirTheta(i) = theta(idx);
end

fig = figure('Name', figName, 'NumberTitle', 'off', 'Color', 'w');
useTiled = exist('tiledlayout', 'file') == 2;
if useTiled
    tl = tiledlayout(fig, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
end

for i = 1:4
    if useTiled, ax = nexttile(tl); else, ax = subplot(1,4,i); end
    img = render_stim_from_ALLCOORDS(ALLCOORDS, RTAB384, dirStim(i));
    imshow(img, 'Parent', ax, 'InitialMagnification', 'fit');
    axis(ax, 'image');
    title(ax, sprintf('%d = %s\nstim %d (%.1f°)', dirTargets(i), dir_label(dirTargets(i)), dirStim(i), dirTheta(i)), ...
        'FontSize', 10);
end

if exist('sgtitle', 'file') == 2
    sgtitle(figName);
end
end

function lbl = dir_label(thetaDeg)
switch mod(thetaDeg, 360)
    case 0
        lbl = 'target R';
    case 90
        lbl = 'target U';
    case 180
        lbl = 'target L';
    case 270
        lbl = 'target D';
    otherwise
        lbl = '?';
end
end

function plot_pref_direction_polar_histogram(ax, prefDeg, ttl)
prefDeg = prefDeg(:);
prefDeg = prefDeg(isfinite(prefDeg));

edgesDeg = linspace(0, 360, 19);
if isempty(prefDeg)
    polarhistogram(ax, deg2rad(0), deg2rad(edgesDeg), ...
        'FaceColor', [0.20 0.55 0.85], 'EdgeColor', 'none', ...
        'FaceAlpha', 0.0, 'EdgeAlpha', 0.0);
else
    polarhistogram(ax, deg2rad(prefDeg), deg2rad(edgesDeg), ...
        'FaceColor', [0.20 0.55 0.85], 'EdgeColor', 'w', 'LineWidth', 0.8);
end

ax.ThetaZeroLocation = 'right';
ax.ThetaDir = 'counterclockwise';
ax.ThetaTick = [0 90 180 270];
ax.ThetaTickLabel = {'0=R','90=U','180=L','270=D'};
ax.RAxisLocation = 135;
ax.GridAlpha = 0.20;
ax.MinorGridAlpha = 0.10;
ax.ThetaColor = [0.25 0.25 0.25];
ax.RColor = [0.25 0.25 0.25];
title(ax, ttl);
end

function make_example_angle_fit_figure(thetaDeg, Ymat, FitStruct, figName, siteSel, harmonic)
if isempty(siteSel)
    return;
end

nShow = numel(siteSel);
fig = figure('Name', figName, 'NumberTitle', 'off', 'Color', 'w');
useTiled = exist('tiledlayout', 'file') == 2;
if useTiled
    tl = tiledlayout(fig, 1, nShow, 'TileSpacing', 'compact', 'Padding', 'compact');
end
if harmonic == 1
    thetaGrid = linspace(0, 360, 361).';
else
    thetaGrid = linspace(0, 180, 181).';
end
for i = 1:nShow
    s = siteSel(i);
    if useTiled
        ax = nexttile(tl);
    else
        ax = subplot(1, nShow, i);
    end

    y = Ymat(s, :).';
    valid = isfinite(y) & isfinite(thetaDeg);
    th = mod(thetaDeg(valid), 360);
    yv = y(valid);
    [thSort, ord] = sort(th);
    ySort = yv(ord);

    scatter(ax, thSort, ySort, 45, 'filled');
    hold(ax, 'on');

    fitY = harmonic_curve(thetaGrid, FitStruct(s), harmonic);
    plot(ax, thetaGrid, fitY, 'k-', 'LineWidth', 1.5);
    xlim(ax, [0 360]);
    xticks(ax, [0 90 180 270 360]);
    xticklabels(ax, {'0=R','90=U','180=L','270=D','360=R'});
    xlabel(ax, 'Attention direction (deg)');
    ylabel(ax, 'Attention-direction effect (spont. SD)');
    title(ax, sprintf('site %d\nVE %.1f%% | eff %.2f | p %.3g', s, ...
        FitStruct(s).r2TrainPct, FitStruct(s).effectRangeObs, FitStruct(s).pValueApprox), ...
        'FontSize', 10);
    grid(ax, 'on');
end
if exist('sgtitle', 'file') == 2
    sgtitle(figName);
end
end

function yHat = harmonic_curve(thetaDeg, F, harmonic)
yHat = F.baseline + F.coefCos * cosd(harmonic * thetaDeg) + F.coefSin * sind(harmonic * thetaDeg);
end

function T = build_fit_table(FitStruct, RFrange, sortField)
n = numel(FitStruct);
statusText = string({FitStruct.status}).';
modelKind = string({FitStruct.modelKind}).';
T = table((1:n).', RFrange(:), modelKind, [FitStruct.nObs].', [FitStruct.amplitude].', ...
    [FitStruct.prefDeg].', [FitStruct.effectRangeObs].', [FitStruct.r2TrainPct].', ...
    [FitStruct.pValueApprox].', statusText, ...
    'VariableNames', {'localSite','globalSiteInR','modelKind','nObs','amplitude', ...
    'prefDeg','effectRangeObs','r2TrainPct','pValueApprox','status'});

ok = isfinite(T.(sortField));
T = sortrows(T(ok, :), sortField, 'descend');
end

function p = fcdf_local(x, df1, df2)
z = (df1 .* x) ./ (df1 .* x + df2);
z = min(max(z, 0), 1);
p = betainc(z, df1/2, df2/2);
end

function d = angdiff_deg_local(a, b)
d = mod((a - b) + 180, 360) - 180;
end
