function OUT = AttendedCentroid_Tuning_IT()
% ATTENDEDCENTROID_TUNING_IT
% Fit signed IT tuning to the centroid of the top/front (attended) capsule.
%
% The predictor is stimulus-level, not quartet-level, because the attended
% capsule can change within a quartet. For each site and stimulus, the
% signed response is normalized as:
%   (muStim - muSpont) / sdSpont
%
% A signed 2D elliptical Gaussian is fit separately in the early and late
% windows:
%   y(x,y) = baseline + amplitude * G_2D(x,y)
%
% Output metrics include:
% - weighted in-sample variance explained (%)
% - approximate weighted F-test p-value versus a constant model
% - fitted effect size over the observed centroid positions

%% Settings
P = struct();
P.Monkey = 1;            % 1 = Nilson, 2 = Figaro
P.minStimuli = 80;       % minimum number of finite stimulus responses required per site
P.makeSummaryFigures = true;
P.makeExampleFigure = true;
P.nExampleSites = 4;     % per window
P.saveResult = true;
P.forceRefit = false;
P.progressEvery = 10;
P.vePlotFloorPct = -100;
P.sigAlpha = 0.05;       % per-site significance threshold on the approximate p-value
P.varFloorFrac = 1e-3;   % floor on stimulus-response variance relative to the median

cfg = config();

assert(exist('lsqcurvefit', 'file') == 2, ...
    'AttendedCentroid_Tuning_IT requires lsqcurvefit (Optimization Toolbox).');

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
    error('AttendedCentroid_Tuning_IT:InvalidMonkey', ...
        'P.Monkey must be 1 (Nilson) or 2 (Figaro).');
end

tallPath = fullfile(cfg.matDir, tallFile);
resp3binPath = fullfile(cfg.matDir, resp3binFile);
outPath = fullfile(cfg.matDir, sprintf('AttendedCentroid_Tuning_IT_weighted_%s.mat', char(monkeySuffix)));

assert(exist(tallPath, 'file') == 2, ...
    'Missing %s. Run Line_Stimuli_IT.m first.', tallPath);
assert(exist(resp3binPath, 'file') == 2, ...
    'Missing %s. Create the 3-bin response summary first.', resp3binPath);

useCached = exist(outPath, 'file') == 2 && ~P.forceRefit;
needsSave = ~useCached;

if useCached
    fprintf('Loading cached IT attended-centroid tuning from %s\n', outPath);
    S = load(outPath);
    assert(isfield(S, 'OUT') && isstruct(S.OUT), ...
        '%s must contain struct OUT.', outPath);
    OUT = S.OUT;
    assert(isfield(OUT, 'StimulusTable') && isfield(OUT, 'FitEarly') && ...
        isfield(OUT, 'FitLate') && isfield(OUT, 'signedStimEarly') && ...
        isfield(OUT, 'signedStimLate') && isfield(OUT, 'varStimEarly') && ...
        isfield(OUT, 'varStimLate') && isfield(OUT, 'RFrange'), ...
        'Cached OUT is missing required fields.');

    StimulusTable = OUT.StimulusTable;
    centroidXY = [StimulusTable.attendedCentroidX StimulusTable.attendedCentroidY];
    FitEarly = OUT.FitEarly;
    FitLate = OUT.FitLate;
    signedStimEarly = OUT.signedStimEarly;
    signedStimLate = OUT.signedStimLate;
    varStimEarly = OUT.varStimEarly;
    varStimLate = OUT.varStimLate;
    RFrange = OUT.RFrange(:);
    nIT = numel(RFrange);
else
    %% Load geometry
    Sgeo = load(tallPath);
    assert(isfield(Sgeo, 'Tall_IT') && isstruct(Sgeo.Tall_IT), ...
        '%s must contain struct Tall_IT.', tallFile);
    assert(isfield(Sgeo, 'ALLCOORDS') && isfield(Sgeo, 'RFrange'), ...
        '%s must contain ALLCOORDS and RFrange.', tallFile);

    Tall_IT = Sgeo.Tall_IT;
    ALLCOORDS = Sgeo.ALLCOORDS;
    RFrange = Sgeo.RFrange(:);
    nIT = numel(RFrange);
    siteRows = (1:nIT).';

    stimNums = arrayfun(@(x) x.stimNum, Tall_IT(:));
    [stimNumsSorted, ordStim] = sort(stimNums(:));
    assert(all(stimNumsSorted(:).' == 1:numel(Tall_IT)), ...
        'Tall_IT.stimNum must cover 1..%d exactly.', numel(Tall_IT));
    Tall_IT = Tall_IT(ordStim);
    nStim = numel(Tall_IT);

    StimulusTable = build_stimulus_table(ALLCOORDS, nStim);
    centroidXY = [StimulusTable.attendedCentroidX StimulusTable.attendedCentroidY];

    %% Load 3-bin responses and localize to IT rows
    R3_full = load_capsules_struct_exclusion_aware(resp3binPath, monkeySuffix, 'cfg', cfg);
    R3 = R3_full;
    R3.meanAct = R3_full.meanAct(RFrange, :, :);
    R3.meanSqAct = R3_full.meanSqAct(RFrange, :, :);
    if ismatrix(R3_full.nTrials) && size(R3_full.nTrials, 1) >= max(RFrange)
        R3.nTrials = R3_full.nTrials(RFrange, :);
    else
        R3.nTrials = R3_full.nTrials;
    end

    assert(size(R3.meanAct,1) == nIT, 'Localized IT response rows do not match Tall_IT.');
    assert(size(R3.meanAct,2) == nStim, 'Localized IT response stimuli do not match Tall_IT.');

    SNR = compute_snr_per_color_sites(R3, Tall_IT, siteRows, 'Verbose', false);
    muSpont = SNR.muSpont(siteRows);
    sdSpont = SNR.sdSpont(siteRows);

    [signedStimEarly, signedStimLate, varStimEarly, varStimLate] = ...
        compute_signed_stim_responses_with_var(R3, muSpont, sdSpont);

    %% Fit signed 2D attended-centroid tuning
    fields = {'status','nObs','baseline','amplitude','muX','muY','sigmaX','sigmaY', ...
        'thetaDeg','rmse','r2TrainPct','effectRangeObs','effectAbsPeakObs', ...
        'fStatApprox','pValueApprox','dfModel','dfError','rssNullWeighted', ...
        'rssFullWeighted'};
    FitEarly = init_fit_struct(nIT, fields);
    FitLate = init_fit_struct(nIT, fields);

    fprintf('Fitting IT attended-centroid tuning for monkey %s (%d sites, %d stimuli)\n', ...
        char(monkeySuffix), nIT, nStim);
    tFit = tic;

    for iSite = 1:nIT
        yEarly = signedStimEarly(iSite, :).';
        yLate = signedStimLate(iSite, :).';
        vEarly = varStimEarly(iSite, :).';
        vLate = varStimLate(iSite, :).';

        FitEarly(iSite) = fit_signed_gaussian2d_weighted(centroidXY, yEarly, vEarly, P);
        FitLate(iSite) = fit_signed_gaussian2d_weighted(centroidXY, yLate, vLate, P);

        if P.progressEvery > 0 && (iSite == 1 || mod(iSite, P.progressEvery) == 0 || iSite == nIT)
            elapsedSec = toc(tFit);
            rate = iSite / max(elapsedSec, eps);
            etaSec = (nIT - iSite) / max(rate, eps);
            fprintf('  Site %d / %d | elapsed %.1fs | ETA %.1fs\n', ...
                iSite, nIT, elapsedSec, etaSec);
        end
    end
end

veEarly = [FitEarly.r2TrainPct].';
veLate = [FitLate.r2TrainPct].';
effEarly = [FitEarly.effectRangeObs].';
effLate = [FitLate.effectRangeObs].';
pEarly = [FitEarly.pValueApprox].';
pLate = [FitLate.pValueApprox].';
isAttendedCentroidTunedEarly = isfinite(pEarly) & (pEarly < P.sigAlpha);
isAttendedCentroidTunedLate = isfinite(pLate) & (pLate < P.sigAlpha);

fprintf('Usable IT sites for attended-centroid fitting (early): %d / %d\n', ...
    nnz(isfinite(veEarly)), nIT);
fprintf('Usable IT sites for attended-centroid fitting (late):  %d / %d\n', ...
    nnz(isfinite(veLate)), nIT);
if any(isfinite(veEarly))
    fprintf('Median attended-centroid variance explained early: %.2f%%\n', median(veEarly(isfinite(veEarly))));
end
if any(isfinite(veLate))
    fprintf('Median attended-centroid variance explained late:  %.2f%%\n', median(veLate(isfinite(veLate))));
end
fprintf('Attended-centroid tuned sites at p < %.3f (early): %d / %d\n', ...
    P.sigAlpha, nnz(isAttendedCentroidTunedEarly), nIT);
fprintf('Attended-centroid tuned sites at p < %.3f (late):  %d / %d\n', ...
    P.sigAlpha, nnz(isAttendedCentroidTunedLate), nIT);

%% Summary figures
if P.makeSummaryFigures
    vePlotEarly = veEarly;
    vePlotLate = veLate;
    nClipEarly = nnz(isfinite(vePlotEarly) & (vePlotEarly < P.vePlotFloorPct));
    nClipLate = nnz(isfinite(vePlotLate) & (vePlotLate < P.vePlotFloorPct));
    vePlotEarly(isfinite(vePlotEarly) & (vePlotEarly < P.vePlotFloorPct)) = P.vePlotFloorPct;
    vePlotLate(isfinite(vePlotLate) & (vePlotLate < P.vePlotFloorPct)) = P.vePlotFloorPct;

    fprintf('Attended-centroid variance explained clipped for display below %.1f%%: early=%d, late=%d\n', ...
        P.vePlotFloorPct, nClipEarly, nClipLate);

    figure('Color', 'w');
    histogram(vePlotEarly(isfinite(vePlotEarly)), 30, 'FaceColor', [0.25 0.45 0.85], 'EdgeColor', 'none');
    hold on;
    histogram(vePlotLate(isfinite(vePlotLate)), 30, 'FaceColor', [0.85 0.35 0.25], 'EdgeColor', 'none');
    xlabel('Variance explained (%)');
    ylabel('N sites');
    title(sprintf('IT attended-centroid tuning variance explained (%s, floor %.0f%%)', ...
        char(monkeySuffix), P.vePlotFloorPct));
    legend('Early', 'Late');
    grid on;

    figure('Color', 'w');
    histogram(effEarly(isfinite(effEarly)), 30, 'FaceColor', [0.25 0.45 0.85], 'EdgeColor', 'none');
    hold on;
    histogram(effLate(isfinite(effLate)), 30, 'FaceColor', [0.85 0.35 0.25], 'EdgeColor', 'none');
    xlabel('Fitted modulation range over observed attended centroids (spont. SD units)');
    ylabel('N sites');
    title(sprintf('IT attended-centroid tuning effect size (%s)', char(monkeySuffix)));
    legend('Early', 'Late');
    grid on;

    figure('Color', 'w');
    scatter([FitLate.muX], [FitLate.muY], 28, veLate, 'filled');
    xlabel('Preferred attended-centroid x');
    ylabel('Preferred attended-centroid y');
    title(sprintf('IT late preferred attended-centroid position (color = %% variance explained, %s)', ...
        char(monkeySuffix)));
    colorbar;
    apply_square_limits(gca, centroidXY);
    grid on;
end

%% Example fits
if P.makeExampleFigure
    make_example_fit_figure(centroidXY, signedStimEarly, FitEarly, ...
        sprintf('IT attended-centroid fits early (%s)', char(monkeySuffix)), P.nExampleSites, RFrange);
    make_example_fit_figure(centroidXY, signedStimLate, FitLate, ...
        sprintf('IT attended-centroid fits late (%s)', char(monkeySuffix)), P.nExampleSites, RFrange);
end

%% Top-site tables
Tearly = build_fit_table(FitEarly, RFrange, 'r2TrainPct');
Tlate = build_fit_table(FitLate, RFrange, 'r2TrainPct');
Tearly.isAttendedCentroidTuned = isAttendedCentroidTunedEarly(Tearly.localSite);
Tlate.isAttendedCentroidTuned = isAttendedCentroidTunedLate(Tlate.localSite);

disp('Top IT attended-centroid sites by variance explained (early):');
disp(Tearly(1:min(10,height(Tearly)), :));
disp('Top IT attended-centroid sites by variance explained (late):');
disp(Tlate(1:min(10,height(Tlate)), :));

%% Pack output
OUT.P = P;
OUT.monkeySuffix = monkeySuffix;
OUT.RFrange = RFrange;
OUT.StimulusTable = StimulusTable;
OUT.signedStimEarly = signedStimEarly;
OUT.signedStimLate = signedStimLate;
OUT.varStimEarly = varStimEarly;
OUT.varStimLate = varStimLate;
OUT.FitEarly = FitEarly;
OUT.FitLate = FitLate;
OUT.pValueApproxEarly = pEarly;
OUT.pValueApproxLate = pLate;
OUT.isAttendedCentroidTunedEarly = isAttendedCentroidTunedEarly;
OUT.isAttendedCentroidTunedLate = isAttendedCentroidTunedLate;
OUT.TableEarly = Tearly;
OUT.TableLate = Tlate;

if needsSave && P.saveResult
    save(outPath, 'OUT', '-v7.3');
    fprintf('Saved IT attended-centroid tuning to %s\n', outPath);
end
end

function StimulusTable = build_stimulus_table(ALLCOORDS, nStim)
stimNum = (1:nStim).';
cueXY = nan(nStim, 2);
tFigXY = nan(nStim, 2);
centroidXY = nan(nStim, 2);

for stimIdx = 1:nStim
    [cueXY(stimIdx,:), tFigXY(stimIdx,:), centroidXY(stimIdx,:)] = ...
        attended_centroid_from_stim(ALLCOORDS, stimIdx);
end

StimulusTable = table(stimNum, cueXY(:,1), cueXY(:,2), tFigXY(:,1), tFigXY(:,2), ...
    centroidXY(:,1), centroidXY(:,2), 'VariableNames', ...
    {'stimNum','cueX','cueY','tFigX','tFigY','attendedCentroidX','attendedCentroidY'});
end

function [cue, tFig, centroid] = attended_centroid_from_stim(ALLCOORDS, stimNum)
f = sprintf('stim_%d', stimNum);
cue = double(ALLCOORDS.(f).s(:))';
tFig = double(ALLCOORDS.(f).t_fig(:))';
centroid = 0.5 * (cue + tFig);
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

function F = fit_signed_gaussian2d_weighted(X, y, varY, P)
fields = {'status','nObs','baseline','amplitude','muX','muY','sigmaX','sigmaY', ...
    'thetaDeg','rmse','r2TrainPct','effectRangeObs','effectAbsPeakObs', ...
    'fStatApprox','pValueApprox','dfModel','dfError','rssNullWeighted', ...
    'rssFullWeighted'};
F = struct();
for i = 1:numel(fields)
    F.(fields{i}) = NaN;
end
F.status = "not_fit";

valid = all(isfinite(X), 2) & isfinite(y) & isfinite(varY) & (varY > 0);
Xv = X(valid, :);
yv = y(valid);
vv = varY(valid);
F.nObs = numel(yv);

if numel(yv) < P.minStimuli
    F.status = "too_few_stimuli";
    return;
end

[pHat, w, ok] = fit_signed_gaussian2d_weighted_core(Xv, yv, vv, P.varFloorFrac);
if ~ok
    F.status = "fit_failed";
    return;
end

yHat = gaussian2d_signed(pHat, Xv);
wMean = sum(w .* yv) / sum(w);
rssNull = sum(w .* (yv - wMean).^2);
rssFull = sum(w .* (yv - yHat).^2);
dfModel = 6;
dfError = numel(yv) - 7;

F.status = "ok";
F.baseline = pHat(1);
F.amplitude = pHat(2);
F.muX = pHat(3);
F.muY = pHat(4);
F.sigmaX = pHat(5);
F.sigmaY = pHat(6);
F.thetaDeg = rad2deg(wrap_to_pi_local(pHat(7)));
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

function [pBest, w, ok] = fit_signed_gaussian2d_weighted_core(X, y, varY, varFloorFrac)
w = [];
xMin = min(X(:,1)); xMax = max(X(:,1));
yMin = min(X(:,2)); yMax = max(X(:,2));
spanX = max(xMax - xMin, 1);
spanY = max(yMax - yMin, 1);

b0 = median(y);
[~, iMax] = max(y);
[~, iMin] = min(y);
rngY = max(range(y), 0.25);

varFloor = max(median(varY(isfinite(varY) & varY > 0)) * varFloorFrac, 1e-6);
varUse = max(varY, varFloor);
w = 1 ./ varUse;
sqrtW = sqrt(w);

lb = [min(y) - 3*rngY, -4*rngY, xMin - spanX, yMin - spanY, 0.25, 0.25, -pi];
ub = [max(y) + 3*rngY,  4*rngY, xMax + spanX, yMax + spanY, 4*spanX, 4*spanY,  pi];

starts = [
    b0, y(iMax) - b0, X(iMax,1), X(iMax,2), 0.5*spanX, 0.5*spanY, 0;
    b0, y(iMin) - b0, X(iMin,1), X(iMin,2), 0.5*spanX, 0.5*spanY, 0;
    mean(y), mean(y - mean(y)), mean(X(:,1)), mean(X(:,2)), 0.7*spanX, 0.7*spanY, 0
    ];

opts = optimoptions('lsqcurvefit', 'Display', 'off', 'MaxFunctionEvaluations', 1e4);
modfun = @(p, Xin) sqrtW .* gaussian2d_signed(p, Xin);
yFit = sqrtW .* y;

bestSse = inf;
pBest = nan(1,7);
ok = false;

for iStart = 1:size(starts,1)
    p0 = starts(iStart,:);
    p0 = min(max(p0, lb), ub);
    try
        [pHat, ~, residual, exitflag] = lsqcurvefit(modfun, p0, X, yFit, lb, ub, opts); %#ok<ASGLU>
        if exitflag <= 0
            continue;
        end
        sse = sum(residual.^2);
        if sse < bestSse
            bestSse = sse;
            pBest = pHat;
            ok = true;
        end
    catch
    end
end
end

function yHat = gaussian2d_signed(p, X)
b = p(1);
a = p(2);
muX = p(3);
muY = p(4);
sx = p(5);
sy = p(6);
th = p(7);

dx = X(:,1) - muX;
dy = X(:,2) - muY;
ct = cos(th);
st = sin(th);

xr =  ct .* dx + st .* dy;
yr = -st .* dx + ct .* dy;

g = exp(-0.5 * ((xr ./ sx).^2 + (yr ./ sy).^2));
yHat = b + a .* g;
end

function make_example_fit_figure(X, Ymat, FitStruct, figName, nExampleSites, globalSites)
ve = [FitStruct.r2TrainPct].';
pval = [FitStruct.pValueApprox].';
statusText = string({FitStruct.status}).';
ok = (statusText == "ok") & isfinite(ve);

ord = find(ok);
if isempty(ord)
    return;
end
[~, idxSort] = sort(ve(ord), 'descend');
siteSel = ord(idxSort(1:min(nExampleSites, numel(idxSort))));

nShow = numel(siteSel);
allY = Ymat(siteSel, :);
allY = allY(isfinite(allY));
if isempty(allY)
    cLim = [-1 1];
else
    cLim = [min(allY) max(allY)];
    if cLim(1) == cLim(2)
        cLim = cLim + 0.5 * [-1 1];
    end
end

figW = max(320 * nShow, 720);
figH = 430;
fig = figure('Name', figName, 'NumberTitle', 'off', 'Color', 'w', ...
    'Units', 'pixels', 'Position', [80 120 figW figH]);

marginL = 0.06;
marginR = 0.03;
marginB = 0.12;
marginT = 0.18;
gap = 0.035;
cbGap = 0.02;
cbW = 0.018;

availW = 1 - marginL - marginR - cbGap - cbW - (nShow - 1) * gap;
availH = 1 - marginB - marginT;
panelPx = min((availW * figW) / nShow, availH * figH);
panelW = panelPx / figW;
panelH = panelPx / figH;
y0 = marginB + 0.5 * (availH - panelH);
x0 = marginL;

axList = gobjects(nShow, 1);
panelPos = zeros(nShow, 4);

for i = 1:nShow
    s = siteSel(i);
    pos = [x0 + (i - 1) * (panelW + gap), y0, panelW, panelH];
    panelPos(i,:) = pos;
    ax = axes('Parent', fig, 'Units', 'normalized', 'Position', pos);
    axList(i) = ax;
    if isprop(ax, 'PositionConstraint')
        ax.PositionConstraint = 'innerposition';
    end

    y = Ymat(s, :).';
    valid = isfinite(y);
    Xv = X(valid, :);
    yv = y(valid);

    scatter(ax, Xv(:,1), Xv(:,2), 55, yv, 'filled');
    hold(ax, 'on');
    colormap(ax, parula);
    caxis(ax, cLim);

    xLim = [min(X(:,1)) max(X(:,1))];
    yLim = [min(X(:,2)) max(X(:,2))];
    [xx, yy] = meshgrid(linspace(xLim(1), xLim(2), 80), linspace(yLim(1), yLim(2), 80));
    zz = gaussian2d_signed([FitStruct(s).baseline, FitStruct(s).amplitude, ...
        FitStruct(s).muX, FitStruct(s).muY, FitStruct(s).sigmaX, ...
        FitStruct(s).sigmaY, deg2rad(FitStruct(s).thetaDeg)], [xx(:) yy(:)]);
    zz = reshape(zz, size(xx));
    contour(ax, xx, yy, zz, 6, 'k-', 'LineWidth', 1);

    plot(ax, FitStruct(s).muX, FitStruct(s).muY, 'kp', 'MarkerSize', 10, 'MarkerFaceColor', 'y');
    xlabel(ax, 'Attended centroid x');
    if i == 1
        ylabel(ax, 'Attended centroid y');
    else
        ylabel(ax, '');
    end
    apply_square_limits(ax, X);
    grid(ax, 'on');
    box(ax, 'on');

    annotation(fig, 'textbox', [pos(1), pos(2) + pos(4) + 0.015, pos(3), 0.055], ...
        'String', sprintf('site %d\nVE %.1f%% | p %.3g', globalSites(s), FitStruct(s).r2TrainPct, pval(s)), ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 10, 'Interpreter', 'none');
end

if ~isempty(axList)
    cb = colorbar(axList(end));
    cb.Label.String = 'Signed response';
    cb.Units = 'normalized';
    cb.Position = [panelPos(end,1) + panelPos(end,3) + cbGap, y0, cbW, panelH];
    axList(end).Position = panelPos(end,:);
end

if exist('sgtitle', 'file') == 2
    sgtitle(figName);
end
end

function apply_square_limits(ax, X)
xMin = min(X(:,1)); xMax = max(X(:,1));
yMin = min(X(:,2)); yMax = max(X(:,2));
xMid = (xMin + xMax) / 2;
yMid = (yMin + yMax) / 2;
halfRange = 0.5 * max(xMax - xMin, yMax - yMin);
if ~isfinite(halfRange) || halfRange <= 0
    halfRange = 1;
end
pad = 0.05 * halfRange;
halfRange = halfRange + pad;
xlim(ax, [xMid - halfRange, xMid + halfRange]);
ylim(ax, [yMid - halfRange, yMid + halfRange]);
set(ax, 'DataAspectRatio', [1 1 1]);
end

function T = build_fit_table(FitStruct, RFrange, sortField)
n = numel(FitStruct);
statusText = string({FitStruct.status}).';
T = table((1:n).', RFrange(:), [FitStruct.nObs].', [FitStruct.amplitude].', ...
    [FitStruct.effectRangeObs].', [FitStruct.r2TrainPct].', ...
    [FitStruct.muX].', [FitStruct.muY].', [FitStruct.sigmaX].', ...
    [FitStruct.sigmaY].', [FitStruct.thetaDeg].', [FitStruct.pValueApprox].', statusText, ...
    'VariableNames', {'localSite','globalSiteInR','nObs','amplitude', ...
    'effectRangeObs','r2TrainPct','muX','muY','sigmaX','sigmaY','thetaDeg','pValueApprox','status'});

ok = isfinite(T.(sortField));
T = sortrows(T(ok, :), sortField, 'descend');
end

function p = fcdf_local(x, df1, df2)
z = (df1 .* x) ./ (df1 .* x + df2);
z = min(max(z, 0), 1);
p = betainc(z, df1/2, df2/2);
end

function th = wrap_to_pi_local(th)
th = mod(th + pi, 2*pi) - pi;
end
