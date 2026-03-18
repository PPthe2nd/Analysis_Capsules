function OUT = CuePosition_Tuning_IT()
% CUEPOSITION_TUNING_IT
% Fit signed cue-position tuning in IT from quartet-averaged responses.
%
% The response for each site and quartet is the signed trial-weighted mean
% response across the 4 stimuli in that quartet, normalized as:
%   (muQuartet - muSpont) / sdSpont
%
% A signed 2D elliptical Gaussian is then fit separately in the early and
% late windows:
%   y(x,y) = baseline + amplitude * G_2D(x,y)
%
% Output metrics include:
% - fitted effect size (observed fitted modulation range)
% - fitted amplitude (signed)
% - weighted in-sample variance explained (%)
% - approximate weighted F-test p-value versus a constant model

%% Settings
P = struct();
P.Monkey = 1;           % 1 = Nilson, 2 = Figaro
P.minQuartets = 20;     % minimum number of finite quartet responses required per site
P.makeSummaryFigures = true;
P.makeExampleFigure = true;
P.nExampleSites = 4;    % per window
P.makeEccentricExampleFigure = true;
P.nEccentricExampleSites = 3;
P.eccentricExampleMinVE = 25;
P.saveResult = true;
P.forceRefit = false;   % if false and a saved OUT exists, reload it instead of refitting
P.progressEvery = 10;   % print progress every N sites during fitting
P.vePlotFloorPct = -100; % histogram display floor; raw values stay unchanged
P.sigAlpha = 0.05;      % per-site significance threshold on the approximate p-value
P.varFloorFrac = 1e-3;  % floor on quartet-response variance relative to the median

cfg = config();

assert(exist('lsqcurvefit', 'file') == 2, ...
    'CuePosition_Tuning_IT requires lsqcurvefit (Optimization Toolbox).');

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
    error('CuePosition_Tuning_IT:InvalidMonkey', ...
        'P.Monkey must be 1 (Nilson) or 2 (Figaro).');
end

tallPath = fullfile(cfg.matDir, tallFile);
resp3binPath = fullfile(cfg.matDir, resp3binFile);
outPath = fullfile(cfg.matDir, sprintf('CuePosition_Tuning_IT_weighted_%s.mat', char(monkeySuffix)));

assert(exist(tallPath, 'file') == 2, ...
    'Missing %s. Run Line_Stimuli_IT.m first.', tallPath);
assert(exist(resp3binPath, 'file') == 2, ...
    'Missing %s. Create the 3-bin response summary first.', resp3binPath);

useCached = exist(outPath, 'file') == 2 && ~P.forceRefit;
needsSave = ~useCached;

if useCached
    fprintf('Loading cached IT cue-position tuning from %s\n', outPath);
    S = load(outPath);
    assert(isfield(S, 'OUT') && isstruct(S.OUT), ...
        '%s must contain struct OUT.', outPath);
    OUT = S.OUT;
    assert(isfield(OUT, 'QuartetTable') && isfield(OUT, 'FitEarly') && ...
        isfield(OUT, 'FitLate') && isfield(OUT, 'signedQuartetEarly') && ...
        isfield(OUT, 'signedQuartetLate') && isfield(OUT, 'varQuartetEarly') && ...
        isfield(OUT, 'varQuartetLate') && isfield(OUT, 'RFrange'), ...
        'Cached OUT is missing required fields.');

    QuartetTable = OUT.QuartetTable;
    cueXY = [QuartetTable.cueX QuartetTable.cueY];
    FitEarly = OUT.FitEarly;
    FitLate = OUT.FitLate;
    signedQuartetEarly = OUT.signedQuartetEarly;
    signedQuartetLate = OUT.signedQuartetLate;
    varQuartetEarly = OUT.varQuartetEarly;
    varQuartetLate = OUT.varQuartetLate;
    RFrange = OUT.RFrange(:);
    nIT = numel(RFrange);
    nQuartets = height(QuartetTable);
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

    [quartetMembers, cueXY, stimRef] = build_quartet_members_and_cues(ALLCOORDS, nStim);
    nQuartets = size(quartetMembers, 1);

    QuartetTable = table((1:nQuartets).', stimRef(:), cueXY(:,1), cueXY(:,2), ...
        'VariableNames', {'quartetIdx','stimRef','cueX','cueY'});

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

    [signedQuartetEarly, signedQuartetLate, varQuartetEarly, varQuartetLate] = ...
        compute_signed_quartet_responses_with_var(R3, quartetMembers, muSpont, sdSpont);

    %% Fit signed 2D cue-position tuning
    Fields = {'status','nObs','baseline','amplitude','muX','muY','sigmaX','sigmaY', ...
        'thetaDeg','rmse','r2TrainPct','effectRangeObs','effectAbsPeakObs', ...
        'fStatApprox','pValueApprox','dfModel','dfError','rssNullWeighted', ...
        'rssFullWeighted'};
    FitEarly = init_fit_struct(nIT, Fields);
    FitLate = init_fit_struct(nIT, Fields);

    fprintf('Fitting IT cue-position tuning for monkey %s (%d sites, %d quartets)\n', ...
        char(monkeySuffix), nIT, nQuartets);
    tFit = tic;

    for iSite = 1:nIT
        yEarly = signedQuartetEarly(iSite, :).';
        yLate = signedQuartetLate(iSite, :).';
        vEarly = varQuartetEarly(iSite, :).';
        vLate = varQuartetLate(iSite, :).';

        FitEarly(iSite) = fit_signed_gaussian2d_weighted(cueXY, yEarly, vEarly, P);
        FitLate(iSite) = fit_signed_gaussian2d_weighted(cueXY, yLate, vLate, P);

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
isCueTunedEarly = isfinite(pEarly) & (pEarly < P.sigAlpha);
isCueTunedLate = isfinite(pLate) & (pLate < P.sigAlpha);

fprintf('Usable IT sites for cue-position fitting (early): %d / %d\n', ...
    nnz(isfinite(veEarly)), nIT);
fprintf('Usable IT sites for cue-position fitting (late):  %d / %d\n', ...
    nnz(isfinite(veLate)), nIT);
if any(isfinite(veEarly))
    fprintf('Median cue-position variance explained early: %.2f%%\n', median(veEarly(isfinite(veEarly))));
end
if any(isfinite(veLate))
    fprintf('Median cue-position variance explained late:  %.2f%%\n', median(veLate(isfinite(veLate))));
end
fprintf('Cue-position tuned sites at p < %.3f (early): %d / %d\n', P.sigAlpha, nnz(isCueTunedEarly), nIT);
fprintf('Cue-position tuned sites at p < %.3f (late):  %d / %d\n', P.sigAlpha, nnz(isCueTunedLate), nIT);
needsSave = needsSave | ~useCached;

%% Summary figures
if P.makeSummaryFigures
    vePlotEarly = veEarly;
    vePlotLate = veLate;
    nClipEarly = nnz(isfinite(vePlotEarly) & (vePlotEarly < P.vePlotFloorPct));
    nClipLate = nnz(isfinite(vePlotLate) & (vePlotLate < P.vePlotFloorPct));
    vePlotEarly(isfinite(vePlotEarly) & (vePlotEarly < P.vePlotFloorPct)) = P.vePlotFloorPct;
    vePlotLate(isfinite(vePlotLate) & (vePlotLate < P.vePlotFloorPct)) = P.vePlotFloorPct;

    fprintf('Cue-position variance explained clipped for display below %.1f%%: early=%d, late=%d\n', ...
        P.vePlotFloorPct, nClipEarly, nClipLate);

    figure;
    histogram(vePlotEarly(isfinite(vePlotEarly)), 30, 'FaceColor', [0.25 0.45 0.85], 'EdgeColor', 'none');
    hold on;
    histogram(vePlotLate(isfinite(vePlotLate)), 30, 'FaceColor', [0.85 0.35 0.25], 'EdgeColor', 'none');
    xlabel('Variance explained (%)');
    ylabel('N sites');
    title(sprintf('IT cue-position tuning variance explained (%s, floor %.0f%%)', ...
        char(monkeySuffix), P.vePlotFloorPct));
    legend('Early','Late');
    grid on;

    figure;
    histogram(effEarly(isfinite(effEarly)), 30, 'FaceColor', [0.25 0.45 0.85], 'EdgeColor', 'none');
    hold on;
    histogram(effLate(isfinite(effLate)), 30, 'FaceColor', [0.85 0.35 0.25], 'EdgeColor', 'none');
    xlabel('Fitted modulation range over observed cue positions');
    ylabel('N sites');
    title(sprintf('IT cue-position tuning effect size (%s)', char(monkeySuffix)));
    legend('Early','Late');
    grid on;

    figure;
    scatter([FitLate.muX], [FitLate.muY], 28, veLate, 'filled');
    xlabel('Preferred cue x');
    ylabel('Preferred cue y');
    title(sprintf('IT late preferred cue position (color = %% variance explained, %s)', char(monkeySuffix)));
    colorbar;
    apply_square_limits(gca, cueXY);
    grid on;
end

%% Example fits
if P.makeExampleFigure
    make_example_fit_figure(cueXY, signedQuartetEarly, FitEarly, ...
        sprintf('IT cue-position fits early (%s)', char(monkeySuffix)), P.nExampleSites);
    make_example_fit_figure(cueXY, signedQuartetLate, FitLate, ...
        sprintf('IT cue-position fits late (%s)', char(monkeySuffix)), P.nExampleSites);
    if P.makeEccentricExampleFigure
        siteSelEcc = select_eccentric_example_sites(cueXY, FitLate, isCueTunedLate, ...
            veLate, P.eccentricExampleMinVE, P.nEccentricExampleSites);
        make_example_fit_figure(cueXY, signedQuartetLate, FitLate, ...
            sprintf('IT cue-position fits late eccentric-preferring (%s)', char(monkeySuffix)), ...
            P.nEccentricExampleSites, siteSelEcc);
    end
end

%% Top-site tables
Tearly = build_fit_table(FitEarly, RFrange, 'r2TrainPct');
Tlate = build_fit_table(FitLate, RFrange, 'r2TrainPct');
Tearly.pValueApprox = pEarly(Tearly.localSite);
Tearly.isCueTuned = isCueTunedEarly(Tearly.localSite);
Tlate.pValueApprox = pLate(Tlate.localSite);
Tlate.isCueTuned = isCueTunedLate(Tlate.localSite);

disp('Top IT cue-position sites by variance explained (early):');
disp(Tearly(1:min(10,height(Tearly)), :));
disp('Top IT cue-position sites by variance explained (late):');
disp(Tlate(1:min(10,height(Tlate)), :));

%% Pack output
OUT.P = P;
OUT.monkeySuffix = monkeySuffix;
OUT.RFrange = RFrange;
OUT.QuartetTable = QuartetTable;
OUT.signedQuartetEarly = signedQuartetEarly;
OUT.signedQuartetLate = signedQuartetLate;
OUT.varQuartetEarly = varQuartetEarly;
OUT.varQuartetLate = varQuartetLate;
OUT.FitEarly = FitEarly;
OUT.FitLate = FitLate;
OUT.pValueApproxEarly = pEarly;
OUT.pValueApproxLate = pLate;
OUT.isCueTunedEarly = isCueTunedEarly;
OUT.isCueTunedLate = isCueTunedLate;
OUT.TableEarly = Tearly;
OUT.TableLate = Tlate;

if needsSave && P.saveResult
    save(outPath, 'OUT', '-v7.3');
    fprintf('Saved IT cue-position tuning to %s\n', outPath);
end

end

function [quartetMembers, cueXY, stimRef] = build_quartet_members_and_cues(ALLCOORDS, nStim)
nQuartets = floor(nStim / 8) * 2;
quartetMembers = zeros(nQuartets, 4);
cueXY = nan(nQuartets, 2);
stimRef = nan(nQuartets, 1);
q = 0;

for base = 0:8:(nStim - 8)
    q = q + 1;
    quartetMembers(q,:) = base + [1 2 5 6];
    [cueXY(q,:), stimRef(q)] = cue_from_quartet(ALLCOORDS, quartetMembers(q,:));

    q = q + 1;
    quartetMembers(q,:) = base + [3 4 7 8];
    [cueXY(q,:), stimRef(q)] = cue_from_quartet(ALLCOORDS, quartetMembers(q,:));
end
end

function [cue, stimRef] = cue_from_quartet(ALLCOORDS, stimQ)
stimRef = stimQ(1);
f0 = sprintf('stim_%d', stimRef);
cue0 = double(ALLCOORDS.(f0).s(:))';
for k = 2:numel(stimQ)
    fk = sprintf('stim_%d', stimQ(k));
    cueK = double(ALLCOORDS.(fk).s(:))';
    assert(norm(cueK - cue0) < 1e-6, ...
        'Cue position differs within quartet [%s].', num2str(stimQ));
end
cue = cue0;
end

function [signedEarly, signedLate, varEarly, varLate] = compute_signed_quartet_responses_with_var(R3, quartetMembers, muSpont, sdSpont)
[nSites, ~, ~] = size(R3.meanAct);
nQuartets = size(quartetMembers, 1);

if isvector(R3.nTrials)
    nTrialsAll = double(R3.nTrials(:));
    perSiteTrials = false;
else
    perSiteTrials = true;
end

muQuartetEarly = nan(nSites, nQuartets);
muQuartetLate = nan(nSites, nQuartets);
varQuartetEarly = nan(nSites, nQuartets);
varQuartetLate = nan(nSites, nQuartets);

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
        stimQ = quartetMembers(qIdx,:);
        nTrQ = nTrSite(stimQ);

        rEarlyQ = rEarly(stimQ);
        rSqEarlyQ = rSqEarly(stimQ);
        idxEarly = isfinite(rEarlyQ) & isfinite(rSqEarlyQ) & isfinite(nTrQ) & (nTrQ > 1);
        if any(idxEarly)
            nUse = nTrQ(idxEarly);
            muUse = rEarlyQ(idxEarly);
            msqUse = rSqEarlyQ(idxEarly);
            varStim = max(0, msqUse - muUse.^2);
            varStim = varStim .* (nUse ./ max(nUse - 1, 1));
            muQuartetEarly(iSite, qIdx) = sum(nUse .* muUse) / sum(nUse);
            varQuartetEarly(iSite, qIdx) = sum(nUse .* varStim) / (sum(nUse)^2);
        end

        rLateQ = rLate(stimQ);
        rSqLateQ = rSqLate(stimQ);
        idxLate = isfinite(rLateQ) & isfinite(rSqLateQ) & isfinite(nTrQ) & (nTrQ > 1);
        if any(idxLate)
            nUse = nTrQ(idxLate);
            muUse = rLateQ(idxLate);
            msqUse = rSqLateQ(idxLate);
            varStim = max(0, msqUse - muUse.^2);
            varStim = varStim .* (nUse ./ max(nUse - 1, 1));
            muQuartetLate(iSite, qIdx) = sum(nUse .* muUse) / sum(nUse);
            varQuartetLate(iSite, qIdx) = sum(nUse .* varStim) / (sum(nUse)^2);
        end
    end
end

signedEarly = bsxfun(@rdivide, bsxfun(@minus, muQuartetEarly, muSpont), sdSpont);
signedLate = bsxfun(@rdivide, bsxfun(@minus, muQuartetLate, muSpont), sdSpont);
varEarly = bsxfun(@rdivide, varQuartetEarly, sdSpont.^2);
varLate = bsxfun(@rdivide, varQuartetLate, sdSpont.^2);

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
Fields = {'status','nObs','baseline','amplitude','muX','muY','sigmaX','sigmaY', ...
    'thetaDeg','rmse','r2TrainPct','effectRangeObs','effectAbsPeakObs', ...
    'fStatApprox','pValueApprox','dfModel','dfError','rssNullWeighted', ...
    'rssFullWeighted'};
F = struct();
for i = 1:numel(Fields)
    F.(Fields{i}) = NaN;
end
F.status = "not_fit";

valid = all(isfinite(X), 2) & isfinite(y) & isfinite(varY) & (varY > 0);
Xv = X(valid, :);
yv = y(valid);
vv = varY(valid);
F.nObs = numel(yv);

if numel(yv) < P.minQuartets
    F.status = "too_few_quartets";
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

function make_example_fit_figure(X, Ymat, FitStruct, figName, nExampleSites, siteSel)
ve = [FitStruct.r2TrainPct].';
pval = [FitStruct.pValueApprox].';
statusText = string({FitStruct.status}).';
ok = (statusText == "ok") & isfinite(ve);
if nargin < 6 || isempty(siteSel)
    ord = find(ok);
    if isempty(ord)
        return;
    end
    [~, idxSort] = sort(ve(ord), 'descend');
    siteSel = ord(idxSort(1:min(nExampleSites, numel(idxSort))));
else
    siteSel = siteSel(:);
    siteSel = siteSel(isfinite(siteSel) & siteSel >= 1 & siteSel <= numel(FitStruct));
    siteSel = unique(siteSel, 'stable');
    siteSel = siteSel(ok(siteSel));
    if isempty(siteSel)
        return;
    end
end

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
    hold on;
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
    xlabel(ax, 'Cue x');
    if i == 1
        ylabel(ax, 'Cue y');
    else
        ylabel(ax, '');
    end
    apply_square_limits(ax, X);
    grid(ax, 'on');
    box(ax, 'on');

    annotation(fig, 'textbox', [pos(1), pos(2) + pos(4) + 0.015, pos(3), 0.055], ...
        'String', sprintf('site %d\nVE %.1f%% | p %.3g', s, FitStruct(s).r2TrainPct, pval(s)), ...
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

function siteSel = select_eccentric_example_sites(X, FitStruct, isSig, ve, minVE, nSelect)
prefObsEcc = preferred_observed_eccentricity(X, FitStruct);
ok = strcmp(string({FitStruct.status}).', "ok") & isSig(:) & isfinite(ve(:)) & ...
    (ve(:) > minVE) & isfinite(prefObsEcc);
idx = find(ok);
if isempty(idx)
    siteSel = [];
    return;
end

sortMat = [prefObsEcc(idx), ve(idx)];
[~, ord] = sortrows(sortMat, [-1 -2]);
siteSel = idx(ord(1:min(nSelect, numel(ord))));
end

function prefObsEcc = preferred_observed_eccentricity(X, FitStruct)
nSites = numel(FitStruct);
prefObsEcc = nan(nSites,1);
ecc = hypot(X(:,1), X(:,2));
for iSite = 1:nSites
    if string(FitStruct(iSite).status) ~= "ok"
        continue;
    end
    yHat = gaussian2d_signed([FitStruct(iSite).baseline, FitStruct(iSite).amplitude, ...
        FitStruct(iSite).muX, FitStruct(iSite).muY, FitStruct(iSite).sigmaX, ...
        FitStruct(iSite).sigmaY, deg2rad(FitStruct(iSite).thetaDeg)], X);
    [~, iMax] = max(yHat);
    prefObsEcc(iSite) = ecc(iMax);
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
