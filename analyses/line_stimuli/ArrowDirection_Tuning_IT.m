function OUT = ArrowDirection_Tuning_IT()
% ARROWDIRECTION_TUNING_IT
% Fit signed IT tuning to the overall arrow direction of each quartet.
%
% The response for each site and quartet is the signed trial-weighted mean
% response across the 4 stimuli in that quartet, normalized as:
%   (muQuartet - muSpont) / sdSpont
%
% Two weighted harmonic models are fit separately in the early and late
% windows:
%   direction: y = b + c*cos(theta) + s*sin(theta)
%   axis:      y = b + c*cos(2*theta) + s*sin(2*theta)
%
% Output metrics include:
% - weighted in-sample variance explained (%)
% - approximate weighted F-test p-value versus a constant model
% - preferred direction / preferred axis and tuning amplitude

%% Settings
P = struct();
P.Monkey = 1;           % 1 = Nilson, 2 = Figaro
P.minQuartets = 20;     % minimum number of finite quartet responses required per site
P.makeSummaryFigures = true;
P.makeExampleFigures = true;
P.makeAngleReferenceFigure = true;
P.nExampleSites = 3;    % per model class
P.saveResult = true;
P.forceRefit = false;
P.progressEvery = 10;
P.vePlotFloorPct = -100;
P.sigAlpha = 0.05;      % per-site significance threshold on the approximate p-value
P.varFloorFrac = 1e-3;  % floor on quartet-response variance relative to the median

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
    error('ArrowDirection_Tuning_IT:InvalidMonkey', ...
        'P.Monkey must be 1 (Nilson) or 2 (Figaro).');
end

tallPath = fullfile(cfg.matDir, tallFile);
resp3binPath = fullfile(cfg.matDir, resp3binFile);
outPath = fullfile(cfg.matDir, sprintf('ArrowDirection_Tuning_IT_%s.mat', char(monkeySuffix)));

assert(exist(tallPath, 'file') == 2, ...
    'Missing %s. Run Line_Stimuli_IT.m first.', tallPath);
assert(exist(resp3binPath, 'file') == 2, ...
    'Missing %s. Create the 3-bin response summary first.', resp3binPath);

useCached = exist(outPath, 'file') == 2 && ~P.forceRefit;
needsSave = ~useCached;

if useCached
    fprintf('Loading cached IT arrow-direction tuning from %s\n', outPath);
    S = load(outPath);
    assert(isfield(S, 'OUT') && isstruct(S.OUT), ...
        '%s must contain struct OUT.', outPath);
    OUT = S.OUT;
    assert(isfield(OUT, 'QuartetTable') && isfield(OUT, 'FitDirectionEarly') && ...
        isfield(OUT, 'FitDirectionLate') && isfield(OUT, 'FitAxisEarly') && ...
        isfield(OUT, 'FitAxisLate') && isfield(OUT, 'signedQuartetEarly') && ...
        isfield(OUT, 'signedQuartetLate') && isfield(OUT, 'RFrange'), ...
        'Cached OUT is missing required fields.');

    QuartetTable = OUT.QuartetTable;
    thetaDeg = QuartetTable.arrowDirDeg;
    FitDirectionEarly = OUT.FitDirectionEarly;
    FitDirectionLate = OUT.FitDirectionLate;
    FitAxisEarly = OUT.FitAxisEarly;
    FitAxisLate = OUT.FitAxisLate;
    signedQuartetEarly = OUT.signedQuartetEarly;
    signedQuartetLate = OUT.signedQuartetLate;
    varQuartetEarly = OUT.varQuartetEarly;
    varQuartetLate = OUT.varQuartetLate;
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

    [quartetMembers, thetaDeg, stimRef, cueXY] = build_quartet_members_and_arrow(ALLCOORDS, nStim);
    nQuartets = size(quartetMembers, 1);

    QuartetTable = table((1:nQuartets).', stimRef(:), cueXY(:,1), cueXY(:,2), thetaDeg(:), ...
        'VariableNames', {'quartetIdx','stimRef','cueX','cueY','arrowDirDeg'});

    %% Load 3-bin responses and localize to IT rows
    R3_full = load_capsules_struct_exclusion_aware(resp3binPath, monkeySuffix, 'cfg', cfg);
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

    %% Fit direction and axis models
    fields = {'status','nObs','baseline','coefCos','coefSin','amplitude','prefDeg', ...
        'rmse','r2TrainPct','effectRangeObs','effectAbsPeakObs','fStatApprox', ...
        'pValueApprox','dfModel','dfError','rssNullWeighted','rssFullWeighted','modelKind'};
    FitDirectionEarly = init_fit_struct(nIT, fields, "direction");
    FitDirectionLate = init_fit_struct(nIT, fields, "direction");
    FitAxisEarly = init_fit_struct(nIT, fields, "axis");
    FitAxisLate = init_fit_struct(nIT, fields, "axis");

    fprintf('Fitting IT arrow-direction tuning for monkey %s (%d sites, %d quartets)\n', ...
        char(monkeySuffix), nIT, nQuartets);
    tFit = tic;

    for iSite = 1:nIT
        yEarly = signedQuartetEarly(iSite, :).';
        yLate = signedQuartetLate(iSite, :).';
        vEarly = varQuartetEarly(iSite, :).';
        vLate = varQuartetLate(iSite, :).';

        FitDirectionEarly(iSite) = fit_weighted_harmonic(thetaDeg, yEarly, vEarly, 1, "direction", P);
        FitDirectionLate(iSite) = fit_weighted_harmonic(thetaDeg, yLate, vLate, 1, "direction", P);
        FitAxisEarly(iSite) = fit_weighted_harmonic(thetaDeg, yEarly, vEarly, 2, "axis", P);
        FitAxisLate(iSite) = fit_weighted_harmonic(thetaDeg, yLate, vLate, 2, "axis", P);

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

veDirEarly = [FitDirectionEarly.r2TrainPct].';
veDirLate = [FitDirectionLate.r2TrainPct].';
veAxisEarly = [FitAxisEarly.r2TrainPct].';
veAxisLate = [FitAxisLate.r2TrainPct].';
effDirEarly = [FitDirectionEarly.effectRangeObs].';
effDirLate = [FitDirectionLate.effectRangeObs].';
effAxisEarly = [FitAxisEarly.effectRangeObs].';
effAxisLate = [FitAxisLate.effectRangeObs].';
pDirEarly = [FitDirectionEarly.pValueApprox].';
pDirLate = [FitDirectionLate.pValueApprox].';
pAxisEarly = [FitAxisEarly.pValueApprox].';
pAxisLate = [FitAxisLate.pValueApprox].';

isDirEarly = isfinite(pDirEarly) & (pDirEarly < P.sigAlpha);
isDirLate = isfinite(pDirLate) & (pDirLate < P.sigAlpha);
isAxisEarly = isfinite(pAxisEarly) & (pAxisEarly < P.sigAlpha);
isAxisLate = isfinite(pAxisLate) & (pAxisLate < P.sigAlpha);

bestModelEarly = select_best_model(veDirEarly, veAxisEarly, isDirEarly, isAxisEarly);
bestModelLate = select_best_model(veDirLate, veAxisLate, isDirLate, isAxisLate);

fprintf('Usable IT sites for arrow-direction fitting (early): %d / %d\n', ...
    nnz(isfinite(veDirEarly) | isfinite(veAxisEarly)), nIT);
fprintf('Usable IT sites for arrow-direction fitting (late):  %d / %d\n', ...
    nnz(isfinite(veDirLate) | isfinite(veAxisLate)), nIT);
fprintf('Direction-tuned sites at p < %.3f | early=%d late=%d\n', ...
    P.sigAlpha, nnz(isDirEarly), nnz(isDirLate));
fprintf('Axis-tuned sites at p < %.3f | early=%d late=%d\n', ...
    P.sigAlpha, nnz(isAxisEarly), nnz(isAxisLate));
fprintf('Best significant model | early: direction=%d axis=%d | late: direction=%d axis=%d\n', ...
    nnz(bestModelEarly=="direction"), nnz(bestModelEarly=="axis"), ...
    nnz(bestModelLate=="direction"), nnz(bestModelLate=="axis"));

%% Summary figures
if P.makeSummaryFigures
    figure('Color', 'w');
    useTiled = exist('tiledlayout', 'file') == 2;
    if useTiled
        tiledlayout(2,2, 'TileSpacing', 'compact', 'Padding', 'compact');
    end

    if useTiled, nexttile; else, subplot(2,2,1); end
    plot_ve_histogram(veDirEarly, veDirLate, P.vePlotFloorPct, ...
        sprintf('Direction VE (%s)', char(monkeySuffix)));

    if useTiled, nexttile; else, subplot(2,2,2); end
    plot_ve_histogram(veAxisEarly, veAxisLate, P.vePlotFloorPct, ...
        sprintf('Axis VE (%s)', char(monkeySuffix)));

    if useTiled, nexttile; else, subplot(2,2,3); end
    histogram(mod([FitDirectionLate(isDirLate).prefDeg], 360), 18, ...
        'FaceColor', [0.20 0.55 0.85], 'EdgeColor', 'none');
    xlabel('Preferred direction (deg)');
    ylabel('N sites');
    title(sprintf('Late significant direction sites (%s)', char(monkeySuffix)));
    xlim([0 360]);
    xticks([0 90 180 270 360]);
    xticklabels({'0=R','90=U','180=L','270=D','360=R'});
    grid on;

    if useTiled, nexttile; else, subplot(2,2,4); end
    histogram(mod([FitAxisLate(isAxisLate).prefDeg], 180), 18, ...
        'FaceColor', [0.85 0.45 0.20], 'EdgeColor', 'none');
    xlabel('Preferred axis (deg)');
    ylabel('N sites');
    title(sprintf('Late significant axis sites (%s)', char(monkeySuffix)));
    xlim([0 180]);
    xticks([0 45 90 135 180]);
    xticklabels({'0=H','45','90=V','135','180=H'});
    grid on;

    figure('Color', 'w');
    useTiled = exist('tiledlayout', 'file') == 2;
    if useTiled
        tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');
    end

    if useTiled, nexttile; else, subplot(1,2,1); end
    plot_effect_histogram(effDirEarly, effDirLate, ...
        sprintf('Direction effect size (%s)', char(monkeySuffix)));

    if useTiled, nexttile; else, subplot(1,2,2); end
    plot_effect_histogram(effAxisEarly, effAxisLate, ...
        sprintf('Axis effect size (%s)', char(monkeySuffix)));
end

if P.makeAngleReferenceFigure
    make_angle_reference_figure(ALLCOORDS, RTAB384, QuartetTable, ...
        sprintf('IT arrow-angle reference (%s)', char(monkeySuffix)));
end

%% Example figures
if P.makeExampleFigures
    siteSelDirEarly = select_best_model_examples(bestModelEarly, veDirEarly, P.nExampleSites, "direction");
    siteSelAxisEarly = select_best_model_examples(bestModelEarly, veAxisEarly, P.nExampleSites, "axis");
    siteSelDir = select_best_model_examples(bestModelLate, veDirLate, P.nExampleSites, "direction");
    siteSelAxis = select_best_model_examples(bestModelLate, veAxisLate, P.nExampleSites, "axis");

    make_example_angle_fit_figure(thetaDeg, signedQuartetEarly, FitDirectionEarly, ...
        sprintf('IT arrow-direction fits early (%s)', char(monkeySuffix)), siteSelDirEarly, 1);
    make_example_angle_fit_figure(thetaDeg, signedQuartetEarly, FitAxisEarly, ...
        sprintf('IT arrow-axis fits early (%s)', char(monkeySuffix)), siteSelAxisEarly, 2);
    make_example_angle_fit_figure(thetaDeg, signedQuartetLate, FitDirectionLate, ...
        sprintf('IT arrow-direction fits late (%s)', char(monkeySuffix)), siteSelDir, 1);
    make_example_angle_fit_figure(thetaDeg, signedQuartetLate, FitAxisLate, ...
        sprintf('IT arrow-axis fits late (%s)', char(monkeySuffix)), siteSelAxis, 2);
end

%% Tables
TdirEarly = build_fit_table(FitDirectionEarly, RFrange, 'r2TrainPct');
TdirLate = build_fit_table(FitDirectionLate, RFrange, 'r2TrainPct');
TaxisEarly = build_fit_table(FitAxisEarly, RFrange, 'r2TrainPct');
TaxisLate = build_fit_table(FitAxisLate, RFrange, 'r2TrainPct');

disp('Top IT arrow-direction sites by variance explained (late):');
disp(TdirLate(1:min(10,height(TdirLate)), :));
disp('Top IT arrow-axis sites by variance explained (late):');
disp(TaxisLate(1:min(10,height(TaxisLate)), :));

%% Pack output
OUT.P = P;
OUT.monkeySuffix = monkeySuffix;
OUT.RFrange = RFrange;
OUT.QuartetTable = QuartetTable;
OUT.signedQuartetEarly = signedQuartetEarly;
OUT.signedQuartetLate = signedQuartetLate;
OUT.varQuartetEarly = varQuartetEarly;
OUT.varQuartetLate = varQuartetLate;
OUT.FitDirectionEarly = FitDirectionEarly;
OUT.FitDirectionLate = FitDirectionLate;
OUT.FitAxisEarly = FitAxisEarly;
OUT.FitAxisLate = FitAxisLate;
OUT.isDirectionTunedEarly = isDirEarly;
OUT.isDirectionTunedLate = isDirLate;
OUT.isAxisTunedEarly = isAxisEarly;
OUT.isAxisTunedLate = isAxisLate;
OUT.bestModelEarly = bestModelEarly;
OUT.bestModelLate = bestModelLate;
OUT.TableDirectionEarly = TdirEarly;
OUT.TableDirectionLate = TdirLate;
OUT.TableAxisEarly = TaxisEarly;
OUT.TableAxisLate = TaxisLate;

if needsSave && P.saveResult
    save(outPath, 'OUT', '-v7.3');
    fprintf('Saved IT arrow-direction tuning to %s\n', outPath);
end
end

function [quartetMembers, thetaDeg, stimRef, cueXY] = build_quartet_members_and_arrow(ALLCOORDS, nStim)
nQuartets = floor(nStim / 8) * 2;
quartetMembers = zeros(nQuartets, 4);
thetaDeg = nan(nQuartets, 1);
stimRef = nan(nQuartets, 1);
cueXY = nan(nQuartets, 2);
q = 0;

for base = 0:8:(nStim - 8)
    q = q + 1;
    quartetMembers(q,:) = base + [1 2 5 6];
    [thetaDeg(q), stimRef(q), cueXY(q,:)] = arrow_from_quartet(ALLCOORDS, quartetMembers(q,:));

    q = q + 1;
    quartetMembers(q,:) = base + [3 4 7 8];
    [thetaDeg(q), stimRef(q), cueXY(q,:)] = arrow_from_quartet(ALLCOORDS, quartetMembers(q,:));
end
end

function [thetaDeg, stimRef, cue] = arrow_from_quartet(ALLCOORDS, stimQ)
stimRef = stimQ(1);
[theta0, cue0] = arrow_direction_from_stim(ALLCOORDS, stimRef);

for k = 2:numel(stimQ)
    [thetaK, cueK] = arrow_direction_from_stim(ALLCOORDS, stimQ(k));
    assert(norm(cueK - cue0) < 1e-6, ...
        'Cue position differs within quartet [%s].', num2str(stimQ));
    assert(abs(angdiff_deg_local(thetaK, theta0)) < 1e-6, ...
        'Arrow direction differs within quartet [%s].', num2str(stimQ));
end

thetaDeg = mod(theta0, 360);
cue = cue0;
end

function [thetaDeg, cue] = arrow_direction_from_stim(ALLCOORDS, stimNum)
f = sprintf('stim_%d', stimNum);
s = double(ALLCOORDS.(f).s(:))';
tFig = double(ALLCOORDS.(f).t_fig(:))';
tBack = double(ALLCOORDS.(f).t_back(:))';

v1 = tFig - s;
v2 = tBack - s;
u1 = v1 / norm(v1);
u2 = v2 / norm(v2);
vBis = u1 + u2;

if norm(vBis) <= 0
    error('ArrowDirection_Tuning_IT:DegenerateBisector', ...
        'Degenerate bisector for stimulus %d.', stimNum);
end

thetaDeg = mod(atan2d(vBis(2), vBis(1)), 360);
cue = s;
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

function bestModel = select_best_model(veDir, veAxis, isDir, isAxis)
n = numel(veDir);
bestModel = strings(n,1);
bestModel(:) = "none";

idxDir = isDir & (~isAxis | (veDir >= veAxis));
idxAxis = isAxis & (~isDir | (veAxis > veDir));
bestModel(idxDir) = "direction";
bestModel(idxAxis) = "axis";
end

function siteSel = select_best_model_examples(bestModel, ve, nSelect, label)
idx = find(bestModel == label & isfinite(ve));
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

histogram(vePlotEarly(isfinite(vePlotEarly)), 30, 'FaceColor', [0.25 0.45 0.85], 'EdgeColor', 'none');
hold on;
histogram(vePlotLate(isfinite(vePlotLate)), 30, 'FaceColor', [0.85 0.35 0.25], 'EdgeColor', 'none');
xlabel('Variance explained (%)');
ylabel('N sites');
title(ttl);
legend('Early','Late');
grid on;
end

function plot_effect_histogram(effEarly, effLate, ttl)
histogram(effEarly(isfinite(effEarly)), 30, 'FaceColor', [0.25 0.45 0.85], 'EdgeColor', 'none');
hold on;
histogram(effLate(isfinite(effLate)), 30, 'FaceColor', [0.85 0.35 0.25], 'EdgeColor', 'none');
xlabel('Fitted modulation range (spont. SD units)');
ylabel('N sites');
title(ttl);
legend('Early','Late');
grid on;
end

function make_angle_reference_figure(ALLCOORDS, RTAB384, QuartetTable, figName)
dirTargets = [0 90 180 270];
axisTargets = [0 90];

theta = double(QuartetTable.arrowDirDeg(:));
stimRef = double(QuartetTable.stimRef(:));

dirStim = nan(size(dirTargets));
dirTheta = nan(size(dirTargets));
for i = 1:numel(dirTargets)
    d = abs(angdiff_deg_local(theta, dirTargets(i)));
    [~, idx] = min(d);
    dirStim(i) = stimRef(idx);
    dirTheta(i) = theta(idx);
end

axisStim = nan(size(axisTargets));
axisTheta = nan(size(axisTargets));
for i = 1:numel(axisTargets)
    d = abs(axisdiff_deg_local(theta, axisTargets(i)));
    [~, idx] = min(d);
    axisStim(i) = stimRef(idx);
    axisTheta(i) = theta(idx);
end

fig = figure('Name', figName, 'NumberTitle', 'off', 'Color', 'w');
useTiled = exist('tiledlayout', 'file') == 2;
if useTiled
    tl = tiledlayout(fig, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
end

for i = 1:4
    if useTiled, ax = nexttile(tl); else, ax = subplot(2,4,i); end
    img = render_stim_from_ALLCOORDS(ALLCOORDS, RTAB384, dirStim(i));
    imshow(img, 'Parent', ax, 'InitialMagnification', 'fit');
    axis(ax, 'image');
    title(ax, sprintf('%d= %s\nstim %d (%.1f°)', dirTargets(i), dir_label(dirTargets(i)), dirStim(i), dirTheta(i)), ...
        'FontSize', 10);
end

for i = 1:2
    if useTiled, ax = nexttile(tl); else, ax = subplot(2,4,4+i); end
    img = render_stim_from_ALLCOORDS(ALLCOORDS, RTAB384, axisStim(i));
    imshow(img, 'Parent', ax, 'InitialMagnification', 'fit');
    axis(ax, 'image');
    title(ax, sprintf('%d= %s axis\nstim %d (%.1f°)', axisTargets(i), axis_label(axisTargets(i)), axisStim(i), axisTheta(i)), ...
        'FontSize', 10);
end

for i = 3:4
    if useTiled, ax = nexttile(tl); else, ax = subplot(2,4,4+i); end
    axis(ax, 'off');
end

if exist('sgtitle', 'file') == 2
    sgtitle(figName);
end
end

function lbl = dir_label(thetaDeg)
switch mod(thetaDeg, 360)
    case 0
        lbl = 'R';
    case 90
        lbl = 'U';
    case 180
        lbl = 'L';
    case 270
        lbl = 'D';
    otherwise
        lbl = '?';
end
end

function lbl = axis_label(thetaDeg)
switch mod(thetaDeg, 180)
    case 0
        lbl = 'H';
    case 90
        lbl = 'V';
    otherwise
        lbl = '?';
end
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
    if harmonic == 1
        th = mod(thetaDeg(valid), 360);
    else
        th = mod(thetaDeg(valid), 180);
    end
    yv = y(valid);
    [thSort, ord] = sort(th);
    ySort = yv(ord);

    scatter(ax, thSort, ySort, 45, 'filled');
    hold(ax, 'on');

    fitY = harmonic_curve(thetaGrid, FitStruct(s), harmonic);
    plot(ax, thetaGrid, fitY, 'k-', 'LineWidth', 1.5);
    if harmonic == 1
        xlim(ax, [0 360]);
        xticks(ax, [0 90 180 270 360]);
        xticklabels(ax, {'0=R','90=U','180=L','270=D','360=R'});
        xlabel(ax, 'Arrow direction (deg)');
    else
        xlim(ax, [0 180]);
        xticks(ax, [0 45 90 135 180]);
        xticklabels(ax, {'0=H','45','90=V','135','180=H'});
        xlabel(ax, 'Arrow axis (deg)');
    end
    ylabel(ax, 'Signed response');
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

function d = axisdiff_deg_local(a, b)
d = mod((a - b) + 90, 180) - 90;
end
