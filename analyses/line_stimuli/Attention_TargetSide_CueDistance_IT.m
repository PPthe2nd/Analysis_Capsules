function OUT = Attention_TargetSide_CueDistance_IT()
% ATTENTION_TARGETSIDE_CUEDISTANCE_IT
% Compare IT attention-direction tuning for quartets where the cue is RF-near
% versus RF-far.
%
% This is a companion analysis to Attention_TargetSide_Tuning_IT. It reuses
% the RF-independent, within-quartet attention-direction contrast:
%   deltaQuartet = mean(resp DirA) - mean(resp DirB)
% where DirA and DirB are the two opposite attended directions within a
% quartet, each averaged over the two color-swapped stimuli.
%
% For each site, quartets are split into cue-near and cue-far halves using a
% balanced site-wise median split on RF-to-cue distance (pixels). The same
% weighted first-harmonic model is then fit separately for near and far
% quartets, and the fitted modulation ranges are compared.

%% Settings
P = struct();
P.Monkey = 1;           % 1 = Nilson, 2 = Figaro
P.minQuartets = 20;     % per split
P.makeSummaryFigure = true;
P.saveResult = true;
P.forceRefit = false;
P.sigAlpha = 0.05;
P.varFloorFrac = 1e-3;

cfg = config();

if P.Monkey == 1
    monkeySuffix = "N";
    tallFile = 'Tall_IT_lines_N.mat';
elseif P.Monkey == 2
    monkeySuffix = "F";
    tallFile = 'Tall_IT_lines_F.mat';
else
    error('Attention_TargetSide_CueDistance_IT:InvalidMonkey', ...
        'P.Monkey must be 1 (Nilson) or 2 (Figaro).');
end

basePath = fullfile(cfg.matDir, sprintf('Attention_TargetSide_Tuning_IT_directiondelta_%s.mat', char(monkeySuffix)));
tallPath = fullfile(cfg.matDir, tallFile);
outPath = fullfile(cfg.matDir, sprintf('Attention_TargetSide_CueDistance_IT_%s.mat', char(monkeySuffix)));

assert(exist(basePath, 'file') == 2, ...
    'Missing %s. Run Attention_TargetSide_Tuning_IT first.', basePath);
assert(exist(tallPath, 'file') == 2, ...
    'Missing %s. Run Line_Stimuli_IT first.', tallPath);

useCached = exist(outPath, 'file') == 2 && ~P.forceRefit;
needsSave = ~useCached;

if useCached
    fprintf('Loading cached IT cue-distance split from %s\n', outPath);
    S = load(outPath);
    assert(isfield(S, 'OUT') && isstruct(S.OUT), '%s must contain OUT.', outPath);
    OUT = S.OUT;
    assert(isfield(OUT, 'FitNearLate') && isfield(OUT, 'FitFarLate') && ...
        isfield(OUT, 'cueDistPx') && isfield(OUT, 'isNear'), ...
        'Cached OUT is missing required fields.');
    return;
end

%% Load base target-direction analysis
Sbase = load(basePath, 'OUT');
assert(isfield(Sbase, 'OUT') && isstruct(Sbase.OUT), '%s must contain OUT.', basePath);
BASE = Sbase.OUT;
requiredBase = {'QuartetTable','deltaQuartetEarly','deltaQuartetLate','varQuartetEarly', ...
    'varQuartetLate','FitEarly','FitLate','RFrange'};
for i = 1:numel(requiredBase)
    assert(isfield(BASE, requiredBase{i}), 'Base OUT is missing field %s.', requiredBase{i});
end

QuartetTable = BASE.QuartetTable;
thetaDeg = double(QuartetTable.targetDirDeg(:));
stimRef = double(QuartetTable.stimRef(:));
deltaQuartetEarly = BASE.deltaQuartetEarly;
deltaQuartetLate = BASE.deltaQuartetLate;
varQuartetEarly = BASE.varQuartetEarly;
varQuartetLate = BASE.varQuartetLate;
RFrange = BASE.RFrange(:);
nIT = numel(RFrange);
nQuartets = height(QuartetTable);

%% Load Tall_IT and cue distances
Stall = load(tallPath, 'Tall_IT');
assert(isfield(Stall, 'Tall_IT') && isstruct(Stall.Tall_IT), '%s must contain Tall_IT.', tallPath);
Tall_IT = Stall.Tall_IT;

cueDistPx = nan(nIT, nQuartets);
for q = 1:nQuartets
    Tq = Tall_IT(stimRef(q)).T;
    assert(height(Tq) >= nIT, 'Tall_IT(%d).T has fewer than %d IT rows.', stimRef(q), nIT);
    cueDistPx(:, q) = double(Tq.r_s(1:nIT));
end

isNear = false(nIT, nQuartets);
isFar = false(nIT, nQuartets);
for iSite = 1:nIT
    [isNear(iSite,:), isFar(iSite,:)] = balanced_median_split_row(cueDistPx(iSite,:));
end

fprintf('Cue-distance split: median near quartets/site = %.1f, median far quartets/site = %.1f\n', ...
    median(sum(isNear, 2)), median(sum(isFar, 2)));

%% Fit near/far separately
fields = {'status','nObs','baseline','coefCos','coefSin','amplitude','prefDeg', ...
    'rmse','r2TrainPct','effectRangeObs','effectAbsPeakObs','fStatApprox', ...
    'pValueApprox','dfModel','dfError','rssNullWeighted','rssFullWeighted','modelKind'};
FitNearEarly = init_fit_struct(nIT, fields, "direction");
FitNearLate = init_fit_struct(nIT, fields, "direction");
FitFarEarly = init_fit_struct(nIT, fields, "direction");
FitFarLate = init_fit_struct(nIT, fields, "direction");

for iSite = 1:nIT
    idxNear = isNear(iSite,:).';
    idxFar = isFar(iSite,:).';

    FitNearEarly(iSite) = fit_weighted_harmonic(thetaDeg(idxNear), ...
        deltaQuartetEarly(iSite, idxNear).', varQuartetEarly(iSite, idxNear).', 1, "direction", P);
    FitNearLate(iSite) = fit_weighted_harmonic(thetaDeg(idxNear), ...
        deltaQuartetLate(iSite, idxNear).', varQuartetLate(iSite, idxNear).', 1, "direction", P);
    FitFarEarly(iSite) = fit_weighted_harmonic(thetaDeg(idxFar), ...
        deltaQuartetEarly(iSite, idxFar).', varQuartetEarly(iSite, idxFar).', 1, "direction", P);
    FitFarLate(iSite) = fit_weighted_harmonic(thetaDeg(idxFar), ...
        deltaQuartetLate(iSite, idxFar).', varQuartetLate(iSite, idxFar).', 1, "direction", P);
end

effNearEarly = [FitNearEarly.effectRangeObs].';
effNearLate = [FitNearLate.effectRangeObs].';
effFarEarly = [FitFarEarly.effectRangeObs].';
effFarLate = [FitFarLate.effectRangeObs].';

okEarly = isfinite(effNearEarly) & isfinite(effFarEarly);
okLate = isfinite(effNearLate) & isfinite(effFarLate);

[pEarly, medDiffEarly] = paired_signrank_local(effNearEarly(okEarly), effFarEarly(okEarly));
[pLate, medDiffLate] = paired_signrank_local(effNearLate(okLate), effFarLate(okLate));

fprintf('Usable near/far effect-size pairs | early=%d late=%d\n', nnz(okEarly), nnz(okLate));
if any(okEarly)
    fprintf('Early effect size | near median=%.3f far median=%.3f near-far median diff=%.3f p=%.3g\n', ...
        median(effNearEarly(okEarly)), median(effFarEarly(okEarly)), medDiffEarly, pEarly);
end
if any(okLate)
    fprintf('Late effect size  | near median=%.3f far median=%.3f near-far median diff=%.3f p=%.3g\n', ...
        median(effNearLate(okLate)), median(effFarLate(okLate)), medDiffLate, pLate);
end

if P.makeSummaryFigure
    figure('Color', 'w');
    useTiled = exist('tiledlayout', 'file') == 2;
    if useTiled
        tiledlayout(2,2, 'TileSpacing', 'compact', 'Padding', 'compact');
    end

    if useTiled, nexttile; else, subplot(2,2,1); end
    scatter_near_far(effNearEarly, effFarEarly, sprintf('Early near vs far (%s)', char(monkeySuffix)));

    if useTiled, nexttile; else, subplot(2,2,2); end
    scatter_near_far(effNearLate, effFarLate, sprintf('Late near vs far (%s)', char(monkeySuffix)));

    if useTiled, nexttile; else, subplot(2,2,3); end
    histogram_diff(effNearEarly, effFarEarly, ...
        sprintf('Early near-far diff | p=%.3g', pEarly));

    if useTiled, nexttile; else, subplot(2,2,4); end
    histogram_diff(effNearLate, effFarLate, ...
        sprintf('Late near-far diff | p=%.3g', pLate));
end

OUT.P = P;
OUT.monkeySuffix = monkeySuffix;
OUT.RFrange = RFrange;
OUT.QuartetTable = QuartetTable;
OUT.cueDistPx = cueDistPx;
OUT.isNear = isNear;
OUT.isFar = isFar;
OUT.FitNearEarly = FitNearEarly;
OUT.FitNearLate = FitNearLate;
OUT.FitFarEarly = FitFarEarly;
OUT.FitFarLate = FitFarLate;
OUT.okEarly = okEarly;
OUT.okLate = okLate;
OUT.pEarly = pEarly;
OUT.pLate = pLate;
OUT.medianDiffEarly = medDiffEarly;
OUT.medianDiffLate = medDiffLate;

if needsSave && P.saveResult
    save(outPath, 'OUT', '-v7.3');
    fprintf('Saved IT cue-distance split analysis to %s\n', outPath);
end
end

function [isNear, isFar] = balanced_median_split_row(d)
d = double(d(:));
isNear = false(size(d));
isFar = false(size(d));
ok = isfinite(d);
idx = find(ok);
if isempty(idx)
    return;
end
[~, ord] = sort(d(ok), 'ascend');
n = numel(idx);
nNear = floor(n / 2);
nearIdx = idx(ord(1:nNear));
farIdx = idx(ord(nNear+1:end));
isNear(nearIdx) = true;
isFar(farIdx) = true;
isNear = isNear(:).';
isFar = isFar(:).';
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
F.prefDeg = mod(atan2d(beta(3), beta(2)), 360);
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

function scatter_near_far(nearVals, farVals, ttl)
ok = isfinite(nearVals) & isfinite(farVals);
nearVals = nearVals(ok);
farVals = farVals(ok);
scatter(nearVals, farVals, 28, 'filled', ...
    'MarkerFaceColor', [0.25 0.45 0.85], 'MarkerFaceAlpha', 0.6);
hold on;
lims = [nearVals(:); farVals(:)];
if isempty(lims)
    lims = [0; 1];
end
lo = min(lims);
hi = max(lims);
if ~(isfinite(lo) && isfinite(hi)) || hi <= lo
    hi = lo + 1;
end
plot([lo hi], [lo hi], 'k--', 'LineWidth', 1);
xlim([lo hi]);
ylim([lo hi]);
xlabel('Near effect size');
ylabel('Far effect size');
title(ttl);
axis square;
grid on;
end

function histogram_diff(nearVals, farVals, ttl)
ok = isfinite(nearVals) & isfinite(farVals);
d = nearVals(ok) - farVals(ok);
if isempty(d)
    d = NaN;
end
histogram(d, 25, 'FaceColor', [0.85 0.35 0.25], 'EdgeColor', 'none');
hold on;
xline(0, 'k--', 'LineWidth', 1);
xlabel('Near - far effect size');
ylabel('N sites');
title(ttl);
grid on;
end

function [pVal, medDiff] = paired_signrank_local(x, y)
medDiff = NaN;
pVal = NaN;
if isempty(x) || isempty(y)
    return;
end
d = x(:) - y(:);
medDiff = median(d);
if exist('signrank', 'file') == 2
    pVal = signrank(x, y);
end
end

function p = fcdf_local(x, df1, df2)
z = (df1 .* x) ./ (df1 .* x + df2);
z = min(max(z, 0), 1);
p = betainc(z, df1/2, df2/2);
end
