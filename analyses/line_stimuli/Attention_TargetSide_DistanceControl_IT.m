function OUT = Attention_TargetSide_DistanceControl_IT(Puser)
% ATTENTION_TARGETSIDE_DISTANCECONTROL_IT
% Test whether the IT target-direction effect survives after accounting for
% RF proximity to the target versus distractor arm.
%
% Uses the same late within-quartet attention-direction contrast as
% Attention_TargetSide_Tuning_IT:
%   deltaQuartet = mean(resp DirA) - mean(resp DirB)
%
% For each site/quartet, this script adds two geometry controls:
%   targetAdvPx : positive when the target arm is closer to the RF than the
%                 distractor arm, aligned to DirA versus DirB
%   meanDistPx  : overall RF-to-object distance
%
% The main test is whether adding directional terms cos(phi), sin(phi)
% improves a weighted regression beyond distance-only predictors.

%% Settings
P = struct();
P.Monkey = 1;          % 1 = Nilson, 2 = Figaro
P.minQuartets = 20;
P.makeSummaryFigure = true;
P.saveResult = true;
P.forceRefit = false;
P.sigAlpha = 0.05;
P.varFloorFrac = 1e-3;

if nargin >= 1 && ~isempty(Puser)
    fn = fieldnames(Puser);
    for i = 1:numel(fn)
        P.(fn{i}) = Puser.(fn{i});
    end
end

cfg = config();

if P.Monkey == 1
    monkeySuffix = "N";
    tallFile = 'Tall_IT_lines_N.mat';
elseif P.Monkey == 2
    monkeySuffix = "F";
    tallFile = 'Tall_IT_lines_F.mat';
else
    error('Attention_TargetSide_DistanceControl_IT:InvalidMonkey', ...
        'P.Monkey must be 1 (Nilson) or 2 (Figaro).');
end

basePath = fullfile(cfg.matDir, sprintf('Attention_TargetSide_Tuning_IT_directiondelta_%s.mat', char(monkeySuffix)));
tallPath = fullfile(cfg.matDir, tallFile);
outPath = fullfile(cfg.matDir, sprintf('Attention_TargetSide_DistanceControl_IT_%s.mat', char(monkeySuffix)));
currentExclusions = site_session_exclusions(monkeySuffix);
hasSessionExclusions = ~isempty(currentExclusions);

assert(exist(tallPath, 'file') == 2, ...
    'Missing %s. Run Line_Stimuli_IT first.', tallPath);

useCached = exist(outPath, 'file') == 2 && ~P.forceRefit;
if useCached
    S = load(outPath, 'OUT');
    if ~session_exclusion_cache_matches(S, monkeySuffix)
        useCached = false;
        if hasSessionExclusions
            fprintf(['Cached IT distance-control analysis does not match the active session exclusions ' ...
                     'for monkey %s; recomputing.\n'], char(monkeySuffix));
        else
            fprintf('Cached IT distance-control analysis is from an exclusion-aware run; recomputing canonical cache.\n');
        end
    else
        fprintf('Loading cached IT distance-control analysis from %s\n', outPath);
        assert(isfield(S, 'OUT') && isstruct(S.OUT), '%s must contain OUT.', outPath);
        requiredCached = {'intPLate','intSig','distGivenDirPLate','uniqueR2Dir','uniqueR2Dist','sharedR2Additive','totalR2Additive','totalR2Interaction','RegressionLate'};
        hasNewFields = all(isfield(S.OUT, requiredCached));
        if hasNewFields
            Rtmp = S.OUT.RegressionLate;
            hasNewFields = ~isempty(Rtmp) && isfield(Rtmp, 'pInteraction') && ...
                isfield(Rtmp, 'partialR2Interaction') && isfield(Rtmp, 'pDistGivenDir') && ...
                isfield(Rtmp, 'uniqueR2Dir') && isfield(Rtmp, 'uniqueR2Dist') && isfield(Rtmp, 'sharedR2Additive');
        end
        if hasNewFields
            OUT = S.OUT;
            if P.makeSummaryFigure
                make_summary_figure(OUT);
            end
            return;
        end
        fprintf('Cached result is missing interaction fields; recomputing.\n');
    end
end

%% Load base analysis
if hasSessionExclusions
    fprintf(['Session exclusions are active for monkey %s; refreshing IT target-side ' ...
             'fits before distance-control regression.\n'], char(monkeySuffix));
    BASE = Attention_TargetSide_Tuning_IT(struct( ...
        'makeSummaryFigures', false, ...
        'makeExampleFigures', false, ...
        'makeAngleReferenceFigure', false));
else
    assert(exist(basePath, 'file') == 2, ...
        'Missing %s. Run Attention_TargetSide_Tuning_IT first.', basePath);
    Sbase = load(basePath, 'OUT');
    assert(isfield(Sbase, 'OUT') && isstruct(Sbase.OUT), '%s must contain OUT.', basePath);
    BASE = Sbase.OUT;
end
requiredBase = {'QuartetTable','deltaQuartetLate','varQuartetLate','FitLate','RFrange'};
for i = 1:numel(requiredBase)
    assert(isfield(BASE, requiredBase{i}), 'Base OUT is missing field %s.', requiredBase{i});
end

QuartetTable = BASE.QuartetTable;
thetaDeg = double(QuartetTable.targetDirDeg(:));
stimRef = double(QuartetTable.stimRef(:));
deltaQuartetLate = double(BASE.deltaQuartetLate);
varQuartetLate = double(BASE.varQuartetLate);
FitLateBase = BASE.FitLate;
RFrange = BASE.RFrange(:);
nIT = numel(RFrange);
nQuartets = height(QuartetTable);

%% Load geometry
Stall = load(tallPath, 'Tall_IT', 'ALLCOORDS', 'RTAB384');
assert(isfield(Stall, 'Tall_IT') && isstruct(Stall.Tall_IT), '%s must contain Tall_IT.', tallPath);
assert(isfield(Stall, 'ALLCOORDS') && isfield(Stall, 'RTAB384'), ...
    '%s must contain ALLCOORDS and RTAB384.', tallPath);
Tall_IT = Stall.Tall_IT;
ALLCOORDS = Stall.ALLCOORDS;
RTAB384 = Stall.RTAB384;

%% Distances from RF to target and distractor arms for every stimulus
[distTargetPx, distDistractorPx] = compute_arm_distances_all_stim(ALLCOORDS, RTAB384, Tall_IT, nIT);

%% Quartet-wise distance predictors aligned to DirA versus DirB
targetAdvPx = nan(nIT, nQuartets);
meanDistPx = nan(nIT, nQuartets);

for q = 1:nQuartets
    pairA = [stimRef(q), stimRef(q) + 4];
    pairB = [stimRef(q) + 1, stimRef(q) + 5];

    dTA = rowmean_omitnan(distTargetPx(:, pairA));
    dDA = rowmean_omitnan(distDistractorPx(:, pairA));
    dTB = rowmean_omitnan(distTargetPx(:, pairB));
    dDB = rowmean_omitnan(distDistractorPx(:, pairB));

    % Positive when DirA has more target-over-distractor proximity than DirB.
    targetAdvPx(:, q) = 0.5 * ((dDA - dTA) - (dDB - dTB));
    % Overall RF-to-object proximity, irrespective of target assignment.
    meanDistPx(:, q) = 0.25 * (dTA + dDA + dTB + dDB);
end

%% Site-wise weighted regressions
RegressionLate = init_regression_struct(nIT);
for iSite = 1:nIT
    RegressionLate(iSite) = fit_direction_distance_control( ...
        thetaDeg, ...
        deltaQuartetLate(iSite, :).', ...
        varQuartetLate(iSite, :).', ...
        targetAdvPx(iSite, :).', ...
        meanDistPx(iSite, :).', ...
        P);
end

origPLate = [FitLateBase.pValueApprox].';
dirOnlyPLate = [RegressionLate.pDirOnly].';
ctrlPLate = [RegressionLate.pDirGivenDist].';
distPLate = [RegressionLate.pDistOnly].';
distGivenDirPLate = [RegressionLate.pDistGivenDir].';
intPLate = [RegressionLate.pInteraction].';
partialR2Dir = [RegressionLate.partialR2Dir].';
partialR2DistGivenDir = [RegressionLate.partialR2DistGivenDir].';
partialR2Interaction = [RegressionLate.partialR2Interaction].';
uniqueR2Dir = [RegressionLate.uniqueR2Dir].';
uniqueR2Dist = [RegressionLate.uniqueR2Dist].';
sharedR2Additive = [RegressionLate.sharedR2Additive].';
totalR2Additive = [RegressionLate.r2Additive].';
totalR2Interaction = [RegressionLate.r2Additive].' + [RegressionLate.partialR2Interaction].' .* (1 - [RegressionLate.r2Additive].');

origSig = isfinite(origPLate) & (origPLate < P.sigAlpha);
dirOnlySig = isfinite(dirOnlyPLate) & (dirOnlyPLate < P.sigAlpha);
ctrlSig = isfinite(ctrlPLate) & (ctrlPLate < P.sigAlpha);
distSig = isfinite(distPLate) & (distPLate < P.sigAlpha);
distGivenDirSig = isfinite(distGivenDirPLate) & (distGivenDirPLate < P.sigAlpha);
intSig = isfinite(intPLate) & (intPLate < P.sigAlpha);
usable = isfinite(ctrlPLate);

fprintf('Late direction sites in base analysis: %d / %d\n', nnz(origSig), nIT);
fprintf('Late direction-significant sites (direction-only model): %d / %d\n', nnz(dirOnlySig), nIT);
fprintf('Late direction sites after distance control: %d / %d\n', nnz(ctrlSig), nIT);
fprintf('Late distance-significant sites (distance-only model): %d / %d\n', nnz(distSig), nIT);
fprintf('Late distance sites after direction control: %d / %d\n', nnz(distGivenDirSig), nIT);
fprintf('Late direction-distance interaction sites: %d / %d\n', nnz(intSig), nIT);
fprintf('Direction survives distance control in %d / %d originally significant late sites\n', ...
    nnz(origSig & ctrlSig), nnz(origSig));
if any(origSig & usable)
    fprintf('Median direction partial R^2 beyond distance, among original late significant sites: %.3f\n', ...
        median(partialR2Dir(origSig & usable), 'omitnan'));
end
if any(distGivenDirSig)
    fprintf('Median distance partial R^2 beyond direction, among distance|direction sites: %.3f\n', ...
        median(partialR2DistGivenDir(distGivenDirSig), 'omitnan'));
end
if any(intSig)
    fprintf('Median interaction partial R^2, among interaction-significant sites: %.3f\n', ...
        median(partialR2Interaction(intSig), 'omitnan'));
end
fprintf('Median additive unique R^2 | direction=%.3f distance=%.3f shared=%.3f\n', ...
    median(uniqueR2Dir(isfinite(uniqueR2Dir)), 'omitnan'), ...
    median(uniqueR2Dist(isfinite(uniqueR2Dist)), 'omitnan'), ...
    median(sharedR2Additive(isfinite(sharedR2Additive)), 'omitnan'));
fprintf('Median total model R^2 | additive=%.3f interaction=%.3f\n', ...
    median(totalR2Additive(isfinite(totalR2Additive)), 'omitnan'), ...
    median(totalR2Interaction(isfinite(totalR2Interaction)), 'omitnan'));

OUT.P = P;
OUT.monkeySuffix = monkeySuffix;
OUT.RFrange = RFrange;
OUT.QuartetTable = QuartetTable;
OUT.distTargetPx = distTargetPx;
OUT.distDistractorPx = distDistractorPx;
OUT.targetAdvPx = targetAdvPx;
OUT.meanDistPx = meanDistPx;
OUT.RegressionLate = RegressionLate;
OUT.origPLate = origPLate;
OUT.dirOnlyPLate = dirOnlyPLate;
OUT.ctrlPLate = ctrlPLate;
OUT.distPLate = distPLate;
OUT.distGivenDirPLate = distGivenDirPLate;
OUT.intPLate = intPLate;
OUT.partialR2Dir = partialR2Dir;
OUT.partialR2DistGivenDir = partialR2DistGivenDir;
OUT.partialR2Interaction = partialR2Interaction;
OUT.uniqueR2Dir = uniqueR2Dir;
OUT.uniqueR2Dist = uniqueR2Dist;
OUT.sharedR2Additive = sharedR2Additive;
OUT.totalR2Additive = totalR2Additive;
OUT.totalR2Interaction = totalR2Interaction;
OUT.origSig = origSig;
OUT.dirOnlySig = dirOnlySig;
OUT.ctrlSig = ctrlSig;
OUT.distSig = distSig;
OUT.distGivenDirSig = distGivenDirSig;
OUT.intSig = intSig;
OUT.siteSessionExclusions = currentExclusions;

if P.saveResult
    save(outPath, 'OUT', '-v7.3');
    fprintf('Saved IT distance-control analysis to %s\n', outPath);
end

if P.makeSummaryFigure
    make_summary_figure(OUT);
end
end

function [distTargetPx, distDistractorPx] = compute_arm_distances_all_stim(ALLCOORDS, RTAB384, Tall_IT, nIT)
nStim = numel(Tall_IT);
distTargetPx = nan(nIT, nStim);
distDistractorPx = nan(nIT, nStim);

fprintf('Computing RF-to-arm distances for %d IT sites across %d stimuli...\n', nIT, nStim);
for stimNum = 1:nStim
    [~, masks] = render_stim_with_masks2(ALLCOORDS, RTAB384, stimNum, 'DrawDots', false);
    dTargetMap = bwdist(masks.figArm);
    dDistrMap = bwdist(masks.backArm);

    T = Tall_IT(stimNum).T;
    assert(height(T) >= nIT, 'Tall_IT(%d).T has fewer than %d IT rows.', stimNum, nIT);
    x = round(double(T.x_px(1:nIT)));
    y = round(double(T.y_px(1:nIT)));
    H = size(dTargetMap, 1);
    W = size(dTargetMap, 2);
    inB = x >= 1 & x <= W & y >= 1 & y <= H;
    if ismember('inBounds', string(T.Properties.VariableNames))
        inB = inB & logical(T.inBounds(1:nIT));
    end

    x = min(max(x, 1), W);
    y = min(max(y, 1), H);
    idx = sub2ind([H W], y, x);

    dT = dTargetMap(idx);
    dD = dDistrMap(idx);
    dT(~inB) = NaN;
    dD(~inB) = NaN;

    distTargetPx(:, stimNum) = dT;
    distDistractorPx(:, stimNum) = dD;

    if mod(stimNum, 50) == 0 || stimNum == 1 || stimNum == nStim
        fprintf('  Stimulus %d / %d\n', stimNum, nStim);
    end
end
end

function m = rowmean_omitnan(X)
good = isfinite(X);
count = sum(good, 2);
X(~good) = 0;
m = sum(X, 2) ./ max(count, 1);
m(count == 0) = NaN;
end

function R = init_regression_struct(nSites)
tmp = struct( ...
    'status', "not_fit", ...
    'nObs', NaN, ...
    'beta0', NaN, ...
    'betaDir0', NaN, ...
    'betaDirCos', NaN, ...
    'betaDirSin', NaN, ...
    'betaCos', NaN, ...
    'betaSin', NaN, ...
    'betaTargetAdv', NaN, ...
    'betaMeanDist', NaN, ...
    'betaAdvCos', NaN, ...
    'betaAdvSin', NaN, ...
    'rssConst', NaN, ...
    'rssDir', NaN, ...
    'rssDist', NaN, ...
    'rssFull', NaN, ...
    'rssInt', NaN, ...
    'pDirOnly', NaN, ...
    'pDistOnly', NaN, ...
    'pDirGivenDist', NaN, ...
    'pDistGivenDir', NaN, ...
    'pInteraction', NaN, ...
    'r2DirModel', NaN, ...
    'r2DistModel', NaN, ...
    'r2Additive', NaN, ...
    'uniqueR2Dir', NaN, ...
    'uniqueR2Dist', NaN, ...
    'sharedR2Additive', NaN, ...
    'partialR2Dist', NaN, ...
    'partialR2Dir', NaN, ...
    'partialR2DistGivenDir', NaN, ...
    'partialR2Interaction', NaN, ...
    'dfErrorFull', NaN);
R = repmat(tmp, nSites, 1);
end

function R = fit_direction_distance_control(thetaDeg, y, varY, targetAdv, meanDist, P)
R = init_regression_struct(1);

valid = isfinite(thetaDeg) & isfinite(y) & isfinite(varY) & (varY > 0) & ...
    isfinite(targetAdv) & isfinite(meanDist);
theta = thetaDeg(valid);
yv = y(valid);
vv = varY(valid);
adv = targetAdv(valid);
md = meanDist(valid);
R.nObs = numel(yv);

if numel(yv) < P.minQuartets
    R.status = "too_few_quartets";
    return;
end

varFloor = max(median(vv(vv > 0), 'omitnan') * P.varFloorFrac, 1e-6);
vv = max(vv, varFloor);
sqrtW = sqrt(1 ./ vv);

advZ = zscore_safe(adv);
mdZ = zscore_safe(md);

Xconst = ones(numel(theta), 1);
Xdir = [ones(numel(theta),1), cosd(theta), sind(theta)];
Xdist = [ones(numel(theta),1), advZ, mdZ];
Xfull = [ones(numel(theta),1), cosd(theta), sind(theta), advZ, mdZ];
Xint = [ones(numel(theta),1), cosd(theta), sind(theta), advZ, mdZ, ...
    advZ .* cosd(theta), advZ .* sind(theta)];

betaConst = weighted_lsq(Xconst, yv, sqrtW);
betaDir = weighted_lsq(Xdir, yv, sqrtW);
betaDist = weighted_lsq(Xdist, yv, sqrtW);
betaFull = weighted_lsq(Xfull, yv, sqrtW);
betaInt = weighted_lsq(Xint, yv, sqrtW);

R.betaDir0 = betaDir(1);
R.betaDirCos = betaDir(2);
R.betaDirSin = betaDir(3);
R.beta0 = betaFull(1);
R.betaCos = betaFull(2);
R.betaSin = betaFull(3);
R.betaTargetAdv = betaFull(4);
R.betaMeanDist = betaFull(5);
R.betaAdvCos = betaInt(6);
R.betaAdvSin = betaInt(7);

rssConst = weighted_rss(Xconst, yv, betaConst, vv);
rssDir = weighted_rss(Xdir, yv, betaDir, vv);
rssDist = weighted_rss(Xdist, yv, betaDist, vv);
rssFull = weighted_rss(Xfull, yv, betaFull, vv);
rssInt = weighted_rss(Xint, yv, betaInt, vv);

R.rssConst = rssConst;
R.rssDir = rssDir;
R.rssDist = rssDist;
R.rssFull = rssFull;
R.rssInt = rssInt;
R.dfErrorFull = numel(yv) - size(Xfull, 2);

[R.pDirOnly, ~] = extra_model_test(rssConst, rssDir, 2, numel(yv) - size(Xdir, 2));
[R.pDistOnly, R.partialR2Dist] = extra_model_test(rssConst, rssDist, 2, numel(yv) - size(Xdist, 2));
[R.pDirGivenDist, R.partialR2Dir] = extra_model_test(rssDist, rssFull, 2, R.dfErrorFull);
[R.pDistGivenDir, R.partialR2DistGivenDir] = extra_model_test(rssDir, rssFull, 2, R.dfErrorFull);
[R.pInteraction, R.partialR2Interaction] = extra_model_test(rssFull, rssInt, 2, numel(yv) - size(Xint, 2));
R.r2DirModel = model_r2(rssConst, rssDir);
R.r2DistModel = model_r2(rssConst, rssDist);
R.r2Additive = model_r2(rssConst, rssFull);
R.uniqueR2Dir = max((rssDist - rssFull) / max(rssConst, eps), 0);
R.uniqueR2Dist = max((rssDir - rssFull) / max(rssConst, eps), 0);
R.sharedR2Additive = R.r2Additive - R.uniqueR2Dir - R.uniqueR2Dist;
R.status = "ok";
end

function beta = weighted_lsq(X, y, sqrtW)
Xw = bsxfun(@times, X, sqrtW);
yw = sqrtW .* y;
beta = Xw \ yw;
end

function rss = weighted_rss(X, y, beta, varY)
yHat = X * beta;
rss = sum((1 ./ varY) .* (y - yHat).^2);
end

function r2 = model_r2(rssNull, rssModel)
r2 = (rssNull - rssModel) / max(rssNull, eps);
end

function [pVal, partialR2] = extra_model_test(rssRed, rssFull, dfExtra, dfError)
pVal = NaN;
partialR2 = NaN;
if ~(isfinite(rssRed) && isfinite(rssFull) && isfinite(dfError) && dfError > 0)
    return;
end

num = max(rssRed - rssFull, 0) / dfExtra;
den = max(rssFull, eps) / dfError;
F = num / max(den, eps);
pVal = 1 - fcdf_local(F, dfExtra, dfError);
partialR2 = max(min((rssRed - rssFull) / max(rssRed, eps), 1), 0);
end

function z = zscore_safe(x)
x = double(x(:));
mu = mean(x, 'omitnan');
sd = std(x, 0, 'omitnan');
if ~(isfinite(sd) && sd > 0)
    z = zeros(size(x));
else
    z = (x - mu) / sd;
end
end

function scatter_pvalues(pOrig, pCtrl, sigAlpha, ttl)
epsP = 1e-6;
x = max(min(double(pOrig(:)), 1), epsP);
y = max(min(double(pCtrl(:)), 1), epsP);
ok = isfinite(x) & isfinite(y);
loglog(x(ok), y(ok), 'o', 'MarkerSize', 5, ...
    'MarkerFaceColor', [0.25 0.55 0.80], 'MarkerEdgeColor', 'none');
hold on;
plot([epsP 1], [epsP 1], 'k--', 'LineWidth', 1);
plot([sigAlpha sigAlpha], [epsP 1], ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
plot([epsP 1], [sigAlpha sigAlpha], ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
xlim([epsP 1]);
ylim([epsP 1]);
xlabel('Original late direction p');
ylabel('Late direction p | distance');
title(ttl);
grid on;
axis square;
end

function scatter_unique_r2(uniqueDir, uniqueDist, sigDir, sigDist, ttl)
x = double(uniqueDir(:));
y = double(uniqueDist(:));
ok = isfinite(x) & isfinite(y);
x = x(ok);
y = y(ok);
isBoth = false(size(x));
if nargin >= 4 && ~isempty(sigDir) && ~isempty(sigDist)
    sDir = logical(sigDir(:));
    sDist = logical(sigDist(:));
    isBoth = sDir(ok) & sDist(ok);
end

plot(x(~isBoth), y(~isBoth), 'o', 'MarkerSize', 5, ...
    'MarkerFaceColor', [0.70 0.70 0.70], 'MarkerEdgeColor', 'none');
hold on;
plot(x(isBoth), y(isBoth), 'o', 'MarkerSize', 6, ...
    'MarkerFaceColor', [0.25 0.55 0.80], 'MarkerEdgeColor', 'none');
mx = max([x(:); y(:); 0.01]);
plot([0 mx], [0 mx], 'k--', 'LineWidth', 1);
xlim([0 mx]);
ylim([0 mx]);
xlabel('Unique direction R^2');
ylabel('Unique distance R^2');
title(ttl);
grid on;
axis square;
end

function make_summary_figure(OUT)
figure('Color', 'w');
useTiled = exist('tiledlayout', 'file') == 2;
if useTiled
    tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');
end

if useTiled, nexttile; else, subplot(1,2,1); end
scatter_unique_r2(OUT.uniqueR2Dir, OUT.uniqueR2Dist, OUT.ctrlSig, OUT.distGivenDirSig, ...
    sprintf('Unique Additive Variance (%s)', char(OUT.monkeySuffix)));

if useTiled, nexttile; else, subplot(1,2,2); end
vals = [nnz(OUT.dirOnlySig), nnz(OUT.distSig), nnz(OUT.ctrlSig), nnz(OUT.distGivenDirSig), nnz(OUT.intSig)];
bar(vals, 'FaceColor', [0.25 0.55 0.80], 'EdgeColor', 'none');
set(gca, 'XTick', 1:5, 'XTickLabel', {'dir model','dist model','dir|dist','dist|dir','dir x dist'});
ylabel('N sites');
title(sprintf('Late regression counts (%s)', char(OUT.monkeySuffix)));
grid on;

txt = sprintf('Median total R^2: additive %.3f | interaction %.3f', ...
    median(OUT.totalR2Additive(isfinite(OUT.totalR2Additive)), 'omitnan'), ...
    median(OUT.totalR2Interaction(isfinite(OUT.totalR2Interaction)), 'omitnan'));
if exist('sgtitle', 'file') == 2
    sgtitle(txt);
else
    annotation(gcf, 'textbox', [0.24 0.93 0.55 0.05], 'String', txt, ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold');
end
end

function p = fcdf_local(x, df1, df2)
z = (df1 .* x) ./ (df1 .* x + df2);
z = min(max(z, 0), 1);
p = betainc(z, df1/2, df2/2);
end
