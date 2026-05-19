function OUT = Show_TargetDirection_Interaction_IT()
% SHOW_TARGETDIRECTION_INTERACTION_IT
% Illustrate strong IT direction-by-distance interaction sites.
%
% For each selected late IT site:
% - plot the quartet attention-direction effect versus direction
% - overlay model-predicted direction curves for low and high targetAdv
% - show how the predicted effect changes with targetAdv at the preferred
%   direction and at the opposite direction

%% Settings
P = struct();
P.Monkey = 1;              % 1 = Nilson, 2 = Figaro
P.nSites = 3;
P.requireDirGivenDist = true;
P.requireInteraction = true;
P.minQuartets = 20;
P.minPrefSeparationDeg = 50;
P.lowQuantilePct = 25;
P.highQuantilePct = 75;
P.prefWindowDeg = 45;
P.makeFigures = true;
P.sigAlpha = 0.05;
P.varFloorFrac = 1e-3;

cfg = config();

if P.Monkey == 1
    monkeySuffix = "N";
elseif P.Monkey == 2
    monkeySuffix = "F";
else
    error('Show_TargetDirection_Interaction_IT:InvalidMonkey', ...
        'P.Monkey must be 1 (Nilson) or 2 (Figaro).');
end

basePath = fullfile(cfg.matDir, sprintf('Attention_TargetSide_Tuning_IT_directiondelta_%s.mat', char(monkeySuffix)));
distPath = fullfile(cfg.matDir, sprintf('Attention_TargetSide_DistanceControl_IT_%s.mat', char(monkeySuffix)));
hasSessionExclusions = ~isempty(site_session_exclusions(monkeySuffix));

if hasSessionExclusions
    fprintf(['Session exclusions are active for monkey %s; refreshing IT target-side ' ...
             'fits before plotting interaction examples.\n'], char(monkeySuffix));
    BASE = Attention_TargetSide_Tuning_IT(struct( ...
        'makeSummaryFigures', false, ...
        'makeExampleFigures', false, ...
        'makeAngleReferenceFigure', false));
    DIST = Attention_TargetSide_DistanceControl_IT(struct('makeSummaryFigure', false));
else
    assert(exist(basePath, 'file') == 2, ...
        'Missing %s. Run Attention_TargetSide_Tuning_IT first.', basePath);
    assert(exist(distPath, 'file') == 2, ...
        'Missing %s. Run Attention_TargetSide_DistanceControl_IT first.', distPath);
    Sbase = load(basePath, 'OUT');
    Sdist = load(distPath, 'OUT');
    BASE = Sbase.OUT;
    DIST = Sdist.OUT;
end

requiredBase = {'FitLate','deltaQuartetLate','varQuartetLate','QuartetTable','RFrange'};
for k = 1:numel(requiredBase)
    assert(isfield(BASE, requiredBase{k}), 'Base OUT missing field %s.', requiredBase{k});
end
requiredDist = {'targetAdvPx','meanDistPx','partialR2Interaction','partialR2Dir', ...
    'partialR2DistGivenDir','intSig','ctrlSig','RegressionLate'};
for k = 1:numel(requiredDist)
    assert(isfield(DIST, requiredDist{k}), 'Distance OUT missing field %s.', requiredDist{k});
end

FitLate = BASE.FitLate;
thetaDeg = double(BASE.QuartetTable.targetDirDeg(:));
deltaLate = double(BASE.deltaQuartetLate);
varLate = double(BASE.varQuartetLate);
targetAdvPx = double(DIST.targetAdvPx);
meanDistPx = double(DIST.meanDistPx);
prefDegAll = mod([FitLate.prefDeg].', 360);
RFrange = BASE.RFrange(:);
nSitesAll = numel(RFrange);

pool = true(nSitesAll, 1);
if P.requireInteraction
    pool = pool & logical(DIST.intSig(:));
end
if P.requireDirGivenDist
    pool = pool & logical(DIST.ctrlSig(:));
end
pool = pool & isfinite(DIST.partialR2Interaction(:)) & isfinite(prefDegAll);
sitePool = find(pool);
assert(~isempty(sitePool), 'No IT sites pass the requested interaction selection.');

siteSel = select_diverse_interaction_sites(sitePool, prefDegAll, DIST.partialR2Interaction(:), P.nSites, P.minPrefSeparationDeg);
assert(~isempty(siteSel), 'Could not select any interaction example sites.');

OUT = struct();
OUT.P = P;
OUT.monkeySuffix = monkeySuffix;
OUT.siteSel = siteSel;
OUT.globalSiteSel = RFrange(siteSel);
OUT.prefDegSel = prefDegAll(siteSel);

for ii = 1:numel(siteSel)
    s = siteSel(ii);
    theta = thetaDeg(:);
    y = deltaLate(s, :).';
    varY = varLate(s, :).';
    adv = targetAdvPx(s, :).';
    md = meanDistPx(s, :).';

    F = fit_interaction_model_local(theta, y, varY, adv, md, P);
    assert(F.nObs >= 20, 'Selected site %d has too few valid quartets.', s);

    qLo = prctile(F.advRaw, P.lowQuantilePct);
    qHi = prctile(F.advRaw, P.highQuantilePct);
    isLow = F.advRaw <= qLo;
    isHigh = F.advRaw >= qHi;

    prefDeg = mod(FitLate(s).prefDeg, 360);
    oppDeg = mod(prefDeg + 180, 360);
    thetaGrid = linspace(0, 360, 361).';
    yLow = predict_interaction_curve(F, thetaGrid, qLo, median(F.mdRaw, 'omitnan'));
    yHigh = predict_interaction_curve(F, thetaGrid, qHi, median(F.mdRaw, 'omitnan'));

    advGrid = linspace(min(F.advRaw), max(F.advRaw), 200).';
    yPref = predict_interaction_curve(F, prefDeg * ones(size(advGrid)), advGrid, median(F.mdRaw, 'omitnan'));
    yOpp = predict_interaction_curve(F, oppDeg * ones(size(advGrid)), advGrid, median(F.mdRaw, 'omitnan'));

    nearPref = abs(local_angdiff_deg(F.theta, prefDeg)) <= P.prefWindowDeg;
    nearOpp = abs(local_angdiff_deg(F.theta, oppDeg)) <= P.prefWindowDeg;

    Site(ii).localSite = s; %#ok<AGROW>
    Site(ii).globalSite = RFrange(s); %#ok<AGROW>
    Site(ii).prefDeg = prefDeg; %#ok<AGROW>
    Site(ii).veLate = FitLate(s).r2TrainPct; %#ok<AGROW>
    Site(ii).partialR2Interaction = DIST.partialR2Interaction(s); %#ok<AGROW>
    Site(ii).partialR2Dir = DIST.partialR2Dir(s); %#ok<AGROW>
    Site(ii).partialR2DistGivenDir = DIST.partialR2DistGivenDir(s); %#ok<AGROW>
    Site(ii).pInteraction = DIST.intPLate(s); %#ok<AGROW>
    Site(ii).theta = F.theta; %#ok<AGROW>
    Site(ii).delta = F.y; %#ok<AGROW>
    Site(ii).targetAdv = F.advRaw; %#ok<AGROW>
    Site(ii).meanDist = F.mdRaw; %#ok<AGROW>
    Site(ii).isLow = isLow; %#ok<AGROW>
    Site(ii).isHigh = isHigh; %#ok<AGROW>
    Site(ii).nearPref = nearPref; %#ok<AGROW>
    Site(ii).nearOpp = nearOpp; %#ok<AGROW>
    Site(ii).lowQuantile = qLo; %#ok<AGROW>
    Site(ii).highQuantile = qHi; %#ok<AGROW>
    Site(ii).thetaGrid = thetaGrid; %#ok<AGROW>
    Site(ii).yLow = yLow; %#ok<AGROW>
    Site(ii).yHigh = yHigh; %#ok<AGROW>
    Site(ii).advGrid = advGrid; %#ok<AGROW>
    Site(ii).yPref = yPref; %#ok<AGROW>
    Site(ii).yOpp = yOpp; %#ok<AGROW>

    if P.makeFigures
        make_interaction_site_figure(Site(ii), char(monkeySuffix), P);
    end
end

OUT.Site = Site;
end

function siteSel = select_diverse_interaction_sites(sitePool, prefDeg, score, nPick, minSepDeg)
sitePool = sitePool(:);
[~, ord] = sort(score(sitePool), 'descend');
siteOrd = sitePool(ord);
siteSel = zeros(0,1);

for i = 1:numel(siteOrd)
    s = siteOrd(i);
    if isempty(siteSel)
        siteSel(end+1,1) = s; %#ok<AGROW>
    else
        d = abs(local_angdiff_deg(prefDeg(s), prefDeg(siteSel)));
        if all(d >= minSepDeg)
            siteSel(end+1,1) = s; %#ok<AGROW>
        end
    end
    if numel(siteSel) >= nPick
        break;
    end
end

if numel(siteSel) < nPick
    remain = setdiff(siteOrd, siteSel, 'stable');
    siteSel = [siteSel; remain(1:min(nPick - numel(siteSel), numel(remain)))]; %#ok<AGROW>
end
end

function F = fit_interaction_model_local(thetaDeg, y, varY, targetAdv, meanDist, P)
valid = isfinite(thetaDeg) & isfinite(y) & isfinite(varY) & (varY > 0) & ...
    isfinite(targetAdv) & isfinite(meanDist);

F = struct('nObs', 0, 'theta', [], 'y', [], 'varY', [], 'advRaw', [], 'mdRaw', [], ...
    'advMu', NaN, 'advSd', NaN, 'mdMu', NaN, 'mdSd', NaN, 'betaInt', nan(7,1));

theta = double(thetaDeg(valid));
yv = double(y(valid));
vv = double(varY(valid));
adv = double(targetAdv(valid));
md = double(meanDist(valid));
F.nObs = numel(yv);
F.theta = theta;
F.y = yv;
F.varY = vv;
F.advRaw = adv;
F.mdRaw = md;
if F.nObs < P.minQuartets
    return;
end

varFloor = max(median(vv(vv > 0), 'omitnan') * P.varFloorFrac, 1e-6);
vv = max(vv, varFloor);
sqrtW = sqrt(1 ./ vv);

[advZ, advMu, advSd] = zscore_safe_local(adv);
[mdZ, mdMu, mdSd] = zscore_safe_local(md);
Xint = [ones(numel(theta),1), cosd(theta), sind(theta), advZ, mdZ, ...
    advZ .* cosd(theta), advZ .* sind(theta)];
betaInt = weighted_lsq_local(Xint, yv, sqrtW);

F.betaInt = betaInt;
F.advMu = advMu;
F.advSd = advSd;
F.mdMu = mdMu;
F.mdSd = mdSd;
end

function yHat = predict_interaction_curve(F, thetaDeg, advRaw, mdRaw)
advZ = normalize_with_stats(double(advRaw(:)), F.advMu, F.advSd);
mdZ = normalize_with_stats(double(mdRaw(:)), F.mdMu, F.mdSd);
thetaDeg = double(thetaDeg(:));
b = F.betaInt;
yHat = b(1) + b(2) * cosd(thetaDeg) + b(3) * sind(thetaDeg) + ...
    b(4) .* advZ + b(5) .* mdZ + ...
    b(6) .* advZ .* cosd(thetaDeg) + b(7) .* advZ .* sind(thetaDeg);
end

function make_interaction_site_figure(S, monkeySuffix, P)
fig = figure('Color', 'w', ...
    'Name', sprintf('IT target-direction interaction site %d', S.globalSite), ...
    'NumberTitle', 'off', 'Position', [80 140 1180 420]);
useTiled = exist('tiledlayout', 'file') == 2;
if useTiled
    tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');
end

if useTiled
    axTmp = nexttile;
else
    axTmp = subplot(1,2,1);
end
axPos = get(axTmp, 'Position');
delete(axTmp);
pax = polaraxes('Parent', fig, 'Position', axPos);
hold(pax, 'on');

allVals = [S.delta(:); S.yLow(:); S.yHigh(:)];
deltaMax = max(abs(allVals(isfinite(allVals))));
if ~(isfinite(deltaMax) && deltaMax > 0)
    deltaMax = 1;
end
r0 = 1.15 * deltaMax;
rAll = r0 + S.delta(:)';
rLow = r0 + S.yLow(:)';
rHigh = r0 + S.yHigh(:)';
thetaAll = deg2rad(S.theta(:)');
thetaGrid = deg2rad(S.thetaGrid(:)');

hAll = polarplot(pax, thetaAll, rAll, 'o', 'MarkerSize', 4.5, ...
    'MarkerFaceColor', [0.78 0.78 0.78], 'MarkerEdgeColor', 'none');
hLowPts = polarplot(pax, thetaAll(S.isLow(:)'), rAll(S.isLow(:)'), 'o', 'MarkerSize', 5.5, ...
    'MarkerFaceColor', [0.18 0.45 0.85], 'MarkerEdgeColor', 'none');
hHighPts = polarplot(pax, thetaAll(S.isHigh(:)'), rAll(S.isHigh(:)'), 'o', 'MarkerSize', 5.5, ...
    'MarkerFaceColor', [0.90 0.45 0.15], 'MarkerEdgeColor', 'none');
hLowFit = polarplot(pax, thetaGrid, rLow, '-', 'Color', [0.18 0.45 0.85], 'LineWidth', 2.2);
hHighFit = polarplot(pax, thetaGrid, rHigh, '-', 'Color', [0.90 0.45 0.15], 'LineWidth', 2.2);
polarplot(pax, thetaGrid, r0 * ones(size(thetaGrid)), '--', ...
    'Color', [0.70 0.70 0.70], 'LineWidth', 1);
polarplot(pax, deg2rad([S.prefDeg S.prefDeg]), [0 r0 + 1.10 * deltaMax], '--', ...
    'Color', [0.35 0.35 0.35], 'LineWidth', 1);
polarplot(pax, deg2rad([mod(S.prefDeg + 180, 360) mod(S.prefDeg + 180, 360)]), [0 r0 + 1.10 * deltaMax], ':', ...
    'Color', [0.45 0.45 0.45], 'LineWidth', 1);

pax.ThetaZeroLocation = 'right';
pax.ThetaDir = 'counterclockwise';
pax.ThetaTick = [0 90 180 270];
pax.RLim = [0 max(r0 + 1.10 * deltaMax, max([rAll(:); rLow(:); rHigh(:)]))];
tickVals = linspace(-deltaMax, deltaMax, 5);
pax.RTick = r0 + tickVals;
pax.RTickLabel = arrayfun(@(v) sprintf('%.2g', v), tickVals, 'UniformOutput', false);
title(pax, sprintf('Site %d | low RF dist (Tar-Dis) q%.0f=%.1f px, high q%.0f=%.1f px', ...
    S.globalSite, P.lowQuantilePct, S.lowQuantile, P.highQuantilePct, S.highQuantile));
legend(pax, [hAll hLowPts hHighPts hLowFit hHighFit], ...
    {'all quartets','low RF dist (Tar-Dis)','high RF dist (Tar-Dis)', ...
    'low RF dist fit','high RF dist fit'}, ...
    'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off');
grid(pax, 'on');

if useTiled, nexttile; else, subplot(1,2,2); end
hold on;
plot(S.targetAdv, S.delta, 'o', 'MarkerSize', 4.5, ...
    'MarkerFaceColor', [0.82 0.82 0.82], 'MarkerEdgeColor', 'none');
plot(S.targetAdv(S.nearPref), S.delta(S.nearPref), 'o', 'MarkerSize', 5.5, ...
    'MarkerFaceColor', [0.85 0.10 0.10], 'MarkerEdgeColor', 'none');
plot(S.targetAdv(S.nearOpp), S.delta(S.nearOpp), 'o', 'MarkerSize', 5.5, ...
    'MarkerFaceColor', [0.15 0.45 0.85], 'MarkerEdgeColor', 'none');
plot(S.advGrid, S.yPref, '-', 'Color', [0.85 0.10 0.10], 'LineWidth', 2.2);
plot(S.advGrid, S.yOpp, '-', 'Color', [0.15 0.45 0.85], 'LineWidth', 2.2);
plot(xlim, [0 0], '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
xlabel('RF dist (Tar-Dis) (px)');
ylabel('\Delta_q (spont. SD units)');
title(sprintf('pref %.0f%c vs opposite %.0f%c', S.prefDeg, char(176), mod(S.prefDeg + 180, 360), char(176)));
legend({ ...
    'all quartets', ...
    sprintf('|dir-pref| \\leq %d%c', P.prefWindowDeg, char(176)), ...
    sprintf('|dir-opp| \\leq %d%c', P.prefWindowDeg, char(176)), ...
    'predicted at pref', ...
    'predicted at opposite'}, ...
    'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off');
grid on;

annotation(fig, 'textbox', [0.16 0.91 0.68 0.07], ...
    'String', sprintf(['IT site %d (%s) | pref %.0f%c | late VE %.1f%% | ' ...
    'unique dir R^2 %.3f | unique dist R^2 %.3f | interaction R^2 %.3f | p_{int}=%.2g'], ...
    S.globalSite, monkeySuffix, S.prefDeg, char(176), ...
    S.veLate, S.partialR2Dir, S.partialR2DistGivenDir, S.partialR2Interaction, S.pInteraction), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'FontSize', 14, 'FontWeight', 'bold');
end

function beta = weighted_lsq_local(X, y, sqrtW)
Xw = bsxfun(@times, X, sqrtW);
yw = sqrtW .* y;
beta = Xw \ yw;
end

function [z, mu, sd] = zscore_safe_local(x)
x = double(x(:));
mu = mean(x, 'omitnan');
sd = std(x, 0, 'omitnan');
if ~(isfinite(sd) && sd > 0)
    z = zeros(size(x));
    sd = NaN;
else
    z = (x - mu) / sd;
end
end

function z = normalize_with_stats(x, mu, sd)
if ~(isfinite(sd) && sd > 0)
    z = zeros(size(x));
else
    z = (x - mu) / sd;
end
end

function d = local_angdiff_deg(a, b)
d = mod((double(a) - double(b)) + 180, 360) - 180;
end
