function G = build_grouped_eccentricity_polygons_allbins(OUT_postAffine, rfEccBySiteIndex, varargin)
% BUILD_GROUPED_ECCENTRICITY_POLYGONS_ALLBINS
% Group valid (site,quartet) combinations by RF eccentricity (equal counts),
% compute group-mean signed time courses, and build per-group polygons for
% target/distractor streams on post-affine coordinates.
%
% Required:
%   OUT_postAffine    : struct with OUT_postAffine.bins from
%                       compute_projected_delta_points_allbins
%   rfEccBySiteIndex  : numeric vector indexed by siteIdx in OUT_postAffine
%
% Name/value:
%   'nGroups'         (default 8)
%   'sigSiteByIndex'  (default []) logical vector over site index
%   'siteWeightsByIndex' (default []) numeric vector over site index
%   'preEndMs'        (default 0)
%   'postStartMs'     (default 300)
%   'preQuantilePct'  (default 95)
%   'cMaxPostPct'     (default 95)
%   'polygonShrink'   (default 0.8)
%   'saveFile'        (default '')
%   'verbose'         (default true)
%
% Output G mirrors build_grouped_alonggc_polygons_allbins, but the grouping
% variable is RF eccentricity instead of along_GC.

p = inputParser;
p.addParameter('nGroups', 8, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('sigSiteByIndex', [], @(x) isempty(x) || (islogical(x) && isvector(x)));
p.addParameter('siteWeightsByIndex', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.addParameter('preEndMs', 0, @(x) isnumeric(x) && isscalar(x));
p.addParameter('postStartMs', 300, @(x) isnumeric(x) && isscalar(x));
p.addParameter('preQuantilePct', 95, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 100);
p.addParameter('cMaxPostPct', 95, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 100);
p.addParameter('polygonShrink', 0.8, @(x) isnumeric(x) && isscalar(x) && isfinite(x));
p.addParameter('saveFile', '', @(x) ischar(x) || isstring(x));
p.addParameter('verbose', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;

assert(isstruct(OUT_postAffine) && isfield(OUT_postAffine, 'bins') && ~isempty(OUT_postAffine.bins), ...
    'OUT_postAffine must contain non-empty bins.');
assert(isnumeric(rfEccBySiteIndex) && isvector(rfEccBySiteIndex) && ~isempty(rfEccBySiteIndex), ...
    'rfEccBySiteIndex must be a non-empty numeric vector.');

B = OUT_postAffine.bins;
nBins = numel(B);
rfEccBySiteIndex = double(rfEccBySiteIndex(:));

for tb = 1:nBins
    need = {'siteIdx','quartetIdx','stream','x_px','y_px','delta','stimIdxSource'};
    for k = 1:numel(need)
        assert(isfield(B(tb), need{k}), 'OUT_postAffine.bins(%d).%s missing.', tb, need{k});
    end
end

if isfield(OUT_postAffine, 'meta') && isfield(OUT_postAffine.meta, 'timeWindows')
    timeWindows = double(OUT_postAffine.meta.timeWindows);
else
    timeWindows = nan(nBins,2);
    for tb = 1:nBins
        timeWindows(tb,:) = double(B(tb).timeWindow);
    end
end

% Collect unique (site,quartet,stream) coordinates from all bins.
siteAll = [];
qAll = [];
streamAll = [];
xAll = [];
yAll = [];
stimAll = [];

for tb = 1:nBins
    s = double(B(tb).siteIdx(:));
    q = double(B(tb).quartetIdx(:));
    st = double(B(tb).stream(:));
    x = double(B(tb).x_px(:));
    y = double(B(tb).y_px(:));
    ss = double(B(tb).stimIdxSource(:));
    ok = isfinite(s) & isfinite(q) & isfinite(st) & isfinite(x) & isfinite(y) & isfinite(ss);
    if ~any(ok)
        continue;
    end

    siteAll = [siteAll; s(ok)]; %#ok<AGROW>
    qAll = [qAll; q(ok)]; %#ok<AGROW>
    streamAll = [streamAll; st(ok)]; %#ok<AGROW>
    xAll = [xAll; x(ok)]; %#ok<AGROW>
    yAll = [yAll; y(ok)]; %#ok<AGROW>
    stimAll = [stimAll; ss(ok)]; %#ok<AGROW>
end

assert(~isempty(siteAll), 'No valid post-affine points found in bins.');

keySQS = [siteAll qAll streamAll];
[uSQS, ~, gSQS] = unique(keySQS, 'rows');
xMean = accumarray(gSQS, xAll, [], @mean);
yMean = accumarray(gSQS, yAll, [], @mean);
stimSrcSQS = accumarray(gSQS, stimAll, [], @first_finite_local);

% Pivot to combo-level rows, with target and distractor coordinates.
keySQ = uSQS(:,1:2);
[uSQ, ~, gSQ] = unique(keySQ, 'rows');
nComb = size(uSQ,1);
xT = nan(nComb,1); yT = nan(nComb,1);
xD = nan(nComb,1); yD = nan(nComb,1);
stimSrc = nan(nComb,1);

for i = 1:size(uSQS,1)
    c = gSQ(i);
    st = uSQS(i,3);
    if st == 1
        xT(c) = xMean(i); yT(c) = yMean(i);
    elseif st == 2
        xD(c) = xMean(i); yD(c) = yMean(i);
    end
    if ~isfinite(stimSrc(c))
        stimSrc(c) = stimSrcSQS(i);
    end
end

siteIdxLocal = round(uSQ(:,1));
rfEcc = nan(nComb,1);
okSite = isfinite(siteIdxLocal) & siteIdxLocal >= 1 & siteIdxLocal <= numel(rfEccBySiteIndex);
rfEcc(okSite) = rfEccBySiteIndex(siteIdxLocal(okSite));

% Optional site-level significance filter.
isSigSite = true(nComb,1);
if ~isempty(opt.sigSiteByIndex)
    sigSiteMask = opt.sigSiteByIndex(:);
    isSigSite = false(nComb,1);
    okSig = isfinite(siteIdxLocal) & siteIdxLocal >= 1 & siteIdxLocal <= numel(sigSiteMask);
    isSigSite(okSig) = sigSiteMask(siteIdxLocal(okSig));
end

% Optional per-site weights (same weight for all quartets of a site).
useSiteWeights = ~isempty(opt.siteWeightsByIndex);
comboWeight = ones(nComb,1);
if useSiteWeights
    siteW = double(opt.siteWeightsByIndex(:));
    comboWeight = zeros(nComb,1);
    okW = isfinite(siteIdxLocal) & siteIdxLocal >= 1 & siteIdxLocal <= numel(siteW);
    comboWeight(okW) = siteW(siteIdxLocal(okW));
    comboWeight(~isfinite(comboWeight) | comboWeight < 0) = 0;
end

isValid = isfinite(rfEcc) & isfinite(xT) & isfinite(yT) & isfinite(xD) & isfinite(yD) & isSigSite;
if useSiteWeights
    isValid = isValid & isfinite(comboWeight) & (comboWeight > 0);
end
uSQ = uSQ(isValid,:);
siteIdxLocal = siteIdxLocal(isValid);
rfEcc = rfEcc(isValid);
xT = xT(isValid); yT = yT(isValid);
xD = xD(isValid); yD = yD(isValid);
stimSrc = round(stimSrc(isValid));
comboWeight = comboWeight(isValid);

nValid = size(uSQ,1);
nGroups = round(opt.nGroups);
assert(nValid >= nGroups, 'Only %d valid combos found; cannot split into %d groups.', nValid, nGroups);

% Equal-count grouping by sorted RF eccentricity.
[~, ord] = sort(rfEcc, 'ascend');
groupIdx = zeros(nValid,1);
for r = 1:nValid
    g = floor((r-1) * nGroups / nValid) + 1;
    groupIdx(ord(r)) = g;
end

groupCounts = accumarray(groupIdx, 1, [nGroups 1], @sum, 0);
groupWeightSum = accumarray(groupIdx, comboWeight, [nGroups 1], @sum, 0);
groupWeightMean = accumarray(groupIdx, comboWeight, [nGroups 1], @mean, NaN);
eccMin = nan(nGroups,1);
eccMax = nan(nGroups,1);
for g = 1:nGroups
    v = rfEcc(groupIdx == g);
    if ~isempty(v)
        eccMin(g) = min(v);
        eccMax(g) = max(v);
    end
end

% Combo delta matrix [nComb x nBins], one signed value per (site,quartet).
deltaMat = nan(nValid, nBins);
for tb = 1:nBins
    s = double(B(tb).siteIdx(:));
    q = double(B(tb).quartetIdx(:));
    d = double(B(tb).delta(:));
    ok = isfinite(s) & isfinite(q) & isfinite(d);
    if ~any(ok)
        continue;
    end

    keyTB = [s(ok) q(ok)];
    dTB = d(ok);
    [uTB, ia] = unique(keyTB, 'rows', 'stable');
    dTB = dTB(ia);

    [tf, loc] = ismember(uTB, uSQ, 'rows');
    if any(tf)
        deltaMat(loc(tf), tb) = dTB(tf);
    end
end

groupMeanSigned = nan(nGroups, nBins);
groupNPerBin = zeros(nGroups, nBins);
groupWeightPerBin = nan(nGroups, nBins);
for g = 1:nGroups
    idx = (groupIdx == g);
    if ~any(idx)
        continue;
    end
    V = deltaMat(idx,:);
    groupNPerBin(g,:) = sum(isfinite(V), 1);
    if useSiteWeights
        W = comboWeight(idx);
        [groupMeanSigned(g,:), groupWeightPerBin(g,:)] = weighted_mean_omitnan_local(V, W);
    else
        groupMeanSigned(g,:) = mean(V, 1, 'omitnan');
    end
end

preMask = timeWindows(:,2) <= opt.preEndMs;
postMask = timeWindows(:,1) >= opt.postStartMs;
assert(any(preMask), 'No pre bins satisfy end<=%.1f ms.', opt.preEndMs);
assert(any(postMask), 'No post bins satisfy start>=%.1f ms.', opt.postStartMs);

preVals = abs(groupMeanSigned(:, preMask));
preVals = preVals(isfinite(preVals));
postVals = abs(groupMeanSigned(:, postMask));
postVals = postVals(isfinite(postVals));

if isempty(preVals)
    thr = NaN;
    preEx = NaN;
else
    thr = prctile(preVals, opt.preQuantilePct);
    preEx = mean(preVals > thr);
end
if isempty(postVals) || ~isfinite(thr)
    postEx = NaN;
else
    postEx = mean(postVals > thr);
end
if isempty(postVals)
    cMax = NaN;
else
    cMax = prctile(postVals, opt.cMaxPostPct);
end

% Build polygons for each group and stream.
groupPolygons = repmat(struct('x', [], 'y', [], 'method', '', 'nPoints', 0), nGroups, 2);
for g = 1:nGroups
    idx = (groupIdx == g);
    groupPolygons(g,1) = make_polygon_local([xT(idx), yT(idx)], opt.polygonShrink);
    groupPolygons(g,2) = make_polygon_local([xD(idx), yD(idx)], opt.polygonShrink);
end

comboTable = table( ...
    uint16(siteIdxLocal), uint16(uSQ(:,2)), uint16(stimSrc), ...
    rfEcc, uint8(groupIdx), comboWeight, ...
    xT, yT, xD, yD, ...
    'VariableNames', {'siteIdx','quartetIdx','stimIdxSource','rfEcc','groupIdx','siteWeight','xT','yT','xD','yD'});

groupSummary = table((1:nGroups)', groupCounts, eccMin, eccMax, eccMin, eccMax, groupWeightSum, groupWeightMean, ...
    'VariableNames', {'groupIdx','nComb','eccMin','eccMax','alongMin','alongMax','sumWeight','meanWeight'});

G = struct();
G.meta = struct( ...
    'created', datestr(now,30), ...
    'nGroups', nGroups, ...
    'nCombValid', nValid, ...
    'nBins', nBins, ...
    'groupingVar', 'rfEccentricity', ...
    'preEndMs', opt.preEndMs, ...
    'postStartMs', opt.postStartMs, ...
    'preQuantilePct', opt.preQuantilePct, ...
    'cMaxPostPct', opt.cMaxPostPct, ...
    'polygonShrink', opt.polygonShrink, ...
    'useSiteWeights', useSiteWeights);
G.timeWindows = timeWindows;
G.preMask = preMask;
G.postMask = postMask;
G.comboTable = comboTable;
G.groupSummary = groupSummary;
G.groupPolygons = groupPolygons;
G.groupMeanSigned = groupMeanSigned;
G.groupNPerBin = groupNPerBin;
G.groupWeightPerBin = groupWeightPerBin;
G.calibration = struct( ...
    'thresholdPreQ', thr, ...
    'preExceedFrac', preEx, ...
    'postExceedFrac', postEx, ...
    'cMaxSuggest', cMax);

wp = comboWeight(comboWeight > 0 & isfinite(comboWeight));
if isempty(wp)
    wmin = NaN; wmed = NaN; wmax = NaN;
else
    wmin = min(wp);
    wmed = median(wp);
    wmax = max(wp);
end
G.weighting = struct( ...
    'enabled', useSiteWeights, ...
    'nCombWeighted', nnz(isfinite(comboWeight) & comboWeight > 0), ...
    'minWeight', wmin, ...
    'medianWeight', wmed, ...
    'maxWeight', wmax);

if opt.verbose
    fprintf('Grouped RF-ecc prep: %d combos -> %d groups | bins=%d\n', nValid, nGroups, nBins);
    fprintf('  group sizes min/max: %d / %d\n', min(groupCounts), max(groupCounts));
    fprintf('  eccentricity ranges: %.3f-%.3f ... %.3f-%.3f\n', ...
        eccMin(1), eccMax(1), eccMin(end), eccMax(end));
    if useSiteWeights
        fprintf('  site weights enabled: min/median/max = %.6g / %.6g / %.6g\n', wmin, wmed, wmax);
    else
        fprintf('  site weights disabled: using unweighted means.\n');
    end
    fprintf('  grouped threshold (pre p%.2f): %.6g | pre exceed %.2f%% | post exceed %.2f%%\n', ...
        opt.preQuantilePct, thr, 100*preEx, 100*postEx);
end

saveFile = char(opt.saveFile);
if ~isempty(saveFile)
    save(saveFile, 'G', '-v7.3');
end

end

function out = first_finite_local(v)
v = v(isfinite(v));
if isempty(v)
    out = NaN;
else
    out = v(1);
end
end

function P = make_polygon_local(pts, shrink)
P = struct('x', [], 'y', [], 'method', 'none', 'nPoints', 0);
if isempty(pts)
    return;
end
pts = pts(all(isfinite(pts),2), :);
pts = unique(pts, 'rows');
n = size(pts,1);
P.nPoints = n;
if n == 0
    return;
elseif n == 1
    P.x = pts(:,1);
    P.y = pts(:,2);
    P.method = 'point';
    return;
elseif n == 2
    P.x = [pts(1,1); pts(2,1); pts(1,1)];
    P.y = [pts(1,2); pts(2,2); pts(1,2)];
    P.method = 'line';
    return;
end

usedBoundary = false;
k = [];
if exist('boundary', 'file') == 2
    try
        k = boundary(pts(:,1), pts(:,2), shrink);
        usedBoundary = true;
    catch
        k = [];
    end
end
if isempty(k) || numel(unique(k)) < 3
    k = convhull(pts(:,1), pts(:,2));
    usedBoundary = false;
end
if k(1) ~= k(end)
    k = [k; k(1)];
end

P.x = pts(k,1);
P.y = pts(k,2);
if usedBoundary
    P.method = 'boundary';
else
    P.method = 'convhull';
end
end

function [mu, wsum] = weighted_mean_omitnan_local(V, w)
[nRows, nCols] = size(V);
mu = nan(1, nCols);
wsum = zeros(1, nCols);
if isempty(V) || isempty(w)
    return;
end
w = double(w(:));
if numel(w) ~= nRows
    error('weighted_mean_omitnan_local:WeightSizeMismatch', ...
        'Weight vector length (%d) must match number of rows in V (%d).', numel(w), nRows);
end

for c = 1:nCols
    x = V(:,c);
    ok = isfinite(x) & isfinite(w) & (w > 0);
    if ~any(ok)
        continue;
    end
    ww = w(ok);
    xx = x(ok);
    den = sum(ww);
    if ~(isfinite(den) && den > 0)
        continue;
    end
    mu(c) = sum(ww .* xx) / den;
    wsum(c) = den;
end
end
