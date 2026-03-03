function G = build_grouped_alonggc_polygons_allbins(OUT_postAffine, Tall_V1, varargin)
% BUILD_GROUPED_ALONGGC_POLYGONS_ALLBINS
% Group significant (site,quartet) combinations by along_GC (equal counts),
% compute group-mean signed time courses, and build per-group polygons for
% target/distractor streams on post-affine coordinates.
%
% Required:
%   OUT_postAffine : struct with OUT_postAffine.bins from
%                    compute_projected_delta_points_allbins
%   Tall_V1        : stimulus RF table struct array
%
% Name/value:
%   'nGroups'         (default 8)
%   'sigSiteByIndex'  (default []) logical vector over absolute site index
%   'siteWeightsByIndex' (default []) numeric vector over absolute site index
%   'preEndMs'        (default 0)
%   'postStartMs'     (default 300)
%   'preQuantilePct'  (default 95)
%   'cMaxPostPct'     (default 95)
%   'polygonShrink'   (default 0.8)
%   'saveFile'        (default '')
%   'verbose'         (default true)
%
% Output G:
%   G.comboTable          one row per valid grouped combo
%   G.groupSummary        per-group counts/ranges
%   G.groupPolygons(g,s)  polygon struct for group g, stream s(1=T,2=D)
%   G.groupMeanSigned     [nGroups x nBins]
%   G.groupNPerBin        [nGroups x nBins]
%   G.timeWindows         [nBins x 2]
%   G.preMask, G.postMask
%   G.calibration         threshold/color suggestions from grouped means

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
B = OUT_postAffine.bins;
nBins = numel(B);

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
stimSrcSQS = accumarray(gSQS, stimAll, [], @first_finite);

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

% along_GC per combo.
along = nan(nComb,1);
for i = 1:nComb
    stimIdx = round(stimSrc(i));
    siteIdx = round(uSQ(i,1));
    if ~isfinite(stimIdx) || ~isfinite(siteIdx) || ...
            stimIdx < 1 || stimIdx > numel(Tall_V1) || siteIdx < 1
        continue;
    end
    if ~isfield(Tall_V1(stimIdx), 'T') || ~istable(Tall_V1(stimIdx).T)
        continue;
    end
    T = Tall_V1(stimIdx).T;
    if siteIdx > height(T) || ~ismember('along_GC', T.Properties.VariableNames)
        continue;
    end
    along(i) = double(T.along_GC(siteIdx));
end

% Optional site-level significance filter.
isSigSite = true(nComb,1);
if ~isempty(opt.sigSiteByIndex)
    sigSiteMask = opt.sigSiteByIndex(:);
    isSigSite = false(nComb,1);
    siteIdx = round(uSQ(:,1));
    okSite = siteIdx >= 1 & siteIdx <= numel(sigSiteMask);
    isSigSite(okSite) = sigSiteMask(siteIdx(okSite));
end

% Optional per-site weights (same weight for all quartets of a site).
useSiteWeights = ~isempty(opt.siteWeightsByIndex);
comboWeight = ones(nComb,1);
if useSiteWeights
    siteW = double(opt.siteWeightsByIndex(:));
    comboWeight = zeros(nComb,1);
    siteIdx = round(uSQ(:,1));
    okSite = isfinite(siteIdx) & siteIdx >= 1 & siteIdx <= numel(siteW);
    comboWeight(okSite) = siteW(siteIdx(okSite));
    comboWeight(~isfinite(comboWeight) | comboWeight < 0) = 0;
end

isValid = isfinite(along) & isfinite(xT) & isfinite(yT) & isfinite(xD) & isfinite(yD) & isSigSite;
if useSiteWeights
    isValid = isValid & isfinite(comboWeight) & (comboWeight > 0);
end
uSQ = uSQ(isValid,:);
along = along(isValid);
xT = xT(isValid); yT = yT(isValid);
xD = xD(isValid); yD = yD(isValid);
stimSrc = round(stimSrc(isValid));
comboWeight = comboWeight(isValid);

nValid = size(uSQ,1);
nGroups = round(opt.nGroups);
assert(nValid >= nGroups, 'Only %d valid combos found; cannot split into %d groups.', nValid, nGroups);

% Equal-count grouping by sorted along_GC.
[~, ord] = sort(along, 'ascend');
groupIdx = zeros(nValid,1);
for r = 1:nValid
    g = floor((r-1) * nGroups / nValid) + 1;
    groupIdx(ord(r)) = g;
end

groupCounts = accumarray(groupIdx, 1, [nGroups 1], @sum, 0);
groupWeightSum = accumarray(groupIdx, comboWeight, [nGroups 1], @sum, 0);
groupWeightMean = accumarray(groupIdx, comboWeight, [nGroups 1], @mean, NaN);
alongMin = nan(nGroups,1);
alongMax = nan(nGroups,1);
for g = 1:nGroups
    v = along(groupIdx == g);
    if ~isempty(v)
        alongMin(g) = min(v);
        alongMax(g) = max(v);
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
        [groupMeanSigned(g,:), groupWeightPerBin(g,:)] = weighted_mean_omitnan(V, W);
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
    groupPolygons(g,1) = make_polygon([xT(idx), yT(idx)], opt.polygonShrink);
    groupPolygons(g,2) = make_polygon([xD(idx), yD(idx)], opt.polygonShrink);
end

comboTable = table( ...
    uint16(uSQ(:,1)), uint16(uSQ(:,2)), uint16(stimSrc), ...
    along, uint8(groupIdx), comboWeight, ...
    xT, yT, xD, yD, ...
    'VariableNames', {'siteIdx','quartetIdx','stimIdxSource','along_GC','groupIdx','siteWeight','xT','yT','xD','yD'});

groupSummary = table((1:nGroups)', groupCounts, alongMin, alongMax, groupWeightSum, groupWeightMean, ...
    'VariableNames', {'groupIdx','nComb','alongMin','alongMax','sumWeight','meanWeight'});

G = struct();
G.meta = struct( ...
    'created', datestr(now,30), ...
    'nGroups', nGroups, ...
    'nCombValid', nValid, ...
    'nBins', nBins, ...
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
    fprintf('Grouped along_GC prep: %d combos -> %d groups | bins=%d\n', nValid, nGroups, nBins);
    fprintf('  group sizes min/max: %d / %d\n', min(groupCounts), max(groupCounts));
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

function out = first_finite(v)
v = v(isfinite(v));
if isempty(v)
    out = NaN;
else
    out = v(1);
end
end

function P = make_polygon(pts, shrink)
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

function [mu, wsum] = weighted_mean_omitnan(V, w)
% Weighted mean per column, ignoring NaNs in V and invalid/nonpositive w.
[nRows, nCols] = size(V);
mu = nan(1, nCols);
wsum = zeros(1, nCols);
if isempty(V) || isempty(w)
    return;
end
w = double(w(:));
if numel(w) ~= nRows
    error('weighted_mean_omitnan:WeightSizeMismatch', ...
        'Weight vector length (%d) must match number of rows in V (%d).', numel(w), nRows);
end
w(~isfinite(w) | w <= 0) = 0;
for c = 1:nCols
    vc = double(V(:,c));
    ok = isfinite(vc) & (w > 0);
    if ~any(ok)
        continue;
    end
    ww = w(ok);
    vv = vc(ok);
    wsum(c) = sum(ww);
    if wsum(c) > 0
        mu(c) = sum(ww .* vv) / wsum(c);
    end
end
end
