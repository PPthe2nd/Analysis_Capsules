function [x, y, vSigned, stream] = smooth_post_affine_bin_knn(binData, K, enforceK, siteWeightsByIndex)
% SMOOTH_POST_AFFINE_BIN_KNN
% KNN smoothing of signed delta values within one post-affine bin.
% If binData.stream exists (1=target,2=distractor), smoothing is done
% within each stream to avoid mixing target/distractor signals.
% Optional siteWeightsByIndex applies fixed per-site weights to neighbor
% averaging (same weight reused across quartets/bins).

if nargin < 3 || isempty(enforceK)
    enforceK = true;
end
if nargin < 4
    siteWeightsByIndex = [];
end

xAll = double(binData.x_px(:));
yAll = double(binData.y_px(:));
dAll = double(binData.delta(:));
if isfield(binData, 'siteIdx')
    siteAll = double(binData.siteIdx(:));
else
    siteAll = nan(size(dAll));
end
if isfield(binData, 'stream')
    streamAll = double(binData.stream(:));
else
    streamAll = ones(size(dAll));
end
if isempty(siteWeightsByIndex)
    wAll = ones(size(dAll));
else
    wBySite = double(siteWeightsByIndex(:));
    wAll = zeros(size(dAll));
    siteRound = round(siteAll);
    okSite = isfinite(siteRound) & siteRound >= 1 & siteRound <= numel(wBySite);
    wAll(okSite) = wBySite(siteRound(okSite));
    wAll(~isfinite(wAll) | wAll < 0) = 0;
end

ok = isfinite(xAll) & isfinite(yAll) & isfinite(dAll) & isfinite(streamAll) & isfinite(wAll);

x = xAll(ok);
y = yAll(ok);
d = dAll(ok);
stream = streamAll(ok);
w = wAll(ok);
nP = numel(d);

if nP == 0
    vSigned = [];
    stream = [];
    return;
end

K = round(double(K));
if K < 1 || ~isfinite(K)
    error('K must be a positive finite integer.');
end

if nP < K
    if enforceK
        error('Frame has only %d points; cannot enforce K=%d.', nP, K);
    end
    K = nP;
end

vSigned = nan(nP,1);
uStream = unique(stream(:))';
for g = uStream
    idxG = find(stream == g);
    if isempty(idxG)
        continue;
    end
    Kg = min(K, numel(idxG));
    if Kg < K && enforceK
        error('Stream %d has only %d points; cannot enforce K=%d.', g, numel(idxG), K);
    end

    xg = x(idxG);
    yg = y(idxG);
    dg = d(idxG);
    wg = w(idxG);
    for ii = 1:numel(idxG)
        dist2 = (xg - xg(ii)).^2 + (yg - yg(ii)).^2;
        [~, ord] = sort(dist2, 'ascend');
        idx = ord(1:Kg);
        ww = wg(idx);
        if any(ww > 0)
            vSigned(idxG(ii)) = sum(ww .* dg(idx)) / sum(ww);
        else
            vSigned(idxG(ii)) = mean(dg(idx));
        end
    end
end

end
