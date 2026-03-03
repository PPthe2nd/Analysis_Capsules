function S = summarize_post_affine_threshold_exceedance(OUT_postAffine, K, threshold, preMask, postMask, enforceK, siteWeightsByIndex)
% SUMMARIZE_POST_AFFINE_THRESHOLD_EXCEEDANCE
% Compute threshold exceedance from post-affine bins using KNN-smoothed abs values.

if nargin < 6 || isempty(enforceK)
    enforceK = true;
end
if nargin < 7
    siteWeightsByIndex = [];
end

assert(isstruct(OUT_postAffine) && isfield(OUT_postAffine, 'bins') && ~isempty(OUT_postAffine.bins), ...
    'OUT_postAffine must contain non-empty bins.');

bins = OUT_postAffine.bins;
nBins = numel(bins);

assert(islogical(preMask) && isvector(preMask) && numel(preMask)==nBins, ...
    'preMask must be a logical vector with one entry per bin.');
assert(islogical(postMask) && isvector(postMask) && numel(postMask)==nBins, ...
    'postMask must be a logical vector with one entry per bin.');
assert(isfinite(threshold) && threshold > 0, 'threshold must be positive and finite.');

fracByBin = nan(nBins,1);
nByBin = zeros(nBins,1);

preVals = [];
postVals = [];

for tb = 1:nBins
    [~, ~, vSigned, stream] = smooth_post_affine_bin_knn(bins(tb), K, enforceK, siteWeightsByIndex);
    if isempty(vSigned)
        continue;
    end

    A = route_by_stream(vSigned, stream);
    A = A(isfinite(A));
    if isempty(A)
        continue;
    end

    fracByBin(tb) = mean(A > threshold);
    nByBin(tb) = numel(A);

    if preMask(tb)
        preVals = [preVals; A]; %#ok<AGROW>
    end
    if postMask(tb)
        postVals = [postVals; A]; %#ok<AGROW>
    end
end

S = struct();
S.K = K;
S.threshold = threshold;
S.fracByBin = fracByBin;
S.nByBin = nByBin;
S.preFrac = nan;
S.postFrac = nan;
S.preN = numel(preVals);
S.postN = numel(postVals);

if ~isempty(preVals)
    S.preFrac = mean(preVals > threshold);
end
if ~isempty(postVals)
    S.postFrac = mean(postVals > threshold);
end

end

function V = route_by_stream(vSigned, stream)
V = zeros(size(vSigned));
isT = (stream == 1);
isD = (stream == 2);
V(isT) = max(0,  vSigned(isT));
V(isD) = max(0, -vSigned(isD));
isOther = ~(isT | isD);
V(isOther) = abs(vSigned(isOther));
end
