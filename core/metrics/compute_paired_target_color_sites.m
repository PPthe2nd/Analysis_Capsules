function TargetColor = compute_paired_target_color_sites(DATA, targetColorByStim, pairsA, pairsB, varargin)
% COMPUTE_PAIRED_TARGET_COLOR_SITES
% Paired precision-weighted target-color metric for an arbitrary site set.
%
% The logic mirrors the paired IT color analysis, but the stimulus label is
% the color of the target capsule rather than the color in the RF.
%
% Inputs:
%   DATA.meanAct   [nSites x 384 x >=3]
%   DATA.meanSqAct [nSites x 384 x >=3]
%   DATA.nTrials   [384 x 1] or [nSites x 384]
%   targetColorByStim [384 x 1] string, each element "yellowArm" or "purple"
%   pairsA / pairsB   complementary pair indices
%
% Outputs:
%   TargetColor.early.* and TargetColor.late.* with fields
%   pairedMuTargetYellow, pairedMuTargetPurple, pairedMeanDiff,
%   pairedWeightedDiff, pairedWeightSum, pairedEffN, pairedT, pairedP,
%   pairedNPairs

p = inputParser;
p.addParameter('WIN_EARLY', 2, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('WIN_LATE', 3, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('useBesselCorrection', true, @(x) islogical(x) && isscalar(x));
p.addParameter('varFloorFrac', 1e-3, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.parse(varargin{:});
opt = p.Results;

assert(isfield(DATA, 'meanAct') && isfield(DATA, 'meanSqAct') && isfield(DATA, 'nTrials'), ...
    'DATA must contain meanAct, meanSqAct, and nTrials.');

[nSitesTotal, nStim, nWin] = size(DATA.meanAct);
assert(all(size(DATA.meanSqAct) == size(DATA.meanAct)), ...
    'DATA.meanSqAct must match DATA.meanAct.');
assert(nStim == numel(targetColorByStim), ...
    'Stimulus count mismatch between responses and targetColorByStim.');
assert(nWin >= max([opt.WIN_EARLY, opt.WIN_LATE]), ...
    'DATA does not contain the requested time windows.');

nTrials = DATA.nTrials;
if isvector(nTrials)
    nTrials = double(nTrials(:));
    assert(numel(nTrials) == nStim, 'DATA.nTrials vector must have %d elements.', nStim);
    perSiteTrials = false;
elseif ismatrix(nTrials) && size(nTrials,1) == nSitesTotal && size(nTrials,2) == nStim
    nTrials = double(nTrials);
    perSiteTrials = true;
else
    error('DATA.nTrials must be a vector with %d elements or a matrix [nSites x %d].', nStim, nStim);
end

fields = ["pairedMuTargetYellow","pairedMuTargetPurple","pairedMeanDiff", ...
    "pairedWeightedDiff","pairedWeightSum","pairedEffN","pairedT","pairedP","pairedNPairs"];
for f = fields
    TargetColor.early.(f) = nan(nSitesTotal,1);
    TargetColor.late.(f)  = nan(nSitesTotal,1);
end

targetColorByStim = string(targetColorByStim(:));

for iSite = 1:nSitesTotal
    mEarly = squeeze(DATA.meanAct(iSite, :, opt.WIN_EARLY));    mEarly = mEarly(:);
    msqEarly = squeeze(DATA.meanSqAct(iSite, :, opt.WIN_EARLY)); msqEarly = msqEarly(:);
    mLate = squeeze(DATA.meanAct(iSite, :, opt.WIN_LATE));      mLate = mLate(:);
    msqLate = squeeze(DATA.meanSqAct(iSite, :, opt.WIN_LATE));  msqLate = msqLate(:);

    pairDiffEarly = nan(numel(pairsA), 1);
    pairVarEarly = nan(numel(pairsA), 1);
    pairYEarly = nan(numel(pairsA), 1);
    pairPEarly = nan(numel(pairsA), 1);
    pairDiffLate = nan(numel(pairsA), 1);
    pairVarLate = nan(numel(pairsA), 1);
    pairYLate = nan(numel(pairsA), 1);
    pairPLate = nan(numel(pairsA), 1);
    nPairEarly = 0;
    nPairLate = 0;

    for k = 1:numel(pairsA)
        a = pairsA(k);
        b = pairsB(k);

        if perSiteTrials
            na = nTrials(iSite, a);
            nb = nTrials(iSite, b);
        else
            na = nTrials(a);
            nb = nTrials(b);
        end

        nEff = min(na, nb);
        if ~(isfinite(nEff) && nEff > 0)
            continue;
        end

        ca = targetColorByStim(a);
        cb = targetColorByStim(b);

        if ca == "yellowArm" && cb == "purple"
            [dEarly, vEarly, muYEarly, muPEarly] = build_paired_diff_local( ...
                mEarly(a), msqEarly(a), na, mEarly(b), msqEarly(b), nb, nEff, opt.useBesselCorrection);
            [dLate, vLate, muYLate, muPLate] = build_paired_diff_local( ...
                mLate(a), msqLate(a), na, mLate(b), msqLate(b), nb, nEff, opt.useBesselCorrection);
        elseif ca == "purple" && cb == "yellowArm"
            [dEarly, vEarly, muYEarly, muPEarly] = build_paired_diff_local( ...
                mEarly(b), msqEarly(b), nb, mEarly(a), msqEarly(a), na, nEff, opt.useBesselCorrection);
            [dLate, vLate, muYLate, muPLate] = build_paired_diff_local( ...
                mLate(b), msqLate(b), nb, mLate(a), msqLate(a), na, nEff, opt.useBesselCorrection);
        else
            dEarly = NaN; vEarly = NaN; muYEarly = NaN; muPEarly = NaN;
            dLate = NaN; vLate = NaN; muYLate = NaN; muPLate = NaN;
        end

        if isfinite(dEarly) && isfinite(vEarly)
            nPairEarly = nPairEarly + 1;
            pairDiffEarly(nPairEarly) = dEarly;
            pairVarEarly(nPairEarly) = vEarly;
            pairYEarly(nPairEarly) = muYEarly;
            pairPEarly(nPairEarly) = muPEarly;
        end
        if isfinite(dLate) && isfinite(vLate)
            nPairLate = nPairLate + 1;
            pairDiffLate(nPairLate) = dLate;
            pairVarLate(nPairLate) = vLate;
            pairYLate(nPairLate) = muYLate;
            pairPLate(nPairLate) = muPLate;
        end
    end

    [pairMuY, pairMuP, pairMeanDiff, pairWeightedDiff, pairWeightSum, pairEffN, pairT, pairP] = ...
        compute_paired_stats_local(pairDiffEarly(1:nPairEarly), pairVarEarly(1:nPairEarly), ...
        pairYEarly(1:nPairEarly), pairPEarly(1:nPairEarly), opt.varFloorFrac);
    TargetColor.early.pairedMuTargetYellow(iSite) = pairMuY;
    TargetColor.early.pairedMuTargetPurple(iSite) = pairMuP;
    TargetColor.early.pairedMeanDiff(iSite) = pairMeanDiff;
    TargetColor.early.pairedWeightedDiff(iSite) = pairWeightedDiff;
    TargetColor.early.pairedWeightSum(iSite) = pairWeightSum;
    TargetColor.early.pairedEffN(iSite) = pairEffN;
    TargetColor.early.pairedT(iSite) = pairT;
    TargetColor.early.pairedP(iSite) = pairP;
    TargetColor.early.pairedNPairs(iSite) = nPairEarly;

    [pairMuY, pairMuP, pairMeanDiff, pairWeightedDiff, pairWeightSum, pairEffN, pairT, pairP] = ...
        compute_paired_stats_local(pairDiffLate(1:nPairLate), pairVarLate(1:nPairLate), ...
        pairYLate(1:nPairLate), pairPLate(1:nPairLate), opt.varFloorFrac);
    TargetColor.late.pairedMuTargetYellow(iSite) = pairMuY;
    TargetColor.late.pairedMuTargetPurple(iSite) = pairMuP;
    TargetColor.late.pairedMeanDiff(iSite) = pairMeanDiff;
    TargetColor.late.pairedWeightedDiff(iSite) = pairWeightedDiff;
    TargetColor.late.pairedWeightSum(iSite) = pairWeightSum;
    TargetColor.late.pairedEffN(iSite) = pairEffN;
    TargetColor.late.pairedT(iSite) = pairT;
    TargetColor.late.pairedP(iSite) = pairP;
    TargetColor.late.pairedNPairs(iSite) = nPairLate;
end
end

function [pairDiff, pairVar, muY, muP] = build_paired_diff_local(muY, msqY, nOrigY, muP, msqP, nOrigP, nEff, useBessel)
pairDiff = NaN;
pairVar = NaN;
if ~(isfinite(muY) && isfinite(msqY) && isfinite(nOrigY) && nOrigY > 1 && ...
        isfinite(muP) && isfinite(msqP) && isfinite(nOrigP) && nOrigP > 1 && ...
        isfinite(nEff) && nEff > 0)
    return;
end

varMeanY = variance_of_balanced_mean_local(muY, msqY, nOrigY, nEff, useBessel);
varMeanP = variance_of_balanced_mean_local(muP, msqP, nOrigP, nEff, useBessel);
if ~(isfinite(varMeanY) && isfinite(varMeanP))
    return;
end

pairDiff = muY - muP;
pairVar = varMeanY + varMeanP;
end

function varMean = variance_of_balanced_mean_local(mu, msq, nOrig, nEff, useBessel)
varMean = NaN;
if ~(isfinite(mu) && isfinite(msq) && isfinite(nOrig) && nOrig > 1 && isfinite(nEff) && nEff > 0)
    return;
end
sampleVar = max(0, msq - mu^2);
if useBessel
    sampleVar = sampleVar * (nOrig / (nOrig - 1));
end
varMean = sampleVar / nEff;
end

function [muYw, muPw, meanDiff, weightedDiff, sumW, effN, tStat, pVal] = ...
        compute_paired_stats_local(pairDiff, pairVar, pairY, pairP, varFloorFrac)
muYw = NaN; muPw = NaN; meanDiff = NaN; weightedDiff = NaN;
sumW = NaN; effN = NaN; tStat = NaN; pVal = NaN;

valid = isfinite(pairDiff) & isfinite(pairVar) & (pairVar >= 0) & isfinite(pairY) & isfinite(pairP);
pairDiff = pairDiff(valid);
pairVar = pairVar(valid);
pairY = pairY(valid);
pairP = pairP(valid);
if numel(pairDiff) < 2
    return;
end

posVar = pairVar(pairVar > 0);
if isempty(posVar)
    varFloor = 1;
else
    varFloor = varFloorFrac * median(posVar);
end
varAdj = max(pairVar, varFloor);
weights = 1 ./ varAdj;
sumW = sum(weights);
if ~(isfinite(sumW) && sumW > 0)
    return;
end

weightedDiff = sum(weights .* pairDiff) / sumW;
muYw = sum(weights .* pairY) / sumW;
muPw = sum(weights .* pairP) / sumW;
meanDiff = mean(pairDiff);

effN = (sumW^2) / sum(weights.^2);
seWeighted = sqrt(1 / sumW);
if ~(isfinite(seWeighted) && seWeighted > 0)
    return;
end
 tStat = weightedDiff / seWeighted;

if isfinite(effN) && (effN > 1)
    pVal = 2 * tcdf(-abs(tStat), effN - 1);
else
    pVal = NaN;
end
end
