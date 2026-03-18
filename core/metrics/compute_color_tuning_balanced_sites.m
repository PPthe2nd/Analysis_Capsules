function ColorTune = compute_color_tuning_balanced_sites(DATA, Tall, siteRows, keepSiteIdx, varargin)
% COMPUTE_COLOR_TUNING_BALANCED_SITES
% Compute balanced yellow-vs-purple tuning metrics for an arbitrary site set.
%
% DATA must contain:
%   DATA.meanAct   [nSites x 384 x >=3]
%   DATA.meanSqAct [nSites x 384 x >=3]
%   DATA.nTrials   [384 x 1] or [nSites x 384]
%
% Tall(stim).T must contain center_color and enough rows for siteRows.
%
% siteRows    : row indices available in DATA/Tall.
% keepSiteIdx : indices within siteRows to actually compute.
%
% Outputs mirror the V1 ColorTune struct and add a paired precision-weighted
% complementary test:
%   ColorTune.early.* and ColorTune.late.* with fields
%   muY, muP, varY, varP, NY, NP, colorIndex, dprime, z, p,
%   pairedMuY, pairedMuP, pairedMeanDiff, pairedWeightedDiff,
%   pairedWeightSum, pairedEffN, pairedT, pairedP, pairedNPairs

p = inputParser;
p.addParameter('WIN_EARLY', 2, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('WIN_LATE', 3, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('useBesselCorrection', true, @(x) islogical(x) && isscalar(x));
p.addParameter('varFloorFrac', 1e-3, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('Verbose', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;

assert(isfield(DATA, 'meanAct') && isfield(DATA, 'meanSqAct') && isfield(DATA, 'nTrials'), ...
    'DATA must contain meanAct, meanSqAct, and nTrials.');

[nSitesTotal, nStim, nWin] = size(DATA.meanAct);
assert(all(size(DATA.meanSqAct) == size(DATA.meanAct)), ...
    'DATA.meanSqAct must match DATA.meanAct.');
assert(nStim == 384, 'Expected 384 stimuli, got %d.', nStim);
assert(nWin >= max([opt.WIN_EARLY, opt.WIN_LATE]), ...
    'DATA does not contain the requested time windows.');
assert(numel(Tall) == nStim, 'Tall has %d entries but DATA has %d stimuli.', numel(Tall), nStim);

siteRows = siteRows(:);
keepSiteIdx = keepSiteIdx(:);
assert(~isempty(siteRows), 'siteRows must not be empty.');
assert(all(siteRows >= 1 & siteRows <= nSitesTotal), ...
    'siteRows must lie within DATA rows 1..%d.', nSitesTotal);
assert(all(keepSiteIdx >= 1 & keepSiteIdx <= numel(siteRows)), ...
    'keepSiteIdx must lie within 1..%d.', numel(siteRows));

keepRows = siteRows(keepSiteIdx);

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

TallStimNums = arrayfun(@(x) x.stimNum, Tall(:));
[sortedStimNums, order] = sort(TallStimNums(:));
assert(all(sortedStimNums(:).' == 1:nStim), 'Tall.stimNum should cover 1..%d exactly.', nStim);
TallSorted = Tall(order);

T0 = TallSorted(1).T;
assert(istable(T0), 'Tall(stim).T must be a table.');
assert(height(T0) >= max(siteRows), ...
    'Tall(stim).T has %d rows; need at least %d.', height(T0), max(siteRows));

varNames = string(T0.Properties.VariableNames);
colIdx = find(varNames == "center_color", 1);
if isempty(colIdx)
    colIdx = find(contains(lower(varNames), "center") & contains(lower(varNames), "color"), 1);
end
assert(~isempty(colIdx), 'Could not find center_color column in Tall(stim).T.');

COL_Y = "yellowArm";
COL_P = "purple";
COL_G = "gray";

nSel = numel(siteRows);
CC = strings(nSel, nStim);
for stim = 1:nStim
    Ttbl = TallSorted(stim).T;
    labs = string(Ttbl{siteRows, colIdx});
    CC(:, stim) = strtrim(labs);
end

pairsA = zeros(192,1);
pairsB = zeros(192,1);
ip = 0;
for stimA = 1:nStim
    pos = mod(stimA - 1, 8) + 1;
    if pos <= 4
        ip = ip + 1;
        pairsA(ip) = stimA;
        pairsB(ip) = stimA + 4;
    end
end
assert(ip == 192, 'Expected 192 complementary pairs.');

ColorTune = struct();
ColorTune.siteRows = siteRows;
ColorTune.keepSites = keepRows;
ColorTune.keepSiteIdx = keepSiteIdx;

fields = ["muY","muP","varY","varP","NY","NP","colorIndex","dprime","z","p", ...
    "pairedMuY","pairedMuP","pairedMeanDiff","pairedWeightedDiff", ...
    "pairedWeightSum","pairedEffN","pairedT","pairedP","pairedNPairs"];
for f = fields
    ColorTune.early.(f) = nan(nSitesTotal,1);
    ColorTune.late.(f)  = nan(nSitesTotal,1);
end

rowToLocal = nan(nSitesTotal,1);
rowToLocal(siteRows) = 1:nSel;

for iSite = keepRows(:).'
    localIdx = rowToLocal(iSite);
    if ~isfinite(localIdx)
        continue;
    end

    mEarly   = squeeze(DATA.meanAct(iSite, :, opt.WIN_EARLY));  mEarly   = mEarly(:);
    msqEarly = squeeze(DATA.meanSqAct(iSite, :, opt.WIN_EARLY)); msqEarly = msqEarly(:);
    mLate    = squeeze(DATA.meanAct(iSite, :, opt.WIN_LATE));   mLate    = mLate(:);
    msqLate  = squeeze(DATA.meanSqAct(iSite, :, opt.WIN_LATE)); msqLate  = msqLate(:);

    sumY_e = 0; sumY2_e = 0; NY_e = 0;
    sumP_e = 0; sumP2_e = 0; NP_e = 0;
    sumY_l = 0; sumY2_l = 0; NY_l = 0;
    sumP_l = 0; sumP2_l = 0; NP_l = 0;
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

        ca = CC(localIdx, a);
        cb = CC(localIdx, b);

        if (ca == COL_G) && (cb == COL_G)
            continue;
        end
        if (ca == "" || ismissing(ca) || cb == "" || ismissing(cb))
            continue;
        end

        if ca == COL_Y
            sumY_e  = sumY_e  + nEff * mEarly(a);
            sumY2_e = sumY2_e + nEff * msqEarly(a);
            NY_e    = NY_e    + nEff;
            sumY_l  = sumY_l  + nEff * mLate(a);
            sumY2_l = sumY2_l + nEff * msqLate(a);
            NY_l    = NY_l    + nEff;
        elseif ca == COL_P
            sumP_e  = sumP_e  + nEff * mEarly(a);
            sumP2_e = sumP2_e + nEff * msqEarly(a);
            NP_e    = NP_e    + nEff;
            sumP_l  = sumP_l  + nEff * mLate(a);
            sumP2_l = sumP2_l + nEff * msqLate(a);
            NP_l    = NP_l    + nEff;
        end

        if cb == COL_Y
            sumY_e  = sumY_e  + nEff * mEarly(b);
            sumY2_e = sumY2_e + nEff * msqEarly(b);
            NY_e    = NY_e    + nEff;
            sumY_l  = sumY_l  + nEff * mLate(b);
            sumY2_l = sumY2_l + nEff * msqLate(b);
            NY_l    = NY_l    + nEff;
        elseif cb == COL_P
            sumP_e  = sumP_e  + nEff * mEarly(b);
            sumP2_e = sumP2_e + nEff * msqEarly(b);
            NP_e    = NP_e    + nEff;
            sumP_l  = sumP_l  + nEff * mLate(b);
            sumP2_l = sumP2_l + nEff * msqLate(b);
            NP_l    = NP_l    + nEff;
        end

        if ca == COL_Y && cb == COL_P
            [dEarly, vEarly, muYEarly, muPEarly] = build_paired_diff(mEarly(a), msqEarly(a), na, ...
                mEarly(b), msqEarly(b), nb, nEff, opt.useBesselCorrection);
            [dLate, vLate, muYLate, muPLate] = build_paired_diff(mLate(a), msqLate(a), na, ...
                mLate(b), msqLate(b), nb, nEff, opt.useBesselCorrection);
        elseif ca == COL_P && cb == COL_Y
            [dEarly, vEarly, muYEarly, muPEarly] = build_paired_diff(mEarly(b), msqEarly(b), nb, ...
                mEarly(a), msqEarly(a), na, nEff, opt.useBesselCorrection);
            [dLate, vLate, muYLate, muPLate] = build_paired_diff(mLate(b), msqLate(b), nb, ...
                mLate(a), msqLate(a), na, nEff, opt.useBesselCorrection);
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

    [muY, varY, muP, varP, NY, NP, cind, dp, z, pVal] = ...
        compute_metrics(sumY_e, sumY2_e, NY_e, sumP_e, sumP2_e, NP_e, opt.useBesselCorrection);
    [pairMuY, pairMuP, pairMeanDiff, pairWeightedDiff, pairWeightSum, pairEffN, pairT, pairP] = ...
        compute_paired_stats(pairDiffEarly(1:nPairEarly), pairVarEarly(1:nPairEarly), ...
        pairYEarly(1:nPairEarly), pairPEarly(1:nPairEarly), opt.varFloorFrac);
    ColorTune.early.muY(iSite) = muY;
    ColorTune.early.varY(iSite) = varY;
    ColorTune.early.muP(iSite) = muP;
    ColorTune.early.varP(iSite) = varP;
    ColorTune.early.NY(iSite) = NY;
    ColorTune.early.NP(iSite) = NP;
    ColorTune.early.colorIndex(iSite) = cind;
    ColorTune.early.dprime(iSite) = dp;
    ColorTune.early.z(iSite) = z;
    ColorTune.early.p(iSite) = pVal;
    ColorTune.early.pairedMuY(iSite) = pairMuY;
    ColorTune.early.pairedMuP(iSite) = pairMuP;
    ColorTune.early.pairedMeanDiff(iSite) = pairMeanDiff;
    ColorTune.early.pairedWeightedDiff(iSite) = pairWeightedDiff;
    ColorTune.early.pairedWeightSum(iSite) = pairWeightSum;
    ColorTune.early.pairedEffN(iSite) = pairEffN;
    ColorTune.early.pairedT(iSite) = pairT;
    ColorTune.early.pairedP(iSite) = pairP;
    ColorTune.early.pairedNPairs(iSite) = nPairEarly;

    [muY, varY, muP, varP, NY, NP, cind, dp, z, pVal] = ...
        compute_metrics(sumY_l, sumY2_l, NY_l, sumP_l, sumP2_l, NP_l, opt.useBesselCorrection);
    [pairMuY, pairMuP, pairMeanDiff, pairWeightedDiff, pairWeightSum, pairEffN, pairT, pairP] = ...
        compute_paired_stats(pairDiffLate(1:nPairLate), pairVarLate(1:nPairLate), ...
        pairYLate(1:nPairLate), pairPLate(1:nPairLate), opt.varFloorFrac);
    ColorTune.late.muY(iSite) = muY;
    ColorTune.late.varY(iSite) = varY;
    ColorTune.late.muP(iSite) = muP;
    ColorTune.late.varP(iSite) = varP;
    ColorTune.late.NY(iSite) = NY;
    ColorTune.late.NP(iSite) = NP;
    ColorTune.late.colorIndex(iSite) = cind;
    ColorTune.late.dprime(iSite) = dp;
    ColorTune.late.z(iSite) = z;
    ColorTune.late.p(iSite) = pVal;
    ColorTune.late.pairedMuY(iSite) = pairMuY;
    ColorTune.late.pairedMuP(iSite) = pairMuP;
    ColorTune.late.pairedMeanDiff(iSite) = pairMeanDiff;
    ColorTune.late.pairedWeightedDiff(iSite) = pairWeightedDiff;
    ColorTune.late.pairedWeightSum(iSite) = pairWeightSum;
    ColorTune.late.pairedEffN(iSite) = pairEffN;
    ColorTune.late.pairedT(iSite) = pairT;
    ColorTune.late.pairedP(iSite) = pairP;
    ColorTune.late.pairedNPairs(iSite) = nPairLate;
end

if opt.Verbose
    fprintf('Done. Example: median d'' early (kept) = %.3f\n', ...
        median(ColorTune.early.dprime(keepRows), 'omitnan'));
end

function [pairDiff, pairVar, muY, muP] = build_paired_diff(muY, msqY, nOrigY, muP, msqP, nOrigP, nEff, useBessel)
pairDiff = NaN;
pairVar = NaN;
if ~(isfinite(muY) && isfinite(msqY) && isfinite(nOrigY) && nOrigY > 1 && ...
        isfinite(muP) && isfinite(msqP) && isfinite(nOrigP) && nOrigP > 1 && ...
        isfinite(nEff) && nEff > 0)
    return;
end

varMeanY = variance_of_balanced_mean(muY, msqY, nOrigY, nEff, useBessel);
varMeanP = variance_of_balanced_mean(muP, msqP, nOrigP, nEff, useBessel);
if ~(isfinite(varMeanY) && isfinite(varMeanP))
    return;
end

pairDiff = muY - muP;
pairVar = varMeanY + varMeanP;
end

function varMean = variance_of_balanced_mean(mu, msq, nOrig, nEff, useBessel)
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
        compute_paired_stats(pairDiff, pairVar, pairY, pairP, varFloorFrac)
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
    varFloor = 1e-6;
else
    varFloor = max(median(posVar) * varFloorFrac, 1e-6);
end
w = 1 ./ max(pairVar, varFloor);
sumW = sum(w);
if ~(isfinite(sumW) && sumW > 0)
    return;
end
meanDiff = mean(pairDiff);
weightedDiff = sum(w .* pairDiff) / sumW;
muYw = sum(w .* pairY) / sumW;
muPw = sum(w .* pairP) / sumW;
effN = (sumW^2) / sum(w.^2);

df = numel(pairDiff) - 1;
resid = pairDiff - weightedDiff;
rssWeighted = sum(w .* (resid .^ 2));
if ~(isfinite(rssWeighted) && df > 0)
    return;
end

sigma2Hat = rssWeighted / df;
seWeighted = sqrt(sigma2Hat / sumW);
if isfinite(seWeighted) && seWeighted > 0
    tStat = weightedDiff / seWeighted;
    pVal = 1 - fcdf_local(tStat.^2, 1, df);
elseif weightedDiff ~= 0
    tStat = sign(weightedDiff) * Inf;
    pVal = 0;
else
    tStat = 0;
    pVal = 1;
end
end

function p = fcdf_local(x, df1, df2)
z = (df1 .* x) ./ (df1 .* x + df2);
z = min(max(z, 0), 1);
p = betainc(z, df1/2, df2/2);
end

end

function [muY,varY,muP,varP,NY,NP,colorIndex,dprime,z,p] = compute_metrics(sumY,sumY2,NY,sumP,sumP2,NP,useBessel)
muY = NaN; varY = NaN; muP = NaN; varP = NaN;
colorIndex = NaN; dprime = NaN; z = NaN; p = NaN;

if NY <= 1 || NP <= 1
    return;
end

muY = sumY / NY;
muP = sumP / NP;

Ex2Y = sumY2 / NY;
Ex2P = sumP2 / NP;

varY = max(0, Ex2Y - muY^2);
varP = max(0, Ex2P - muP^2);

if useBessel
    varY = varY * (NY / (NY - 1));
    varP = varP * (NP / (NP - 1));
end

denom = abs((muY + muP) / 2);
if denom > 0
    colorIndex = (muY - muP) / denom;
end

denomDP = sqrt(0.5 * (varY + varP));
if denomDP > 0
    dprime = (muY - muP) / denomDP;
end

denomZ = sqrt(varY / NY + varP / NP);
if denomZ > 0
    z = (muY - muP) / denomZ;
    p = 2 * normcdf(-abs(z), 0, 1);
end
end
