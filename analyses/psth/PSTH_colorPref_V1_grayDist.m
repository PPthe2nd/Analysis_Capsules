%% PSTH_colorPref_V1_grayDist50_perStim.m
% Balanced time course: preferred vs non-preferred color (10 ms bins)
% + Adds an UNBALANCED gray (background) trace with a *per-stimulus* distance criterion:
%
%   - Preferred / Nonpreferred curves:
%       Selected ONLY by CI (+ optional SNR) and computed with BALANCED complementary pairs.
%
%   - Gray/background curve:
%       For each site, average activity ONLY across those stimuli where:
%           (RF is on gray) AND (dist_to_nearest_color_px >= minDistThr).
%       A site is excluded from the gray curve ONLY if it has zero qualifying gray stimuli.
%       (This has ZERO influence on the preferred/nonpreferred curves.)
%
% Requires (in workspace):
%   R.meanAct     [nCh x 384 x nBins]
%   R.meanSqAct   [nCh x 384 x nBins]
%   R.nTrials     [384 x 1]
%   R.timeWindows [nBins x 2]
%   Tall_V1(384x1).T : table with columns:
%       - center_color (string/cellstr)
%       - dist_to_nearest_color_px (numeric)  [only meaningful for gray RFs]
%   ColorTune.early.colorIndex or ColorTune.late.colorIndex (512x1)
%   ciThr  (scalar)
% Optional (if useSNRgate=true):
%   bestSNR (512x1), SNRthr (scalar)
%
% Output:
%   Figure with mean±SEM for preferred, nonpreferred, and gray(background) curves.

%% -----------------------
% SETTINGS
% -----------------------
useSNRgate = true;            % set false if you do not want SNR gating

COL_Y = "yellowArm";
COL_P = "purple";
COL_G = "gray";

%% -----------------------
% BASIC DIMS
% -----------------------
[nCh, nStim2, nBins] = size(R.meanAct);
nStim = nStim2;
assert(nStim == 384, 'Expected 384 stimuli.');
assert(size(R.timeWindows,1) == nBins, 'R.timeWindows rows must equal #bins.');

nTrialsRaw = double(R.nTrials);
if isvector(nTrialsRaw)
    nTrials = nTrialsRaw(:);
    perSiteTrials = false;
elseif ismatrix(nTrialsRaw) && size(nTrialsRaw,2) == nStim
    nTrials = [];
    perSiteTrials = true;
else
    error('R.nTrials must be a stimulus vector or a site-by-stimulus matrix.');
end

% V1 sites: rows 1..512 in R, and rows 1..512 in Tall_V1(stim).T
v1Sites = 1:512;

% time axis (bin centers)
tCenters = mean(R.timeWindows, 2);   % nBins x 1

%% -----------------------
% Choose which colorIndex to use
% -----------------------
switch lower(useColorIndexFrom)
    case "early"
        ci = ColorTune.early.colorIndex(1:512);  % V1 rows only
    case "late"
        ci = ColorTune.late.colorIndex(1:512);   % V1 rows only
    otherwise
        error('useColorIndexFrom must be "early" or "late".');
end
ci = ci(:);

% Preferred color per site: CI = (Y-P)/... ; positive -> prefers yellow
prefIsYellow = ci > 0;

%% -----------------------
% Initial selection for preferred/nonpreferred curves (CI + optional SNR)
% -----------------------
keep = isfinite(ci) & (abs(ci) > ciThr);

if useSNRgate
    snrThr = SNRthr; % expects SNRthr and bestSNR in workspace
    bestSNR_v1 = bestSNR(1:512);
    keep = keep & isfinite(bestSNR_v1(:)) & (bestSNR_v1(:) > snrThr);
end

keepColor = keep;
keepSitesColor = find(keepColor);

fprintf('Color selection: %d / 512 sites (abs(CI)>%.2f)%s\n', ...
    numel(keepSitesColor), ciThr, ternary(useSNRgate, sprintf(' & bestSNR>%.2f', snrThr), ''));

%% -----------------------
% Sort Tall_V1 by stimNum (safety)
% -----------------------
TallStimNums = arrayfun(@(x) x.stimNum, Tall_V1(:));
[sortedStimNums, order] = sort(TallStimNums(:));
assert(all(sortedStimNums(:).' == 1:nStim), 'Tall_V1.stimNum should cover 1..384.');
TallSorted = Tall_V1(order);

% Identify required columns in the table
T0 = TallSorted(1).T;
assert(istable(T0), 'Tall_V1(stim).T must be a table.');
vn = string(T0.Properties.VariableNames);

ccIdx = find(vn == "center_color", 1);
if isempty(ccIdx)
    ccIdx = find(contains(lower(vn), "center") & contains(lower(vn), "color"), 1);
end
assert(~isempty(ccIdx), 'Could not find center_color column.');

distIdx = find(vn == "dist_to_nearest_color_px", 1);
assert(~isempty(distIdx), 'Could not find dist_to_nearest_color_px column.');

% Build matrices: CC(site,stim) and Dist(site,stim)
CC   = strings(512, nStim);
Dist = nan(512, nStim);

for stim = 1:nStim
    Ti = TallSorted(stim).T;
    CC(:,stim)   = strtrim(string(Ti{:, ccIdx}));
    Dist(:,stim) = double(Ti{:, distIdx});
end

%% -----------------------
% Per-stimulus distance-qualified gray mask
% -----------------------
% qualGray(site,stim) true iff that site is on gray AND has distance>=minDistThr
qualGray = (CC == COL_G) & isfinite(Dist) & (Dist >= minDistThr);

% Sites that will contribute to the gray curve: must be in keepColor, and have >=1 qualifying gray stimulus
nQualGray = sum(qualGray, 2);
keepGraySite = keepColor & (nQualGray > 0);
keepSitesGray = find(keepGraySite);

fprintf('Gray selection (for gray curve ONLY): %d / %d color-selected sites have >=1 gray stim with dist>=%.1f px\n', ...
    numel(keepSitesGray), numel(keepSitesColor), minDistThr);

%% -----------------------
% Complement pairs within each block of 8: (1<->5),(2<->6),(3<->7),(4<->8), repeat
% We'll store only "first half" (pos 1..4) and pair with +4
% -----------------------
pairsA = zeros(nStim/2,1);
pairsB = zeros(nStim/2,1);
k = 0;
for a = 1:nStim
    pos = mod(a-1, 8) + 1;
    if pos <= 4
        k = k + 1;
        pairsA(k) = a;
        pairsB(k) = a + 4;
    end
end
nPairs = k;  % should be 192

%% -----------------------
% Compute per-site time courses
% -----------------------
muPref = nan(numel(keepSitesColor), nBins);
muNon  = nan(numel(keepSitesColor), nBins);
muGray = nan(numel(keepSitesColor), nBins);  % NaN for sites not in keepGraySite or with no qualifying trials

for ii = 1:numel(keepSitesColor)
    site = keepSitesColor(ii);  % 1..512 index into Tall / V1 list
    ch   = v1Sites(site);       % row in R.meanAct (here identity)

    if perSiteTrials
        nTrSite = double(nTrialsRaw(ch,:)).';
    else
        nTrSite = nTrials;
    end

    % Balanced accumulators for Y/P
    sumY = zeros(nBins,1); NY = 0;
    sumP = zeros(nBins,1); NP = 0;

    % Gray accumulators (unbalanced; but only over qualifying stimuli)
    sumG = zeros(nBins,1); NG = 0;
    doGray = keepGraySite(site);

    for ip = 1:nPairs
        a = pairsA(ip);
        b = pairsB(ip);

        na = nTrSite(a); nb = nTrSite(b);

        % time courses for these stimuli (nBins x 1)
        ma = squeeze(R.meanAct(ch,a,:));
        mb = squeeze(R.meanAct(ch,b,:));

        ca = CC(site,a);
        cb = CC(site,b);

        % ---- Gray accumulation: per-stimulus qualified mask (NO pair balancing) ----
        if doGray
            if qualGray(site, a) && (na > 0)
                sumG = sumG + na * ma;
                NG = NG + na;
            end
            if qualGray(site, b) && (nb > 0)
                sumG = sumG + nb * mb;
                NG = NG + nb;
            end
        end

        % ---- Balanced Y/P accumulation (pair-balanced) ----
        nEff = min(na, nb);
        if nEff <= 0
            continue;
        end

        % stimulus a
        if ca == COL_Y
            sumY = sumY + nEff * ma; NY = NY + nEff;
        elseif ca == COL_P
            sumP = sumP + nEff * ma; NP = NP + nEff;
        end

        % stimulus b
        if cb == COL_Y
            sumY = sumY + nEff * mb; NY = NY + nEff;
        elseif cb == COL_P
            sumP = sumP + nEff * mb; NP = NP + nEff;
        end
    end

    % Store gray mean if any qualifying trials
    if doGray && (NG > 0)
        muGray(ii,:) = (sumG / NG).';
    end

    % Store preferred/nonpreferred if both colors have data
    if (NY < 1) || (NP < 1)
        continue;
    end

    muY = (sumY / NY).'; % 1 x nBins
    muP = (sumP / NP).'; % 1 x nBins

    if prefIsYellow(site)
        muPref(ii,:) = muY;
        muNon(ii,:)  = muP;
    else
        muPref(ii,:) = muP;
        muNon(ii,:)  = muY;
    end
end

% Drop sites that ended up NaN everywhere for pref or nonpref (this does NOT consider gray)
good = any(isfinite(muPref), 2) & any(isfinite(muNon), 2);

muPref = muPref(good,:);
muNon  = muNon(good,:);
muGray = muGray(good,:);

fprintf('After requiring usable balanced Y/P trials: N=%d sites (pref/nonpref curves)\n', size(muPref,1));
fprintf('Gray curve uses subset (via NaNs): N=%d sites with any qualifying gray trials\n', sum(any(isfinite(muGray),2)));

%% -----------------------
% Plot mean ± SEM over sites
% -----------------------
mPref = mean(muPref, 1, 'omitnan');
mNon  = mean(muNon,  1, 'omitnan');
mGray = mean(muGray, 1, 'omitnan');

semPref = std(muPref, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(muPref),1));
semNon  = std(muNon,  0, 1, 'omitnan') ./ sqrt(sum(isfinite(muNon),1));
semGray = std(muGray, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(muGray),1));

figure; hold on;

plot(tCenters, mPref, 'LineWidth', 2);
plot(tCenters, mNon,  'LineWidth', 2);
plot(tCenters, mGray, 'LineWidth', 2);

% Shaded error (simple patch; neutral color)
fill([tCenters; flipud(tCenters)], [ (mPref-semPref)'; flipud((mPref+semPref)') ], ...
     'k', 'FaceAlpha', 0.12, 'EdgeColor', 'none');
fill([tCenters; flipud(tCenters)], [ (mNon-semNon)';  flipud((mNon+semNon)')  ], ...
     'k', 'FaceAlpha', 0.06, 'EdgeColor', 'none');
fill([tCenters; flipud(tCenters)], [ (mGray-semGray)'; flipud((mGray+semGray)') ], ...
     'k', 'FaceAlpha', 0.04, 'EdgeColor', 'none');

xline(0,'k-');
xlabel('Time from stimulus onset (ms)');
ylabel('Mean response (a.u.)');

ttl = sprintf('Balanced Pref vs Nonpref + Gray(dist>=%.0f px) | abs(CI)>%.2f', minDistThr, ciThr);
if useSNRgate
    ttl = ttl + ", bestSNR>" + snrThr;
end
title(ttl);

legend('Preferred','Nonpreferred',sprintf('Gray (dist >= %.0f px)', minDistThr),'Location','best');
grid on;

%% Helper: inline ternary
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
