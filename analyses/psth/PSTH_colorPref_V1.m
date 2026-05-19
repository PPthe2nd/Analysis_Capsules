%% Balanced time course: preferred vs non-preferred color (10 ms bins)
% + Adds an UNBALANCED gray (background) trace:
%   For each selected site, we average activity across ALL stimuli where that
%   site's RF is on gray, regardless of complementary pair logic.
%
% Requires:
%   R.meanAct     [nCh x 384 x nBins]    (you have 70 bins)
%   R.meanSqAct   [nCh x 384 x nBins]    (for variance/SEM if needed later)
%   R.nTrials     [384 x 1]
%   R.timeWindows [nBins x 2]            (e.g. -200..500 in 10 ms steps)
%   Tall_V1(384x1).T is a table with center_color per site (512 rows)
%   ColorTune.early.colorIndex or ColorTune.late.colorIndex (512x1)
%   bestSNR (512x1) if you enable SNR gating
%
% Notes:
% - Preferred/nonpreferred traces are computed with BALANCED complementary pairs.
% - Gray trace ignores pair balancing; it's a descriptive baseline for "RF on background".

% -----------------------
% SETTINGS
% -----------------------
useColorIndexFrom = "early";      % "early" or "late"
snrThr = SNRthr;                  % optional extra gate: require bestSNR > snrThr
useSNRgate = true;

COL_Y = "yellowArm";
COL_P = "purple";
COL_G = "gray";

% -----------------------
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

% -----------------------
% Choose which colorIndex to use for tuning selection + preference
% -----------------------
switch lower(useColorIndexFrom)
    case "early"
        ci = ColorTune.early.colorIndex(:);  % 512x1
    case "late"
        ci = ColorTune.late.colorIndex(:);   % 512x1
    otherwise
        error('useColorIndexFrom must be "early" or "late".');
end

% Optional: also require SNR threshold
if useSNRgate
    keep = isfinite(ci) & (abs(ci) > ciThr) & isfinite(bestSNR(:)) & (bestSNR(:) > snrThr);
else
    keep = isfinite(ci) & (abs(ci) > ciThr);
end

keepSites = find(keep);   % indices 1..512
fprintf('Selected %d / 512 sites (abs(CI)>%.2f)%s\n', numel(keepSites), ciThr, ...
    ternary(useSNRgate, sprintf(' & bestSNR>%.2f', snrThr), ''));

% Preferred color per site
prefIsYellow = ci > 0;    % CI = (Y-P)/... ; positive -> prefers yellow

% -----------------------
% Sort Tall_V1 by stimNum (safety)
% -----------------------
TallStimNums = arrayfun(@(x) x.stimNum, Tall_V1(:));
[sortedStimNums, order] = sort(TallStimNums(:));
assert(all(sortedStimNums(:).' == 1:nStim), 'Tall_V1.stimNum should cover 1..384.');
TallSorted = Tall_V1(order);

% Identify center_color column in the table
T0 = TallSorted(1).T;
assert(istable(T0), 'Tall_V1(stim).T must be a table.');
vn = string(T0.Properties.VariableNames);
colIdx = find(vn == "center_color", 1);
if isempty(colIdx)
    colIdx = find(contains(lower(vn), "center") & contains(lower(vn), "color"), 1);
end
assert(~isempty(colIdx), 'Could not find center_color column.');

% Build CC(site,stim): 512 x 384 string matrix
CC = strings(512, nStim);
for stim = 1:nStim
    labs = string(TallSorted(stim).T{:, colIdx});   % 512x1
    CC(:,stim) = strtrim(labs);
end

% -----------------------
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

% -----------------------
% Compute per-site time courses:
%   - Balanced muYellow(t), muPurple(t) using complementary pairs (min trial count per pair)
%   - Unbalanced muGray(t) using ALL stimuli where RF is on gray (uses full trial counts)
% Then map Y/P to preferred/nonpreferred per site.
% -----------------------
muPref = nan(numel(keepSites), nBins);
muNon  = nan(numel(keepSites), nBins);
muGray = nan(numel(keepSites), nBins);

for ii = 1:numel(keepSites)
    site = keepSites(ii);     % 1..512 index into Tall / V1 list
    ch   = v1Sites(site);     % row in R.meanAct (here identity)

    if perSiteTrials
        nTrSite = double(nTrialsRaw(ch,:)).';
    else
        nTrSite = nTrials;
    end

    % Balanced accumulators over time bins (Y/P)
    sumY  = zeros(nBins,1);  sumY2 = zeros(nBins,1);  NY = 0;
    sumP  = zeros(nBins,1);  sumP2 = zeros(nBins,1);  NP = 0;

    % Unbalanced gray accumulator (G)
    sumG  = zeros(nBins,1);  sumG2 = zeros(nBins,1);  NG = 0;

    for ip = 1:nPairs
        a = pairsA(ip);
        b = pairsB(ip);

        na = nTrSite(a); nb = nTrSite(b);

        % time courses for these stimuli (nBins x 1)
        ma  = squeeze(R.meanAct(ch,a,:));
        mb  = squeeze(R.meanAct(ch,b,:));
        m2a = squeeze(R.meanSqAct(ch,a,:));
        m2b = squeeze(R.meanSqAct(ch,b,:));

        ca = CC(site,a);
        cb = CC(site,b);

        % ---- Unbalanced gray accumulation (ignores pair balancing) ----
        if (ca == COL_G) && (na > 0)
            sumG  = sumG  + na * ma;
            sumG2 = sumG2 + na * m2a;
            NG = NG + na;
        end
        if (cb == COL_G) && (nb > 0)
            sumG  = sumG  + nb * mb;
            sumG2 = sumG2 + nb * m2b;
            NG = NG + nb;
        end

        % ---- Balanced Y/P accumulation (pair-balanced) ----
        nEff = min(na, nb);
        if nEff <= 0
            continue;
        end

        % Add balanced contributions for stimulus a
        if ca == COL_Y
            sumY  = sumY  + nEff * ma;
            sumY2 = sumY2 + nEff * m2a;
            NY = NY + nEff;
        elseif ca == COL_P
            sumP  = sumP  + nEff * ma;
            sumP2 = sumP2 + nEff * m2a;
            NP = NP + nEff;
        end

        % Add balanced contributions for stimulus b
        if cb == COL_Y
            sumY  = sumY  + nEff * mb;
            sumY2 = sumY2 + nEff * m2b;
            NY = NY + nEff;
        elseif cb == COL_P
            sumP  = sumP  + nEff * mb;
            sumP2 = sumP2 + nEff * m2b;
            NP = NP + nEff;
        end
    end

    % store gray if available
    if NG > 0
        muGray(ii,:) = (sumG / NG).';
    end

    % store pref/nonpref only if both colors have data
    if (NY < 1) || (NP < 1)
        continue;
    end

    muY = (sumY / NY).';   % 1 x nBins
    muP = (sumP / NP).';   % 1 x nBins

    if prefIsYellow(site)
        muPref(ii,:) = muY;
        muNon(ii,:)  = muP;
    else
        muPref(ii,:) = muP;
        muNon(ii,:)  = muY;
    end
end

% Drop sites that ended up NaN everywhere for pref or nonpref (gray does NOT gate selection)
good = any(isfinite(muPref), 2) & any(isfinite(muNon), 2);
muPref = muPref(good,:);
muNon  = muNon(good,:);
muGray = muGray(good,:);
fprintf('After requiring usable balanced Y/P trials: N=%d sites\n', size(muPref,1));

% -----------------------
% Plot mean ± SEM over sites
% -----------------------
mPref = mean(muPref, 1, 'omitnan');
mNon  = mean(muNon,  1, 'omitnan');
mGray = mean(muGray, 1, 'omitnan');

semPref = std(muPref, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(muPref),1));
semNon  = std(muNon,  0, 1, 'omitnan') ./ sqrt(sum(isfinite(muNon),1));
semGray = std(muGray, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(muGray),1));

figure; hold on;

% mean lines
plot(tCenters, mPref, 'LineWidth', 2);
plot(tCenters, mNon,  'LineWidth', 2);
plot(tCenters, mGray, 'LineWidth', 2);

% shaded error (simple patch)
fill([tCenters; flipud(tCenters)], [ (mPref-semPref)'; flipud((mPref+semPref)') ], ...
     'k', 'FaceAlpha', 0.12, 'EdgeColor', 'none');
fill([tCenters; flipud(tCenters)], [ (mNon-semNon)';  flipud((mNon+semNon)')  ], ...
     'k', 'FaceAlpha', 0.06, 'EdgeColor', 'none');
fill([tCenters; flipud(tCenters)], [ (mGray-semGray)'; flipud((mGray+semGray)') ], ...
     'k', 'FaceAlpha', 0.04, 'EdgeColor', 'none');

xline(0,'k-');
xlabel('Time from stimulus onset (ms)');
ylabel('Mean response (a.u.)');
title(sprintf('Balanced Pref vs Nonpref + Gray RF (abs(CI)>%.2f, N=%d)', ciThr, size(muPref,1)));
legend('Preferred','Nonpreferred','Gray (RF on background)','Location','best');
grid on;

%% Helper (inline ternary)
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
