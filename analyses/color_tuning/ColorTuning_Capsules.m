%% Color tuning with complementary balancing (V1, 512 sites)

% Here is the prompt:

% Now it is time to compute the color tuning for those sites with a sufficient SNR. For this we are going to look into the second and third time window, 
% i.e. 40-240ms and 300-500ms. We can use the fact that of the 384 stimuli the V1 receptive fields (RFs) fall either on the yellow capsule, on the purple capsule,
% or on the grey background.
% We are going to disregard responses to the grey background. 
% There is one important aspect in this stimulus set: the 384 stimuli come in groups of 8. Stimulus 1, 2, 5 and 6 have the same global configuration and the 
% same is true for stimulus 3, 4, 7 and 8. Here stimulus 1 (3) is the same as 5 (7), except that all pixels that are yellow in 1 (3) are purple in 5 (7) 
% and vice versa. The same is true for stimulus 2 (4) and 6 (8). I therefore call 1 and 5 complementary, and the same is true for 2-6, 3-7 and 4-8. 
% Complementary stimuli look the same the background is identical but the yellow and purple switch. This pattern repeats for stimuli 9-16, 17-24, etc., etc. 

% The number of trials per stimulus is relatively small, less than 10. But there are many stimuli (out of the 384) for which the RF of a recording sites 
% falls on a non-background pixel that is yellow or purple. 
% To avoid measuring tuning to other aspects than color, we will need to include the same number of trials for complementary stimuli. 
% For example if I have 4 trials for stimulus 1 and 6 trials for stimulus 5, I can maximally select 4 trials for both complementary stimuli. In this example 
% I would have to select randomly 4 of the 6 trials with stimulus 5.  
% Once we have a balanced stimulus set we can use all the remaining trials to compute:
% 1)	Colorindex, which is (yellow-purple)/modulus(average)
% 2)	D-prime, which can be based on the fact that we also stored the average squared response for every stimulus
% 3)	Significance of the difference, e.g. based on z-score

% Requires:
%   R.meanAct   [>=512 x 384 x 3]
%   R.meanSqAct [>=512 x 384 x 3]
%   R.nTrials   [384 x 1]
%   Tall_V1(384x1).T is a 512x18 table with variable center_color (as before)
%   SNR.bestSNR [512 x 1] from previous step

% -----------------------
% Settings
% -----------------------
keepSites = find(bestSNR > SNRthr);   % indices 1..512
fprintf('Keeping %d / 512 sites (bestSNR > %.2f)\n', numel(keepSites), SNRthr);

WIN_EARLY = 2; % 40-240 ms
WIN_LATE  = 3; % 300-500 ms

COL_Y = "yellowArm";
COL_P = "purple";
COL_G = "gray";

useBesselCorrection = true; % variance within stimulus: n/(n-1) correction
alpha = 0.05;              % significance threshold for paired histogram overlay

cfg = config();
legacyColorTunePath = fullfile(cfg.matDir, 'ColorTune_balanced_V1.mat');
loadedLegacyColorTune = false;

nStim = size(R.meanAct, 2);
assert(nStim == 384, 'Expected 384 stimuli.');

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

% -----------------------
% Sort Tall_V1 by stimNum (safety)
% -----------------------
TallStimNums = arrayfun(@(x) x.stimNum, Tall_V1(:));
[sortedStimNums, order] = sort(TallStimNums(:));
assert(all(sortedStimNums(:).' == 1:nStim), 'Tall_V1.stimNum should cover 1..%d.', nStim);
TallSorted = Tall_V1(order);

% -----------------------
% Extract center_color column index
% -----------------------
T0 = TallSorted(1).T;
assert(istable(T0), 'Tall_V1(stim).T must be a table.');
varNames = string(T0.Properties.VariableNames);
colIdx = find(varNames == "center_color", 1);
if isempty(colIdx)
    colIdx = find(contains(lower(varNames), "center") & contains(lower(varNames), "color"), 1);
end
assert(~isempty(colIdx), 'Could not find center_color column in Tall_V1(stim).T.');

% Build label matrix CC: [512 x 384]
CC = strings(512, nStim);
for stim = 1:nStim
    Ttbl = TallSorted(stim).T;
    labs = string(Ttbl{:, colIdx});  % 512x1
    CC(:,stim) = strtrim(labs);
end

% -----------------------
% Complement pairing within each block of 8
% Pair map: (1<->5), (2<->6), (3<->7), (4<->8), then repeats
% We'll iterate stimA = 1..384 and only handle pairs where stimA is the "first" (within-block 1..4).
% -----------------------
pairsA = []; pairsB = [];
for stimA = 1:nStim
    pos = mod(stimA-1, 8) + 1; % 1..8 within block
    if pos <= 4
        stimB = stimA + 4;
        pairsA(end+1,1) = stimA; %#ok<AGROW>
        pairsB(end+1,1) = stimB; %#ok<AGROW>
    end
end
nPairs = numel(pairsA);
assert(nPairs == 192, 'Expected 192 complementary pairs.');

% -----------------------
% Helper to accumulate balanced sums for a given site & window
% We accumulate for Yellow and Purple separately:
%   sumX  = Σ nEff * mean
%   sumX2 = Σ nEff * meanSq
%   N     = Σ nEff
% -----------------------
nSites = 512;
v1RowsInR = (1:512);

% Output arrays (512x1), set NaN for non-kept sites
if exist(legacyColorTunePath, 'file') == 2
    Slegacy = load(legacyColorTunePath, 'ColorTune');
    if isfield(Slegacy, 'ColorTune') && isstruct(Slegacy.ColorTune) && ...
            isfield(Slegacy.ColorTune, 'early') && isfield(Slegacy.ColorTune, 'late') && ...
            isfield(Slegacy.ColorTune.early, 'colorIndex') && isfield(Slegacy.ColorTune.late, 'colorIndex')
        ColorTune = Slegacy.ColorTune;
        loadedLegacyColorTune = true;
        fprintf('Loaded historical V1 color tuning from %s\n', legacyColorTunePath);
    end
end

if ~loadedLegacyColorTune
    ColorTune = struct();
    ColorTune.keepSites = keepSites;
    ColorTune.thr = SNRthr;

    % For each window store mu/var/N for Y and P
    fields = ["muY","muP","varY","varP","NY","NP","colorIndex","dprime","z","p"];
    for f = fields
        ColorTune.early.(f) = nan(nSites,1);
        ColorTune.late.(f)  = nan(nSites,1);
    end

    % -----------------------
    % Main loop over sites (only compute kept ones)
    % -----------------------
    for iSite = keepSites(:).'
        s = v1RowsInR(iSite); % row in R

        if perSiteTrials
            nTrSite = double(nTrialsRaw(s,:)).';
        else
            nTrSite = nTrials;
        end

        % Pre-fetch responses for speed
        mEarly  = squeeze(R.meanAct(s,:,WIN_EARLY)).';   % [384 x 1]
        msqEarly= squeeze(R.meanSqAct(s,:,WIN_EARLY)).'; % [384 x 1]
        mLate   = squeeze(R.meanAct(s,:,WIN_LATE)).';
        msqLate = squeeze(R.meanSqAct(s,:,WIN_LATE)).';

        % Accumulators for early
        sumY_e = 0; sumY2_e = 0; NY_e = 0;
        sumP_e = 0; sumP2_e = 0; NP_e = 0;

        % Accumulators for late
        sumY_l = 0; sumY2_l = 0; NY_l = 0;
        sumP_l = 0; sumP2_l = 0; NP_l = 0;

        for ip = 1:nPairs
            a = pairsA(ip);
            b = pairsB(ip);

            na = nTrSite(a);
            nb = nTrSite(b);
            nEff = min(na, nb);

            if nEff <= 0
                continue;
            end

            % Exclude background: if RF on gray in BOTH, skip.
            % (In your design, if it's gray in one it should be gray in complement too.)
            ca = CC(iSite,a);
            cb = CC(iSite,b);

            if (ca == COL_G) && (cb == COL_G)
                continue;
            end

            % Safety: if either is missing/empty, skip
            if (ca == "" || ismissing(ca) || cb == "" || ismissing(cb))
                continue;
            end

            % Early window: stimulus a contributes to color ca, stimulus b contributes to color cb
            % Add nEff "balanced" trials from each
            if ca == COL_Y
                sumY_e  = sumY_e  + nEff * mEarly(a);
                sumY2_e = sumY2_e + nEff * msqEarly(a);
                NY_e    = NY_e    + nEff;
            elseif ca == COL_P
                sumP_e  = sumP_e  + nEff * mEarly(a);
                sumP2_e = sumP2_e + nEff * msqEarly(a);
                NP_e    = NP_e    + nEff;
            end

            if cb == COL_Y
                sumY_e  = sumY_e  + nEff * mEarly(b);
                sumY2_e = sumY2_e + nEff * msqEarly(b);
                NY_e    = NY_e    + nEff;
            elseif cb == COL_P
                sumP_e  = sumP_e  + nEff * mEarly(b);
                sumP2_e = sumP2_e + nEff * msqEarly(b);
                NP_e    = NP_e    + nEff;
            end

            % Late window
            if ca == COL_Y
                sumY_l  = sumY_l  + nEff * mLate(a);
                sumY2_l = sumY2_l + nEff * msqLate(a);
                NY_l    = NY_l    + nEff;
            elseif ca == COL_P
                sumP_l  = sumP_l  + nEff * mLate(a);
                sumP2_l = sumP2_l + nEff * msqLate(a);
                NP_l    = NP_l    + nEff;
            end

            if cb == COL_Y
                sumY_l  = sumY_l  + nEff * mLate(b);
                sumY2_l = sumY2_l + nEff * msqLate(b);
                NY_l    = NY_l    + nEff;
            elseif cb == COL_P
                sumP_l  = sumP_l  + nEff * mLate(b);
                sumP2_l = sumP2_l + nEff * msqLate(b);
                NP_l    = NP_l    + nEff;
            end
        end

        % Compute stats for early + late
        [muY, varY, muP, varP, NY, NP, cind, dp, z, p] = compute_metrics(sumY_e,sumY2_e,NY_e,sumP_e,sumP2_e,NP_e,useBesselCorrection);
        ColorTune.early.muY(iSite) = muY;  ColorTune.early.varY(iSite) = varY;  ColorTune.early.NY(iSite) = NY;
        ColorTune.early.muP(iSite) = muP;  ColorTune.early.varP(iSite) = varP;  ColorTune.early.NP(iSite) = NP;
        ColorTune.early.colorIndex(iSite) = cind;
        ColorTune.early.dprime(iSite)     = dp;
        ColorTune.early.z(iSite)          = z;
        ColorTune.early.p(iSite)          = p;

        [muY, varY, muP, varP, NY, NP, cind, dp, z, p] = compute_metrics(sumY_l,sumY2_l,NY_l,sumP_l,sumP2_l,NP_l,useBesselCorrection);
        ColorTune.late.muY(iSite) = muY;   ColorTune.late.varY(iSite) = varY;   ColorTune.late.NY(iSite) = NY;
        ColorTune.late.muP(iSite) = muP;   ColorTune.late.varP(iSite) = varP;   ColorTune.late.NP(iSite) = NP;
        ColorTune.late.colorIndex(iSite) = cind;
        ColorTune.late.dprime(iSite)     = dp;
        ColorTune.late.z(iSite)          = z;
        ColorTune.late.p(iSite)          = p;
    end

    fprintf('Done. Example: median d'' early (kept) = %.3f\n', median(ColorTune.early.dprime(keepSites), 'omitnan'));
end

if ~isfield(ColorTune, 'keepSites') || isempty(ColorTune.keepSites)
    ColorTune.keepSites = keepSites;
end
if ~isfield(ColorTune, 'thr') || isempty(ColorTune.thr)
    ColorTune.thr = SNRthr;
end
if loadedLegacyColorTune && ~isequal(ColorTune.keepSites(:), keepSites(:))
    fprintf(['Historical V1 color file uses %d keep-sites; current workspace has %d. ' ...
        'Plotting will follow the historical keep-sites so the color-index histogram stays stable.\n'], ...
        numel(ColorTune.keepSites), numel(keepSites));
end

plotKeepSites = ColorTune.keepSites(:);

% Add paired precision-weighted significance using the shared helper while
% preserving the historical pooled colorIndex/dprime/z outputs above.
siteRows = (1:nSites).';
ColorTunePair = compute_color_tuning_balanced_sites(R, Tall_V1, siteRows, plotKeepSites, ...
    'WIN_EARLY', WIN_EARLY, 'WIN_LATE', WIN_LATE, ...
    'useBesselCorrection', useBesselCorrection, 'Verbose', true);
pairedFields = ["pairedMuY","pairedMuP","pairedMeanDiff","pairedWeightedDiff", ...
    "pairedWeightSum","pairedEffN","pairedT","pairedP","pairedNPairs"];
for f = pairedFields
    ColorTune.early.(f) = ColorTunePair.early.(f);
    ColorTune.late.(f) = ColorTunePair.late.(f);
end

ciEarlyAll = ColorTune.early.colorIndex(plotKeepSites);
pEarlyAll = ColorTune.early.p(plotKeepSites);
keepEarly = isfinite(ciEarlyAll);
ciEarly = ciEarlyAll(keepEarly);
sigEarlyMask = keepEarly & isfinite(pEarlyAll) & (pEarlyAll < alpha);
ciEarlySig = ciEarlyAll(sigEarlyMask);

ciLateAll = ColorTune.late.colorIndex(plotKeepSites);
pLateAll = ColorTune.late.p(plotKeepSites);
keepLate = isfinite(ciLateAll);
ciLate = ciLateAll(keepLate);
sigLateMask = keepLate & isfinite(pLateAll) & (pLateAll < alpha);
ciLateSig = ciLateAll(sigLateMask);

fprintf('Significant V1 color tuning by pooled test (p < %.2f): early=%d late=%d\n', ...
    alpha, nnz(sigEarlyMask), nnz(sigLateMask));

figure('Color','w');
useTiled = exist('tiledlayout', 'file') == 2;
if useTiled
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
end

if useTiled, nexttile; else, subplot(1, 2, 1); end
hEarly = histogram(ciEarly, 30, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
hold on;
histogram(ciEarlySig, 'BinEdges', hEarly.BinEdges, 'FaceColor', [0.8 0 0], 'EdgeColor', 'none');
xline(0, 'k-');
xlabel('Color Index (yellow - purple)');
ylabel('Number of sites');
title('Early');
legend('All sites', sprintf('Significant (p<%.2f)', alpha));
grid on;

if useTiled, nexttile; else, subplot(1, 2, 2); end
hLate = histogram(ciLate, 30, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
hold on;
histogram(ciLateSig, 'BinEdges', hLate.BinEdges, 'FaceColor', [0.8 0 0], 'EdgeColor', 'none');
xline(0, 'k-');
xlabel('Color Index (yellow - purple)');
ylabel('Number of sites');
title('Late');
legend('All sites', sprintf('Significant (p<%.2f)', alpha));
grid on;

if exist('sgtitle', 'file') == 2
    sgtitle(sprintf('V1 Color tuning, SNR>%.2f', SNRthr));
end

% Optional save
% save('ColorTune_balanced_V1.mat','ColorTune');

%% ---- local helper function ----
function [muY,varY,muP,varP,NY,NP,colorIndex,dprime,z,p] = compute_metrics(sumY,sumY2,NY,sumP,sumP2,NP,useBessel)
    % Defaults
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
        varY = varY * (NY/(NY-1));
        varP = varP * (NP/(NP-1));
    end

    % 1) Color index: (Y - P) / |(Y + P)/2|
    denom = abs((muY + muP)/2);
    if denom > 0
        colorIndex = (muY - muP) / denom;
    end

    % 2) d-prime: (Y - P) / sqrt(0.5*(varY + varP))
    denomDP = sqrt(0.5*(varY + varP));
    if denomDP > 0
        dprime = (muY - muP) / denomDP;
    end

    % 3) z-score for difference of means (Welch-style)
    denomZ = sqrt(varY/NY + varP/NP);
    if denomZ > 0
        z = (muY - muP) / denomZ;
        p = 2 * normcdf(-abs(z), 0, 1);
    end
end
