%% Compute per-site SNR for V1 sites (1:512) for yellow vs purple, early vs late
% Tall_V1(stim).T is a 512x18 TABLE.

% -----------------------
% USER SETTINGS
% -----------------------
v1Sites = 1:512;
WIN_SPONT = 1;                % -200..0 ms
WIN_EARLY = 2;                % 40..240 ms
WIN_LATE  = 3;                % 300..500 ms

COL_YELLOW = "yellowArm";
COL_PURPLE = "purple";
COL_GRAY   = "gray";

minTrialsPerStim = 1;
minTotalTrialsPerColor = 1;
useBesselCorrection = true;

% -----------------------
% BASIC CHECKS
% -----------------------
nStim = size(R.meanAct, 2);
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

[nSitesTotal, nStim2, nWin] = size(R.meanAct);
assert(nStim2 == nStim, 'R.meanAct second dim (%d) != stimulus count (%d).', nStim2, nStim);
assert(all(size(R.meanSqAct) == size(R.meanAct)), 'R.meanSqAct must match size of R.meanAct.');
assert(nWin >= 3, 'Expected >=3 windows.');
assert(nSitesTotal >= 512, 'R.meanAct has only %d sites (<512).', nSitesTotal);
assert(numel(Tall_V1) == nStim, 'Tall_V1 has %d entries but R has %d stimuli.', numel(Tall_V1), nStim);

% Sort Tall_V1 by stimNum so TallSorted(stimNum) corresponds to stimNum
TallStimNums = arrayfun(@(x) x.stimNum, Tall_V1(:));
[sortedStimNums, order] = sort(TallStimNums(:));
assert(all(sortedStimNums(:).' == 1:nStim), 'Tall_V1.stimNum should cover 1..%d exactly.', nStim);
TallSorted = Tall_V1(order);

% Confirm T is table 512x18
T0 = TallSorted(1).T;
assert(istable(T0), 'Tall_V1(stim).T must be a table.');
assert(all(size(T0) == [512 19]), 'Expected T to be 512x18, got %s.', mat2str(size(T0)));

nSites = numel(v1Sites); % 512

% -----------------------
% PREALLOCATE
% -----------------------
muSpont = nan(nSites,1);
sdSpont = nan(nSites,1);
nSpontTotal = nan(nSites,1);

isYellow = false(nSites, nStim);
isPurple = false(nSites, nStim);

muYellowEarly = nan(nSites,1);
muYellowLate  = nan(nSites,1);
muPurpleEarly = nan(nSites,1);
muPurpleLate  = nan(nSites,1);

nYellowTrials = zeros(nSites,1);
nPurpleTrials = zeros(nSites,1);

% -----------------------
% 1) NOISE: pooled spontaneous mean and SD per V1 site across all trials/stimuli
% -----------------------
for i = 1:nSites
    s = v1Sites(i);
    if perSiteTrials
        nTrSite = double(nTrialsRaw(s,:)).';
    else
        nTrSite = nTrials;
    end
    stimOk = (nTrSite >= minTrialsPerStim);

    mu_i  = squeeze(R.meanAct(s,:,WIN_SPONT)).';    % [nStim x 1]
    msq_i = squeeze(R.meanSqAct(s,:,WIN_SPONT)).';  % [nStim x 1]

    valid = stimOk & isfinite(mu_i) & isfinite(msq_i);

    N = sum(nTrSite(valid));
    nSpontTotal(i) = N;
    if N <= 1, continue; end

    mu = sum(nTrSite(valid) .* mu_i(valid)) / N;
    muSpont(i) = mu;

    Ex2 = sum(nTrSite(valid) .* msq_i(valid)) / N;
    varPop = max(0, Ex2 - mu^2);

    if useBesselCorrection
        varEst = varPop * (N/(N-1));
    else
        varEst = varPop;
    end
    sdSpont(i) = sqrt(varEst);
end

% -----------------------
% 2) COLOR MASKS: extract center_color from TallSorted(stim).T table
% -----------------------
% Find the column that contains center_color
varNames = string(T0.Properties.VariableNames);
colIdx = find(varNames == "center_color", 1);

if isempty(colIdx)
    % Sometimes people use different naming, try a fuzzy match
    colIdx = find(contains(lower(varNames), "center") & contains(lower(varNames), "color"), 1);
end
assert(~isempty(colIdx), 'Could not find a variable named center_color (or similar) in Tall_V1(stim).T.');

for stimNum = 1:nStim
    Ttbl = TallSorted(stimNum).T;         % 512x19 table
    cc = Ttbl{:, colIdx};                 % column data (could be cellstr/categorical/string)
    labs = string(cc);                    % 512x1

    % (Optional) clean up: sometimes trailing spaces etc.
    labs = strtrim(labs);

    isYellow(:,stimNum) = (labs == COL_YELLOW);
    isPurple(:,stimNum) = (labs == COL_PURPLE);
    % gray excluded implicitly
end

% -----------------------
% 3) TRIAL-WEIGHTED MEAN RESPONSES per site/color/window
% -----------------------
for i = 1:nSites
    s = v1Sites(i);
    if perSiteTrials
        nTrSite = double(nTrialsRaw(s,:)).';
    else
        nTrSite = nTrials;
    end
    stimOk = (nTrSite >= minTrialsPerStim);

    rEarly = squeeze(R.meanAct(s,:,WIN_EARLY)).';  % [nStim x 1]
    rLate  = squeeze(R.meanAct(s,:,WIN_LATE)).';

    idxY = stimOk & isYellow(i,:).';
    NY = sum(nTrSite(idxY));
    nYellowTrials(i) = NY;
    if NY >= minTotalTrialsPerColor
        muYellowEarly(i) = sum(nTrSite(idxY) .* rEarly(idxY)) / NY;
        muYellowLate(i)  = sum(nTrSite(idxY) .* rLate(idxY))  / NY;
    end

    idxP = stimOk & isPurple(i,:).';
    NP = sum(nTrSite(idxP));
    nPurpleTrials(i) = NP;
    if NP >= minTotalTrialsPerColor
        muPurpleEarly(i) = sum(nTrSite(idxP) .* rEarly(idxP)) / NP;
        muPurpleLate(i)  = sum(nTrSite(idxP) .* rLate(idxP))  / NP;
    end
end

% -----------------------
% 4) SNR definition and outputs
% -----------------------
SNR = struct();
SNR.v1Sites = v1Sites(:);

SNR.muSpont = muSpont;
SNR.sdSpont = sdSpont;
SNR.nSpontTrials = nSpontTotal;

SNR.muYellowEarly = muYellowEarly;
SNR.muYellowLate  = muYellowLate;
SNR.muPurpleEarly = muPurpleEarly;
SNR.muPurpleLate  = muPurpleLate;

SNR.nYellowTrials = nYellowTrials;
SNR.nPurpleTrials = nPurpleTrials;

SNR.yellowEarly = (muYellowEarly - muSpont) ./ sdSpont;
SNR.yellowLate  = (muYellowLate  - muSpont) ./ sdSpont;
SNR.purpleEarly = (muPurpleEarly - muSpont) ./ sdSpont;
SNR.purpleLate  = (muPurpleLate  - muSpont) ./ sdSpont;

badNoise = ~isfinite(sdSpont) | (sdSpont <= 0);
SNR.yellowEarly(badNoise) = NaN;
SNR.yellowLate(badNoise)  = NaN;
SNR.purpleEarly(badNoise) = NaN;
SNR.purpleLate(badNoise)  = NaN;

SNR.yellowEarly(nYellowTrials < minTotalTrialsPerColor) = NaN;
SNR.yellowLate(nYellowTrials  < minTotalTrialsPerColor) = NaN;
SNR.purpleEarly(nPurpleTrials < minTotalTrialsPerColor) = NaN;
SNR.purpleLate(nPurpleTrials  < minTotalTrialsPerColor) = NaN;

% -----------------------
% 5) QUICK SUMMARIES
% -----------------------
fprintf('\n--- SNR summary (V1 sites only) ---\n');
fprintf('Using Tall_V1(stim).T column: %s\n', T0.Properties.VariableNames{colIdx});
fprintf('Valid noise SD: %d / %d\n', sum(isfinite(sdSpont) & sdSpont>0), nSites);
fprintf('Median spont SD: %.4f\n', median(sdSpont(isfinite(sdSpont) & sdSpont>0)));
fprintf('Sites with >=%d yellow pooled trials: %d\n', minTotalTrialsPerColor, sum(nYellowTrials>=minTotalTrialsPerColor));
fprintf('Sites with >=%d purple pooled trials: %d\n', minTotalTrialsPerColor, sum(nPurpleTrials>=minTotalTrialsPerColor));

% quick look at strongest yellowEarly
[~,ordY] = sort(SNR.yellowEarly, 'descend', 'MissingPlacement','last');
kshow = min(10, nSites);
disp(table((1:kshow).', ordY(1:kshow), v1Sites(ordY(1:kshow)).', SNR.yellowEarly(ordY(1:kshow)), nYellowTrials(ordY(1:kshow)), ...
    'VariableNames', {'Rank','IdxInV1Vector','SiteNumberInR','SNR_yellowEarly','nYellowTrials'}));

% save('SNR_V1_byColor_byWindow.mat','SNR');
