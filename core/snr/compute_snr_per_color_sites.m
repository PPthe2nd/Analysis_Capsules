function SNR = compute_snr_per_color_sites(DATA, Tall, siteRows, varargin)
% COMPUTE_SNR_PER_COLOR_SITES
% Compute per-site color-specific SNR for an arbitrary set of site rows.
%
% DATA must contain:
%   DATA.meanAct    [nSites x nStim x >=3]
%   DATA.meanSqAct  [nSites x nStim x >=3]
%   DATA.nTrials    [1 x nStim] or [nSites x nStim]
%
% Tall is a struct array with Tall(stim).T tables containing center_color.
%
% siteRows are row indices into DATA.meanAct and the Tall(stim).T tables.
%
% Output fields match the SNR struct expected by attention_modulation_V1_3bin.

p = inputParser;
p.addParameter('WIN_SPONT', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('WIN_EARLY', 2, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('WIN_LATE',  3, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('minTrialsPerStim', 1, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('minTotalTrialsPerColor', 1, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('useBesselCorrection', true, @(x) islogical(x) && isscalar(x));
p.addParameter('Verbose', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;

assert(isfield(DATA, 'meanAct') && isfield(DATA, 'meanSqAct') && isfield(DATA, 'nTrials'), ...
    'DATA must contain meanAct, meanSqAct, and nTrials.');

[nSitesTotal, nStim, nWin] = size(DATA.meanAct);
assert(all(size(DATA.meanSqAct) == size(DATA.meanAct)), ...
    'DATA.meanSqAct must match DATA.meanAct.');
assert(nWin >= max([opt.WIN_SPONT, opt.WIN_EARLY, opt.WIN_LATE]), ...
    'DATA does not contain the requested time windows.');
assert(numel(Tall) == nStim, 'Tall has %d entries but DATA has %d stimuli.', numel(Tall), nStim);

siteRows = siteRows(:);
assert(~isempty(siteRows), 'siteRows must not be empty.');
assert(all(siteRows >= 1 & siteRows <= nSitesTotal), ...
    'siteRows must lie within DATA.meanAct rows 1..%d.', nSitesTotal);

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

COL_YELLOW = "yellowArm";
COL_PURPLE = "purple";

muSpont = nan(nSitesTotal,1);
sdSpont = nan(nSitesTotal,1);
nSpontTotal = nan(nSitesTotal,1);

muYellowEarly = nan(nSitesTotal,1);
muYellowLate  = nan(nSitesTotal,1);
muPurpleEarly = nan(nSitesTotal,1);
muPurpleLate  = nan(nSitesTotal,1);

nYellowTrials = zeros(nSitesTotal,1);
nPurpleTrials = zeros(nSitesTotal,1);

nSel = numel(siteRows);
isYellow = false(nSel, nStim);
isPurple = false(nSel, nStim);

for stimNum = 1:nStim
    Ttbl = TallSorted(stimNum).T;
    labs = string(Ttbl{siteRows, colIdx});
    labs = strtrim(labs);
    isYellow(:, stimNum) = (labs == COL_YELLOW);
    isPurple(:, stimNum) = (labs == COL_PURPLE);
end

for i = 1:nSel
    s = siteRows(i);

    if perSiteTrials
        nTr = nTrials(s, :).';
    else
        nTr = nTrials;
    end

    mu_i  = squeeze(DATA.meanAct(s,:,opt.WIN_SPONT));
    msq_i = squeeze(DATA.meanSqAct(s,:,opt.WIN_SPONT));
    mu_i = mu_i(:);
    msq_i = msq_i(:);

    valid = isfinite(mu_i) & isfinite(msq_i) & isfinite(nTr) & (nTr >= opt.minTrialsPerStim);
    N = sum(nTr(valid));
    nSpontTotal(s) = N;
    if N > 1
        mu = sum(nTr(valid) .* mu_i(valid)) / N;
        muSpont(s) = mu;

        Ex2 = sum(nTr(valid) .* msq_i(valid)) / N;
        varPop = max(0, Ex2 - mu.^2);
        if opt.useBesselCorrection
            varEst = varPop * (N / (N - 1));
        else
            varEst = varPop;
        end
        sdSpont(s) = sqrt(varEst);
    end

    rEarly = squeeze(DATA.meanAct(s,:,opt.WIN_EARLY));
    rLate  = squeeze(DATA.meanAct(s,:,opt.WIN_LATE));
    rEarly = rEarly(:);
    rLate = rLate(:);

    idxY = isYellow(i,:).' & isfinite(rEarly) & isfinite(rLate) & isfinite(nTr) & (nTr >= opt.minTrialsPerStim);
    NY = sum(nTr(idxY));
    nYellowTrials(s) = NY;
    if NY >= opt.minTotalTrialsPerColor
        muYellowEarly(s) = sum(nTr(idxY) .* rEarly(idxY)) / NY;
        muYellowLate(s)  = sum(nTr(idxY) .* rLate(idxY))  / NY;
    end

    idxP = isPurple(i,:).' & isfinite(rEarly) & isfinite(rLate) & isfinite(nTr) & (nTr >= opt.minTrialsPerStim);
    NP = sum(nTr(idxP));
    nPurpleTrials(s) = NP;
    if NP >= opt.minTotalTrialsPerColor
        muPurpleEarly(s) = sum(nTr(idxP) .* rEarly(idxP)) / NP;
        muPurpleLate(s)  = sum(nTr(idxP) .* rLate(idxP))  / NP;
    end
end

SNR = struct();
SNR.siteRows = siteRows;
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

tooFewY = nYellowTrials < opt.minTotalTrialsPerColor;
tooFewP = nPurpleTrials < opt.minTotalTrialsPerColor;
SNR.yellowEarly(tooFewY) = NaN;
SNR.yellowLate(tooFewY)  = NaN;
SNR.purpleEarly(tooFewP) = NaN;
SNR.purpleLate(tooFewP)  = NaN;

if opt.Verbose
    goodNoise = isfinite(sdSpont(siteRows)) & sdSpont(siteRows) > 0;
    fprintf('\n--- SNR summary (selected site rows) ---\n');
    fprintf('Using Tall(stim).T column: %s\n', T0.Properties.VariableNames{colIdx});
    fprintf('Valid noise SD: %d / %d\n', sum(goodNoise), nSel);
    if any(goodNoise)
        fprintf('Median spont SD: %.4f\n', median(sdSpont(siteRows(goodNoise))));
    end
    fprintf('Sites with >=%d yellow pooled trials: %d\n', ...
        opt.minTotalTrialsPerColor, sum(nYellowTrials(siteRows) >= opt.minTotalTrialsPerColor));
    fprintf('Sites with >=%d purple pooled trials: %d\n', ...
        opt.minTotalTrialsPerColor, sum(nPurpleTrials(siteRows) >= opt.minTotalTrialsPerColor));
end

end
