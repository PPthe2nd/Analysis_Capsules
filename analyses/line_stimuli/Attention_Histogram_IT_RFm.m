function OUT = Attention_Histogram_IT_RFm(Monkey, Opts)
% ATTENTION_HISTOGRAM_IT_RFM
% Histogram of IT attention indices in the line task using RF.m-center
% Tall_IT geometry, gated by attention-validity, object assignment, and
% early quartet responsiveness.

if nargin < 1 || isempty(Monkey)
    Monkey = 1; % 1 = Nilson, 2 = Figaro
end
if nargin < 2 || isempty(Opts)
    Opts = struct();
end

Opts = normalize_opts_local(Opts);
cfg = config();

if Monkey == 1
    monkeySuffix = "N";
    tallFile = 'Tall_IT_lines_N.mat';
    resp3binFile = 'SNR_capsules_N_d12.mat';
elseif Monkey == 2
    monkeySuffix = "F";
    tallFile = 'Tall_IT_lines_F.mat';
    resp3binFile = 'SNR_capsules_F_d12.mat';
else
    error('Attention_Histogram_IT_RFm:InvalidMonkey', ...
        'Monkey must be 1 (Nilson) or 2 (Figaro).');
end

tallPath = fullfile(cfg.matDir, tallFile);
resp3binPath = fullfile(cfg.matDir, resp3binFile);

assert(exist(tallPath, 'file') == 2, ...
    'Missing %s. Run Line_Stimuli_IT.m first.', tallPath);
assert(exist(resp3binPath, 'file') == 2, ...
    'Missing %s. Create the 3-bin response summary first.', resp3binPath);

Sgeo = load(tallPath);
assert(isfield(Sgeo, 'Tall_IT') && isstruct(Sgeo.Tall_IT), ...
    '%s must contain struct Tall_IT.', tallFile);
assert(isfield(Sgeo, 'RFrange') && ~isempty(Sgeo.RFrange), ...
    '%s must contain non-empty RFrange.', tallFile);

Tall_IT = Sgeo.Tall_IT;
RFrange = Sgeo.RFrange(:);
nIT = numel(RFrange);
siteRows = (1:nIT).';

stimNums = arrayfun(@(x) x.stimNum, Tall_IT(:));
[stimNumsSorted, ordStim] = sort(stimNums(:));
assert(all(stimNumsSorted(:).' == 1:numel(Tall_IT)), ...
    'Tall_IT.stimNum must cover 1..%d exactly.', numel(Tall_IT));
Tall_IT = Tall_IT(ordStim);
nStim = numel(Tall_IT);

Sresp = load(resp3binPath);
assert(isfield(Sresp, 'R') && isstruct(Sresp.R), ...
    '%s must contain struct R.', resp3binFile);

R3_full = Sresp.R;
R3 = localize_response_rows_local(R3_full, RFrange);
assert(size(R3.meanAct,1) == nIT, 'Localized IT response rows do not match Tall_IT.');
assert(size(R3.meanAct,2) == nStim, 'Localized IT response stimuli do not match Tall_IT.');

SNRnorm = compute_snr_per_color_sites(R3, Tall_IT, siteRows, 'Verbose', false);
optsTD = struct('v1Sites', siteRows, 'timeIdx', Opts.timeIdx3, ...
    'excludeOverlap', Opts.excludeOverlap, 'verbose', false, 'epsDen', Opts.epsDen);
ATT = attention_modulation_V1_3bin(R3, Tall_IT, SNRnorm, optsTD);

denCheck = abs((ATT.muT + ATT.muD) / 2);
denCheck(denCheck < Opts.epsDen) = Opts.epsDen;
idxCheck = (ATT.muT - ATT.muD) ./ denCheck;
okCheck = isfinite(ATT.index) & isfinite(idxCheck);
if any(okCheck)
    maxErr = max(abs(ATT.index(okCheck) - idxCheck(okCheck)));
    assert(maxErr < 1e-10, ...
        'Attention index mismatch: max abs difference %.3g.', maxErr);
end

[hasObjectRF, nObjectStim, topKAbsQuartetEarly] = compute_it_responsive_gate_local( ...
    Tall_IT, R3, SNRnorm, Opts.MinObjectStim, Opts.TopKQuartets);

keepValid = ATT.validSite(:);
keepObject = hasObjectRF(:);
keepResponsive = isfinite(topKAbsQuartetEarly) & (topKAbsQuartetEarly > Opts.RespThr);
keep = keepValid & keepObject & keepResponsive;

isSig = keep & isfinite(ATT.pValueTD(:)) & (ATT.pValueTD(:) < Opts.pThresh);
Nmatch = ATT.wY(:) + ATT.wP(:);

idxAll = ATT.index(:);
vAllRaw = idxAll(keep & isfinite(idxAll));
vSigRaw = idxAll(isSig & isfinite(idxAll));
vAll = max(min(vAllRaw, Opts.idxClip), -Opts.idxClip);
vSig = max(min(vSigRaw, Opts.idxClip), -Opts.idxClip);

fprintf('IT attention histogram RF.m (%s)\n', char(monkeySuffix));
fprintf('Attention-valid sites: %d / %d\n', nnz(keepValid), nIT);
fprintf('With object RF (>= %d object stim): %d / %d\n', Opts.MinObjectStim, nnz(keepObject), nIT);
fprintf('Early quartet-responsive (top-%d abs > %.2f): %d / %d\n', ...
    Opts.TopKQuartets, Opts.RespThr, nnz(keepResponsive), nIT);
fprintf('Kept for histogram (valid + object RF + early responsive): %d / %d\n', nnz(keep), nIT);
fprintf('Attention-significant within kept (pTD < %.3f): %d / %d kept\n', ...
    Opts.pThresh, nnz(isSig), nnz(keep));
if any(keep)
    fprintf('Median matched trial count Nmatch among kept sites: %.1f\n', median(Nmatch(keep), 'omitnan'));
    fprintf('Median early top-%d quartet response among kept sites: %.3f\n', ...
        Opts.TopKQuartets, median(topKAbsQuartetEarly(keep), 'omitnan'));
end
fprintf('Display clipped at +/-%.1f only: all=%d, significant=%d values clipped\n', ...
    Opts.idxClip, nnz(abs(vAllRaw) > Opts.idxClip), nnz(abs(vSigRaw) > Opts.idxClip));

fig = [];
ax = [];
if Opts.PlotFigure
    fig = figure('Color', 'w', 'Name', sprintf('IT attention index RF.m (%s)', char(monkeySuffix)), ...
        'NumberTitle', 'off', 'Tag', 'IT_attention_hist_RFm');
    ax = axes('Parent', fig);
    hold(ax, 'on');
    histogram(ax, vAll, Opts.nHistBins, ...
        'FaceColor', [0.80 0.80 0.80], 'EdgeColor', 'none');
    histogram(ax, vSig, Opts.nHistBins, ...
        'FaceColor', [0.85 0.20 0.20], 'EdgeColor', 'none');
    xlabel(ax, 'Attention index');
    ylabel(ax, 'Number of sites');
    title(ax, sprintf(['IT Attention index RF.m (%s, valid + objectRF + early responsive, ' ...
        'pTD < %.3f)'], char(monkeySuffix), Opts.pThresh));
    legend(ax, ...
        sprintf('All kept sites (N=%d)', numel(vAll)), ...
        sprintf('Significant sites (N=%d)', numel(vSig)), ...
        'Location', 'best');
    grid(ax, 'on');
end

OUT = struct();
OUT.Monkey = Monkey;
OUT.monkeySuffix = monkeySuffix;
OUT.tallPath = tallPath;
OUT.resp3binPath = resp3binPath;
OUT.filters = Opts;
OUT.ATT = ATT;
OUT.keepValid = keepValid;
OUT.keepObject = keepObject;
OUT.keepResponsive = keepResponsive;
OUT.keep = keep;
OUT.isSig = isSig;
OUT.Nmatch = Nmatch;
OUT.nObjectStim = nObjectStim;
OUT.topKAbsQuartetEarly = topKAbsQuartetEarly;
OUT.idxAllRaw = vAllRaw;
OUT.idxSigRaw = vSigRaw;
OUT.idxClip = Opts.idxClip;
OUT.figure = fig;
OUT.axes = ax;
end

function Opts = normalize_opts_local(Opts)
defaults = struct();
defaults.timeIdx3 = 3;
defaults.pThresh = 0.05;
defaults.idxClip = 2;
defaults.excludeOverlap = true;
defaults.nHistBins = 30;
defaults.epsDen = 1e-6;
defaults.MinObjectStim = 1;
defaults.RespThr = 0.7;
defaults.TopKQuartets = 5;
defaults.PlotFigure = true;

fn = fieldnames(defaults);
for i = 1:numel(fn)
    if ~isfield(Opts, fn{i}) || isempty(Opts.(fn{i}))
        Opts.(fn{i}) = defaults.(fn{i});
    end
end
end

function R_loc = localize_response_rows_local(R_full, siteGlobal)
R_loc = R_full;
R_loc.meanAct = R_full.meanAct(siteGlobal, :, :);
R_loc.meanSqAct = R_full.meanSqAct(siteGlobal, :, :);
if ismatrix(R_full.nTrials) && size(R_full.nTrials,1) >= max(siteGlobal)
    R_loc.nTrials = R_full.nTrials(siteGlobal, :);
else
    R_loc.nTrials = R_full.nTrials;
end
end

function [hasObjectRF, nObjectStim, topKAbsQuartetEarly] = compute_it_responsive_gate_local( ...
        Tall_IT, R3, SNR, MinObjectStim, TopKQuartets)
nIT = size(R3.meanAct, 1);
nStim = numel(Tall_IT);
siteRows = (1:nIT).';

nTargetStim = zeros(nIT, 1);
nDistrStim = zeros(nIT, 1);
for stim = 1:nStim
    T = Tall_IT(stim).T;
    assign = string(T.assignment(siteRows));
    nTargetStim = nTargetStim + (assign == "target");
    nDistrStim = nDistrStim + (assign == "distractor");
end
nObjectStim = nTargetStim + nDistrStim;
hasObjectRF = nObjectStim >= MinObjectStim;

if isvector(R3.nTrials)
    nTrialsAll = double(R3.nTrials(:));
    perSiteTrials = false;
else
    perSiteTrials = true;
end

muSpont = SNR.muSpont(siteRows);
sdSpont = SNR.sdSpont(siteRows);
badNoise = ~isfinite(sdSpont) | (sdSpont <= 0);

nQuartets = floor(nStim / 8) * 2;
quartetMembers = zeros(nQuartets, 4);
q = 0;
for base = 0:8:(nStim - 8)
    q = q + 1;
    quartetMembers(q,:) = base + [1 2 5 6];
    q = q + 1;
    quartetMembers(q,:) = base + [3 4 7 8];
end

muQuartetEarly = nan(nIT, nQuartets);
for iSite = 1:nIT
    if perSiteTrials
        nTrSite = double(R3.nTrials(iSite, :)).';
    else
        nTrSite = nTrialsAll;
    end

    rEarly = squeeze(R3.meanAct(iSite, :, 2)).';
    for qIdx = 1:nQuartets
        stimQ = quartetMembers(qIdx, :);
        nTrQ = nTrSite(stimQ);
        rEarlyQ = rEarly(stimQ);
        idxEarly = isfinite(rEarlyQ) & isfinite(nTrQ) & (nTrQ > 0);
        if any(idxEarly)
            muQuartetEarly(iSite, qIdx) = sum(nTrQ(idxEarly) .* rEarlyQ(idxEarly)) / sum(nTrQ(idxEarly));
        end
    end
end

signedQuartetEarly = bsxfun(@rdivide, bsxfun(@minus, muQuartetEarly, muSpont), sdSpont);
signedQuartetEarly(badNoise, :) = NaN;

topKAbsQuartetEarly = nan(nIT, 1);
for iSite = 1:nIT
    vals = abs(signedQuartetEarly(iSite, :));
    vals = vals(isfinite(vals));
    if isempty(vals)
        continue;
    end
    vals = sort(vals, 'descend');
    kUse = min(TopKQuartets, numel(vals));
    topKAbsQuartetEarly(iSite) = mean(vals(1:kUse));
end
end
