% First-pass V4 analysis for the line-task stimulus set.
%
% This script intentionally keeps the scope narrow:
% - geometry QC from Tall_V4
% - onset-driven screening via color-specific SNR
% - attention-modulation screening in the 3-bin data
% - basic plots and site-count summaries
%
% Deferred for later because they are riskier or more V1-specific:
% - color-tuning rescue beyond the current balanced ColorTune analysis in V4
% - replacing center-based RF assignment with overlap-based inclusion
% - refactoring V1 and V4 scripts into one shared pipeline

Monkey = 1; % 1 = Nilson, 2 = Figaro
ExampleStimulus = 38; % stimulus used for quick overlay/QC plots
SNRthr = 0.7; % minimum onset SNR for the primary visually driven site set
pTDthr = 0.05; % p-value threshold for attention-modulation rescue
NminMatched = 20; % minimum matched T/D trial weight after color balancing for attention rescue
MinObjectStim = 1; % site must fall on target/distractor in at least this many stimuli

cfg = config();

if Monkey == 1
    monkeySuffix = "N";
    tallFile = 'Tall_V4_lines_N.mat';
    resp3binFile = 'SNR_capsules_N_d12.mat';
elseif Monkey == 2
    monkeySuffix = "F";
    tallFile = 'Tall_V4_lines_F.mat';
    resp3binFile = 'SNR_capsules_F_d12.mat';
else
    error('Analyse_Line_Stimuli_V4:InvalidMonkey', ...
        'Monkey must be 1 (Nilson) or 2 (Figaro).');
end

tallPath = fullfile(cfg.matDir, tallFile);
resp3binPath = fullfile(cfg.matDir, resp3binFile);

assert(exist(tallPath, 'file') == 2, ...
    'Missing %s. Run Line_Stimuli_V4.m first.', tallPath);
assert(exist(resp3binPath, 'file') == 2, ...
    'Missing %s. Create the 3-bin response summary for this monkey first.', resp3binPath);

Sgeo = load(tallPath);
assert(isfield(Sgeo, 'Tall_V4') && isstruct(Sgeo.Tall_V4), ...
    '%s must contain struct Tall_V4.', tallFile);
assert(isfield(Sgeo, 'ALLCOORDS') && isfield(Sgeo, 'RTAB384'), ...
    '%s must contain ALLCOORDS and RTAB384.', tallFile);
assert(isfield(Sgeo, 'RFrange') && ~isempty(Sgeo.RFrange), ...
    '%s must contain non-empty RFrange.', tallFile);

Tall_V4 = Sgeo.Tall_V4;
ALLCOORDS = Sgeo.ALLCOORDS;
RTAB384 = Sgeo.RTAB384;
RFrange = Sgeo.RFrange(:);

nV4 = numel(RFrange);
siteRows = (1:nV4).';

% Geometry summary across all stimuli.
nTargetStim = zeros(nV4,1);
nDistrStim  = zeros(nV4,1);
nBackStim   = zeros(nV4,1);
nOverlapStim = zeros(nV4,1);

for stim = 1:numel(Tall_V4)
    T = Tall_V4(stim).T;
    assert(height(T) >= nV4, 'Tall_V4(%d).T has too few rows.', stim);

    assign = string(T.assignment(siteRows));
    nTargetStim  = nTargetStim  + (assign == "target");
    nDistrStim   = nDistrStim   + (assign == "distractor");
    nBackStim    = nBackStim    + (assign == "background");
    nOverlapStim = nOverlapStim + (T.overlap(siteRows) ~= 0);
end

nObjectStim = nTargetStim + nDistrStim;
hasObjectRF = nObjectStim >= MinObjectStim;

fprintf('V4 geometry summary for monkey %s\n', char(monkeySuffix));
fprintf('Sites with >=%d target/distractor assignments: %d / %d\n', ...
    MinObjectStim, nnz(hasObjectRF), nV4);
fprintf('Median #object assignments per V4 site: %.1f\n', median(nObjectStim));
fprintf('Median #overlap assignments per V4 site: %.1f\n', median(nOverlapStim));

% Load 3-bin response moments and localize them to the V4 rows used in Tall_V4.
R3_full = load_capsules_struct_exclusion_aware(resp3binPath, monkeySuffix, 'cfg', cfg);
R3 = R3_full;
R3.meanAct = R3_full.meanAct(RFrange, :, :);
R3.meanSqAct = R3_full.meanSqAct(RFrange, :, :);
if ismatrix(R3_full.nTrials) && size(R3_full.nTrials,1) >= max(RFrange)
    R3.nTrials = R3_full.nTrials(RFrange, :);
else
    R3.nTrials = R3_full.nTrials;
end

% Onset-response screening using the same color-specific SNR definition as V1,
% but computed on the local V4 rows.
SNR = compute_snr_per_color_sites(R3, Tall_V4, siteRows, 'Verbose', true);
SNRmat = [SNR.yellowEarly(siteRows), SNR.yellowLate(siteRows), ...
          SNR.purpleEarly(siteRows), SNR.purpleLate(siteRows)];
[bestSNR, bestIdx] = max(SNRmat, [], 2, 'omitnan'); %#ok<ASGLU>

% Attention-modulation screening in the late window.
optsTD = struct('v1Sites', siteRows, 'timeIdx', 3, 'excludeOverlap', true, 'verbose', false);
OUTtd = attention_modulation_V1_3bin(R3, Tall_V4, SNR, optsTD);

matchedN = OUTtd.wY + OUTtd.wP;
isOnsetDriven = hasObjectRF & isfinite(bestSNR) & (bestSNR > SNRthr);
isAttentionMod = hasObjectRF & isfinite(OUTtd.pValueTD) & (OUTtd.pValueTD < pTDthr) & (matchedN >= NminMatched);

keepSitesPrimary = find(isOnsetDriven);
keepSitesExpanded = find(isOnsetDriven | isAttentionMod);

fprintf('Onset-driven V4 sites (bestSNR > %.2f): %d / %d\n', SNRthr, nnz(isOnsetDriven), nV4);
fprintf('Attention-modulated-only V4 sites (pTD < %.3f, matchedN >= %d): %d\n', ...
    pTDthr, NminMatched, nnz(isAttentionMod & ~isOnsetDriven));
fprintf('Expanded V4 set (onset OR attention): %d / %d\n', numel(keepSitesExpanded), nV4);

% Quick look at strongest V4 sites.
[~, ordBest] = sort(bestSNR, 'descend', 'MissingPlacement', 'last');
kshow = min(10, nV4);
disp(table((1:kshow).', ordBest(1:kshow), RFrange(ordBest(1:kshow)), ...
    bestSNR(ordBest(1:kshow)), matchedN(ordBest(1:kshow)), ...
    'VariableNames', {'Rank','LocalV4Row','GlobalSiteInR','bestSNR','matchedN'}));

% Summary plots.
figure;
histogram(bestSNR(isfinite(bestSNR)), 30);
xlabel('Max SNR');
ylabel('N sites');
title(sprintf('V4 max SNR per site (%s, N=%d)', char(monkeySuffix), sum(isfinite(bestSNR))));
grid on;

h = plot_projected_RFs_on_example_stim(Tall_V4, ALLCOORDS, RTAB384, ExampleStimulus, ...
    'MarkerSize', 4, 'Alpha', 0.15); %#ok<NASGU>
title(sprintf('All projected V4 RFs on stim %d', ExampleStimulus));

if ~isempty(keepSitesExpanded)
    h = plot_projected_RFs_on_example_stim(Tall_V4, ALLCOORDS, RTAB384, ExampleStimulus, ...
        'MarkerSize', 5, 'Alpha', 0.25, 'SiteIdx', keepSitesExpanded); %#ok<NASGU>
    title(sprintf('Selected V4 RFs on stim %d (onset OR attention)', ExampleStimulus));
end

xAlong = [];
for stim = 1:numel(Tall_V4)
    T = Tall_V4(stim).T;
    xAlong = [xAlong; T.along_GC(keepSitesExpanded)]; %#ok<AGROW>
end
xAlong = xAlong(isfinite(xAlong));

figure;
histogram(xAlong, 'BinMethod', 'fd');
xlabel('along\_GC');
ylabel('Count');
title(sprintf('along\\_GC for selected V4 sites (%s, N=%d values)', char(monkeySuffix), numel(xAlong)));
grid on;

V4Select = struct();
V4Select.Monkey = Monkey;
V4Select.monkeySuffix = monkeySuffix;
V4Select.RFrange = RFrange;
V4Select.nTargetStim = nTargetStim;
V4Select.nDistrStim = nDistrStim;
V4Select.nBackStim = nBackStim;
V4Select.nOverlapStim = nOverlapStim;
V4Select.nObjectStim = nObjectStim;
V4Select.hasObjectRF = hasObjectRF;
V4Select.SNR = SNR;
V4Select.bestSNR = bestSNR;
V4Select.OUTtd = OUTtd;
V4Select.keepSitesPrimary = keepSitesPrimary;
V4Select.keepSitesExpanded = keepSitesExpanded;
