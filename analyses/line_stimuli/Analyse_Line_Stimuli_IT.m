% First-pass IT analysis for the line-task stimulus set.
%
% This script intentionally keeps the scope narrow:
% - geometry QC from Tall_IT
% - signed responsiveness across all stimuli in the standard 3-bin windows
% - signed responsiveness per quartet
% - quick plots and site-count summaries that keep suppressive sites visible

Monkey = 1; % 1 = Nilson, 2 = Figaro
ExampleStimulus = 38; % stimulus used for quick overlay/QC plots
MinObjectStim = 1; % site must fall on target/distractor in at least this many stimuli
RespThr = 0.7; % response-magnitude threshold for summary counts only
TopKQuartets = 5; % responsiveness score uses mean abs(response) of the top-k quartets

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
    error('Analyse_Line_Stimuli_IT:InvalidMonkey', ...
        'Monkey must be 1 (Nilson) or 2 (Figaro).');
end

tallPath = fullfile(cfg.matDir, tallFile);
resp3binPath = fullfile(cfg.matDir, resp3binFile);

assert(exist(tallPath, 'file') == 2, ...
    'Missing %s. Run Line_Stimuli_IT.m first.', tallPath);
assert(exist(resp3binPath, 'file') == 2, ...
    'Missing %s. Create the 3-bin response summary for this monkey first.', resp3binPath);

Sgeo = load(tallPath);
assert(isfield(Sgeo, 'Tall_IT') && isstruct(Sgeo.Tall_IT), ...
    '%s must contain struct Tall_IT.', tallFile);
assert(isfield(Sgeo, 'ALLCOORDS') && isfield(Sgeo, 'RTAB384'), ...
    '%s must contain ALLCOORDS and RTAB384.', tallFile);
assert(isfield(Sgeo, 'RFrange') && ~isempty(Sgeo.RFrange), ...
    '%s must contain non-empty RFrange.', tallFile);

Tall_IT = Sgeo.Tall_IT;
ALLCOORDS = Sgeo.ALLCOORDS;
RTAB384 = Sgeo.RTAB384;
RFrange = Sgeo.RFrange(:);

nIT = numel(RFrange);
siteRows = (1:nIT).';

% Sort by stimNum so geometry and responses align explicitly to 1..nStim.
stimNums = arrayfun(@(x) x.stimNum, Tall_IT(:));
[stimNumsSorted, ordStim] = sort(stimNums(:));
assert(all(stimNumsSorted(:).' == 1:numel(Tall_IT)), ...
    'Tall_IT.stimNum must cover 1..%d exactly.', numel(Tall_IT));
Tall_IT = Tall_IT(ordStim);
nStim = numel(Tall_IT);

% Geometry summary across all stimuli.
nTargetStim = zeros(nIT,1);
nDistrStim  = zeros(nIT,1);
nBackStim   = zeros(nIT,1);
nOverlapStim = zeros(nIT,1);

for stim = 1:nStim
    T = Tall_IT(stim).T;
    assert(height(T) >= nIT, 'Tall_IT(%d).T has too few rows.', stim);

    assign = string(T.assignment(siteRows));
    nTargetStim  = nTargetStim  + (assign == "target");
    nDistrStim   = nDistrStim   + (assign == "distractor");
    nBackStim    = nBackStim    + (assign == "background");
    nOverlapStim = nOverlapStim + (T.overlap(siteRows) ~= 0);
end

nObjectStim = nTargetStim + nDistrStim;
hasObjectRF = nObjectStim >= MinObjectStim;

fprintf('IT geometry summary for monkey %s\n', char(monkeySuffix));
fprintf('Sites with >=%d target/distractor assignments: %d / %d\n', ...
    MinObjectStim, nnz(hasObjectRF), nIT);
fprintf('Median #object assignments per IT site: %.1f\n', median(nObjectStim));
fprintf('Median #overlap assignments per IT site: %.1f\n', median(nOverlapStim));

% Load 3-bin response moments and localize them to the IT rows used in Tall_IT.
R3_full = load_capsules_struct_exclusion_aware(resp3binPath, monkeySuffix, 'cfg', cfg);
R3 = R3_full;
R3.meanAct = R3_full.meanAct(RFrange, :, :);
R3.meanSqAct = R3_full.meanSqAct(RFrange, :, :);
if ismatrix(R3_full.nTrials) && size(R3_full.nTrials,1) >= max(RFrange)
    R3.nTrials = R3_full.nTrials(RFrange, :);
else
    R3.nTrials = R3_full.nTrials;
end

assert(size(R3.meanAct,1) == nIT, 'Localized IT response rows do not match Tall_IT.');
assert(size(R3.meanAct,2) == nStim, 'Localized IT response stimuli do not match Tall_IT.');

% Reuse the standard by-color SNR helper for spontaneous mean/noise.
SNR = compute_snr_per_color_sites(R3, Tall_IT, siteRows, 'Verbose', true);

% Signed response across all stimuli.
muAllEarly = nan(nIT,1);
muAllLate = nan(nIT,1);
signedAllEarly = nan(nIT,1);
signedAllLate = nan(nIT,1);

if isvector(R3.nTrials)
    nTrialsAll = double(R3.nTrials(:));
    perSiteTrials = false;
else
    perSiteTrials = true;
end

for iSite = 1:nIT
    if perSiteTrials
        nTr = double(R3.nTrials(iSite, :)).';
    else
        nTr = nTrialsAll;
    end

    rEarly = squeeze(R3.meanAct(iSite,:,2)).';
    rLate  = squeeze(R3.meanAct(iSite,:,3)).';

    idxEarly = isfinite(rEarly) & isfinite(nTr) & (nTr > 0);
    idxLate  = isfinite(rLate)  & isfinite(nTr) & (nTr > 0);

    if any(idxEarly)
        muAllEarly(iSite) = sum(nTr(idxEarly) .* rEarly(idxEarly)) / sum(nTr(idxEarly));
    end
    if any(idxLate)
        muAllLate(iSite) = sum(nTr(idxLate) .* rLate(idxLate)) / sum(nTr(idxLate));
    end
end

muSpont = SNR.muSpont(siteRows);
sdSpont = SNR.sdSpont(siteRows);
badNoise = ~isfinite(sdSpont) | (sdSpont <= 0);

signedAllEarly = (muAllEarly - muSpont) ./ sdSpont;
signedAllLate  = (muAllLate  - muSpont) ./ sdSpont;
signedAllEarly(badNoise) = NaN;
signedAllLate(badNoise) = NaN;

% Per-quartet signed responses using the same quartet logic as the V1/V4 attention code.
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
muQuartetLate  = nan(nIT, nQuartets);
signedQuartetEarly = nan(nIT, nQuartets);
signedQuartetLate  = nan(nIT, nQuartets);

for iSite = 1:nIT
    if perSiteTrials
        nTrSite = double(R3.nTrials(iSite, :)).';
    else
        nTrSite = nTrialsAll;
    end

    rEarly = squeeze(R3.meanAct(iSite,:,2)).';
    rLate  = squeeze(R3.meanAct(iSite,:,3)).';

    for qIdx = 1:nQuartets
        stimQ = quartetMembers(qIdx,:);

        nTrQ = nTrSite(stimQ);
        rEarlyQ = rEarly(stimQ);
        rLateQ = rLate(stimQ);

        idxEarly = isfinite(rEarlyQ) & isfinite(nTrQ) & (nTrQ > 0);
        idxLate  = isfinite(rLateQ)  & isfinite(nTrQ) & (nTrQ > 0);

        if any(idxEarly)
            muQuartetEarly(iSite, qIdx) = sum(nTrQ(idxEarly) .* rEarlyQ(idxEarly)) / sum(nTrQ(idxEarly));
        end
        if any(idxLate)
            muQuartetLate(iSite, qIdx) = sum(nTrQ(idxLate) .* rLateQ(idxLate)) / sum(nTrQ(idxLate));
        end
    end
end

signedQuartetEarly = bsxfun(@rdivide, bsxfun(@minus, muQuartetEarly, muSpont), sdSpont);
signedQuartetLate  = bsxfun(@rdivide, bsxfun(@minus, muQuartetLate,  muSpont), sdSpont);
signedQuartetEarly(badNoise, :) = NaN;
signedQuartetLate(badNoise, :) = NaN;

[bestPosQuartetEarly, bestPosQuartetIdxEarly] = max(signedQuartetEarly, [], 2, 'omitnan');
[bestNegQuartetEarly, bestNegQuartetIdxEarly] = min(signedQuartetEarly, [], 2, 'omitnan');
[bestPosQuartetLate, bestPosQuartetIdxLate] = max(signedQuartetLate, [], 2, 'omitnan');
[bestNegQuartetLate, bestNegQuartetIdxLate] = min(signedQuartetLate, [], 2, 'omitnan');

bestAbsQuartetEarly = max(abs(signedQuartetEarly), [], 2, 'omitnan');
bestAbsQuartetLate  = max(abs(signedQuartetLate),  [], 2, 'omitnan');

topKAbsQuartetEarly = nan(nIT,1);
topKAbsQuartetLate = nan(nIT,1);
for iSite = 1:nIT
    valsEarly = abs(signedQuartetEarly(iSite, :));
    valsEarly = valsEarly(isfinite(valsEarly));
    if ~isempty(valsEarly)
        valsEarly = sort(valsEarly, 'descend');
        kEarly = min(TopKQuartets, numel(valsEarly));
        topKAbsQuartetEarly(iSite) = mean(valsEarly(1:kEarly));
    end

    valsLate = abs(signedQuartetLate(iSite, :));
    valsLate = valsLate(isfinite(valsLate));
    if ~isempty(valsLate)
        valsLate = sort(valsLate, 'descend');
        kLate = min(TopKQuartets, numel(valsLate));
        topKAbsQuartetLate(iSite) = mean(valsLate(1:kLate));
    end
end

quartetRespScore = max(topKAbsQuartetEarly, topKAbsQuartetLate);
isResponsiveQuartet = hasObjectRF & isfinite(quartetRespScore) & (quartetRespScore > RespThr);

isExcEarly = isfinite(signedAllEarly) & (signedAllEarly > RespThr);
isSuppEarly = isfinite(signedAllEarly) & (signedAllEarly < -RespThr);
isExcLate = isfinite(signedAllLate) & (signedAllLate > RespThr);
isSuppLate = isfinite(signedAllLate) & (signedAllLate < -RespThr);

fprintf('All-stim IT response counts in SNR units (threshold %.2f)\n', RespThr);
fprintf('Early excitatory: %d | Early suppressive: %d\n', nnz(isExcEarly), nnz(isSuppEarly));
fprintf('Late excitatory:  %d | Late suppressive:  %d\n', nnz(isExcLate), nnz(isSuppLate));
fprintf('Quartet-responsive IT sites (top-%d abs quartet score > %.2f): %d / %d\n', ...
    TopKQuartets, RespThr, nnz(isResponsiveQuartet), nIT);

[~, ordLatePos] = sort(signedAllLate, 'descend', 'MissingPlacement', 'last');
[~, ordLateNeg] = sort(signedAllLate, 'ascend', 'MissingPlacement', 'last');
kshow = min(10, nIT);

disp('Top late excitatory IT sites (all stimuli):');
disp(table((1:kshow).', ordLatePos(1:kshow), RFrange(ordLatePos(1:kshow)), ...
    signedAllLate(ordLatePos(1:kshow)), bestPosQuartetLate(ordLatePos(1:kshow)), ...
    bestPosQuartetIdxLate(ordLatePos(1:kshow)), ...
    'VariableNames', {'Rank','LocalITRow','GlobalSiteInR','signedLateAll','bestLateQuartet','bestLateQuartetIdx'}));

disp('Top late suppressive IT sites (all stimuli):');
disp(table((1:kshow).', ordLateNeg(1:kshow), RFrange(ordLateNeg(1:kshow)), ...
    signedAllLate(ordLateNeg(1:kshow)), bestNegQuartetLate(ordLateNeg(1:kshow)), ...
    bestNegQuartetIdxLate(ordLateNeg(1:kshow)), ...
    'VariableNames', {'Rank','LocalITRow','GlobalSiteInR','signedLateAll','worstLateQuartet','worstLateQuartetIdx'}));

% Summary plots.
figure;
histogram(signedAllEarly(isfinite(signedAllEarly)), 30, 'FaceColor', [0.25 0.45 0.85], 'EdgeColor', 'none');
hold on;
histogram(signedAllLate(isfinite(signedAllLate)), 30, 'FaceColor', [0.85 0.35 0.25], 'EdgeColor', 'none');
xlabel('Signed response ((mu - spont) / sdSpont)');
ylabel('N sites');
title(sprintf('IT all-stim signed responses (%s)', char(monkeySuffix)));
legend('Early','Late');
grid on;

figure;
scatter(signedAllEarly, signedAllLate, 24, nObjectStim, 'filled');
xlabel('Signed early response');
ylabel('Signed late response');
title(sprintf('IT response summary by site (%s)', char(monkeySuffix)));
cb = colorbar; %#ok<NASGU>
yline(0, '--k');
xline(0, '--k');
grid on;

figure;
histogram(bestPosQuartetLate(isfinite(bestPosQuartetLate)), 30, 'FaceColor', [0.80 0.20 0.20], 'EdgeColor', 'none');
hold on;
histogram(bestNegQuartetLate(isfinite(bestNegQuartetLate)), 30, 'FaceColor', [0.20 0.35 0.80], 'EdgeColor', 'none');
xlabel('Signed late response per quartet');
ylabel('N sites');
title(sprintf('IT best and worst late quartet responses (%s)', char(monkeySuffix)));
legend('Best positive quartet','Best negative quartet');
grid on;

figure;
% Distribution of naive GC positions for response-selected IT combinations.
xAlong = [];
for stim = 1:nStim
    T = Tall_IT(stim).T;
    xAlong = [xAlong; T.along_GC(isResponsiveQuartet)]; %#ok<AGROW>
end
xAlong = xAlong(isfinite(xAlong));

histogram(xAlong, 'BinMethod', 'fd');
xlabel('along\_GC');
ylabel('Count');
title(sprintf('IT along\\_GC for responsive sites (%s, top-%d > %.2f, N=%d values)', ...
    char(monkeySuffix), TopKQuartets, RespThr, numel(xAlong)));
grid on;

h = plot_projected_RFs_on_example_stim(Tall_IT, ALLCOORDS, RTAB384, ExampleStimulus, ...
    'MarkerSize', 4, 'Alpha', 0.15); %#ok<NASGU>
title(sprintf('All projected IT RFs on stim %d', ExampleStimulus));

if any(hasObjectRF)
    h = plot_projected_RFs_on_example_stim(Tall_IT, ALLCOORDS, RTAB384, ExampleStimulus, ...
        'MarkerSize', 5, 'Alpha', 0.25, 'SiteIdx', find(hasObjectRF)); %#ok<NASGU>
    title(sprintf('IT RFs with object assignments on stim %d', ExampleStimulus));
end

ITSummary = struct();
ITSummary.Monkey = Monkey;
ITSummary.monkeySuffix = monkeySuffix;
ITSummary.RFrange = RFrange;
ITSummary.nTargetStim = nTargetStim;
ITSummary.nDistrStim = nDistrStim;
ITSummary.nBackStim = nBackStim;
ITSummary.nOverlapStim = nOverlapStim;
ITSummary.nObjectStim = nObjectStim;
ITSummary.hasObjectRF = hasObjectRF;
ITSummary.TopKQuartets = TopKQuartets;
ITSummary.topKAbsQuartetEarly = topKAbsQuartetEarly;
ITSummary.topKAbsQuartetLate = topKAbsQuartetLate;
ITSummary.quartetRespScore = quartetRespScore;
ITSummary.isResponsiveQuartet = isResponsiveQuartet;
ITSummary.xAlong = xAlong;
ITSummary.SNR = SNR;
ITSummary.muAllEarly = muAllEarly;
ITSummary.muAllLate = muAllLate;
ITSummary.signedAllEarly = signedAllEarly;
ITSummary.signedAllLate = signedAllLate;
ITSummary.quartetMembers = quartetMembers;
ITSummary.muQuartetEarly = muQuartetEarly;
ITSummary.muQuartetLate = muQuartetLate;
ITSummary.signedQuartetEarly = signedQuartetEarly;
ITSummary.signedQuartetLate = signedQuartetLate;
ITSummary.bestPosQuartetEarly = bestPosQuartetEarly;
ITSummary.bestNegQuartetEarly = bestNegQuartetEarly;
ITSummary.bestPosQuartetLate = bestPosQuartetLate;
ITSummary.bestNegQuartetLate = bestNegQuartetLate;
ITSummary.bestPosQuartetIdxEarly = bestPosQuartetIdxEarly;
ITSummary.bestNegQuartetIdxEarly = bestNegQuartetIdxEarly;
ITSummary.bestPosQuartetIdxLate = bestPosQuartetIdxLate;
ITSummary.bestNegQuartetIdxLate = bestNegQuartetIdxLate;
ITSummary.bestAbsQuartetEarly = bestAbsQuartetEarly;
ITSummary.bestAbsQuartetLate = bestAbsQuartetLate;
ITSummary.RespThr = RespThr;
