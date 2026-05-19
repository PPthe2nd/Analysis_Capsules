function OUT = Inspect_Prestim_TargetSide_Trials_IT()
% INSPECT_PRESTIM_TARGETSIDE_TRIALS_IT
% Trial-level diagnostic for the prestimulus offset seen in the IT
% distance-based target-side summaries.
%
% This function reuses the geometric target-closer / distractor-closer
% stimulus classes from Attention_TargetSide_Timecourse_IT, but goes back to the
% original single-trial normMUA file. For a few representative IT sites, it
% plots prestimulus trial means for the two condition classes so we can
% inspect whether the offset reflects a broad shift or a small number of
% odd trials.
%
% Notes:
% - Trial filtering matches the cached response files: correct trials only,
%   days 1 and 2.
% - The y-axis is in normMUA units from the original trial file, not in the
%   spontaneous-SD units used by the stimulus-averaged summaries.

%% Settings
P = struct();
P.Monkey = 1;                        % 1 = Nilson, 2 = Figaro
P.GroupModes = {'distance_sig', 'distance_only'};
P.nExamplesPerGroup = 2;
P.selectionQuantiles = [0.35 0.75];  % representative prestim offsets within each group
P.preWindowMs = [-200 0];
P.onlyCorrect = true;
P.correctCol = 9;
P.correctVal = 1;
P.days = [1 2];
P.dayCol = 11;
P.nHistBins = 28;
P.plotFigure = true;
P.saveResult = true;

cfg = config();

if P.Monkey == 1
    monkeySuffix = "N";
    monkeyFolder = 'Mr Nilson';
else
    monkeySuffix = "F";
    monkeyFolder = 'Figaro';
end

timecoursePath = fullfile(cfg.resultsDir, sprintf('Attention_TargetSide_Timecourse_IT_%s_d12.mat', char(monkeySuffix)));
distPath = fullfile(cfg.matDir, sprintf('Attention_TargetSide_DistanceControl_IT_%s.mat', char(monkeySuffix)));
basePath = fullfile(cfg.matDir, sprintf('Attention_TargetSide_Tuning_IT_directiondelta_%s.mat', char(monkeySuffix)));
trialRespPath = fullfile(cfg.dataRoot, monkeyFolder, 'ObjAtt_lines_normMUA.mat');
trialInfoPath = fullfile(cfg.dataRoot, monkeyFolder, 'ObjAtt_lines_MUA_trials.mat');
outPath = fullfile(cfg.resultsDir, sprintf('Inspect_Prestim_TargetSide_Trials_IT_%s.mat', char(monkeySuffix)));
hasSessionExclusions = ~isempty(site_session_exclusions(monkeySuffix));

assert(exist(trialRespPath, 'file') == 2, 'Missing %s.', trialRespPath);
assert(exist(trialInfoPath, 'file') == 2, 'Missing %s.', trialInfoPath);

%% Load cached analysis outputs
if hasSessionExclusions
    fprintf(['Session exclusions are active for monkey %s; refreshing IT target-side ' ...
             'outputs before prestim trial diagnostics.\n'], char(monkeySuffix));
    TIME = Attention_TargetSide_Timecourse_IT(struct('plotFigure', false, 'plotDifferenceFigure', false));
    DIST = Attention_TargetSide_DistanceControl_IT(struct('makeSummaryFigure', false));
    BASE = Attention_TargetSide_Tuning_IT(struct( ...
        'makeSummaryFigures', false, ...
        'makeExampleFigures', false, ...
        'makeAngleReferenceFigure', false));
else
    assert(exist(timecoursePath, 'file') == 2, ...
        'Missing %s. Run Attention_TargetSide_Timecourse_IT first.', timecoursePath);
    assert(exist(distPath, 'file') == 2, ...
        'Missing %s. Run Attention_TargetSide_DistanceControl_IT first.', distPath);
    assert(exist(basePath, 'file') == 2, ...
        'Missing %s. Run Attention_TargetSide_Tuning_IT first.', basePath);
    St = load(timecoursePath, 'OUT');
    Sd = load(distPath, 'OUT');
    Sb = load(basePath, 'OUT');
    assert(isfield(St, 'OUT') && isstruct(St.OUT), '%s must contain OUT.', timecoursePath);
    assert(isfield(Sd, 'OUT') && isstruct(Sd.OUT), '%s must contain OUT.', distPath);
    assert(isfield(Sb, 'OUT') && isstruct(Sb.OUT), '%s must contain OUT.', basePath);
    TIME = St.OUT;
    DIST = Sd.OUT;
    BASE = Sb.OUT;
end

RFrange = TIME.RFrange(:);
nIT = numel(RFrange);
assert(numel(DIST.RegressionLate) == nIT, 'Distance output does not match IT site count.');
assert(size(DIST.targetAdvPx, 1) == nIT, 'targetAdvPx rows do not match IT site count.');

tCenters = TIME.tCenters(:);
preMaskSummary = tCenters >= P.preWindowMs(1) & tCenters < P.preWindowMs(2);
assert(any(preMaskSummary), 'No summary prestim bins found in [%d %d).', P.preWindowMs(1), P.preWindowMs(2));

sitePrestimDist = mean( ...
    TIME.tcDistTargetCloser(:, preMaskSummary) - TIME.tcDistDistractorCloser(:, preMaskSummary), ...
    2, 'omitnan');
stimRef = double(BASE.QuartetTable.stimRef(:));
nQuartets = height(BASE.QuartetTable);
pairA = [stimRef, stimRef + 4];
pairB = [stimRef + 1, stimRef + 5];

%% Example-site selection from cached site-level prestim offsets
groupDefs = make_group_defs(TIME);
exampleRows = struct([]);
for g = 1:numel(P.GroupModes)
    groupName = P.GroupModes{g};
    [groupLabel, mask] = resolve_group_mask(groupDefs, groupName);
    eligible = find(mask(:) & isfinite(sitePrestimDist(:)));
    assert(~isempty(eligible), 'No eligible IT sites found for group %s.', groupName);

    qVec = choose_quantiles(P.selectionQuantiles, P.nExamplesPerGroup);
    picked = pick_sites_by_quantile(eligible, sitePrestimDist, qVec);
    for j = 1:numel(picked)
        s = picked(j);
        row = struct();
        row.groupName = string(groupName);
        row.groupLabel = string(groupLabel);
        row.siteLocal = s;
        row.siteGlobal = RFrange(s);
        row.summaryPrestimDiff = sitePrestimDist(s);
        exampleRows = [exampleRows; row]; %#ok<AGROW>
    end
end

%% Load raw trial metadata
m1 = matfile(trialRespPath);
m2 = matfile(trialInfoPath);
ALLMAT = m2.ALLMAT;
tb = double(m2.tb);
tb = tb(:);
assert(size(ALLMAT, 2) >= max([1 P.correctCol P.dayCol]), ...
    'ALLMAT does not contain the required columns.');

stimPerTrial = double(ALLMAT(:, 1));
trialInclude = true(size(stimPerTrial));
if P.onlyCorrect
    trialInclude = trialInclude & (double(ALLMAT(:, P.correctCol)) == P.correctVal);
end
if ~isempty(P.days)
    trialInclude = trialInclude & ismember(double(ALLMAT(:, P.dayCol)), P.days(:));
end

preIdx = find(tb >= P.preWindowMs(1) & tb < P.preWindowMs(2));
assert(~isempty(preIdx), 'No raw prestim samples found in [%d %d).', P.preWindowMs(1), P.preWindowMs(2));

%% Build example-site trial summaries
Examples = repmat(init_example_struct(), numel(exampleRows), 1);
for i = 1:numel(exampleRows)
    sLocal = exampleRows(i).siteLocal;
    sGlobal = exampleRows(i).siteGlobal;

    [prefStim, nonStim] = build_distance_stim_sets(sLocal, DIST, pairA, pairB, nQuartets);

    prefTrialMask = trialInclude & ismember(stimPerTrial, prefStim);
    nonTrialMask = trialInclude & ismember(stimPerTrial, nonStim);
    prefTrials = find(prefTrialMask);
    nonTrials = find(nonTrialMask);

    xAll = squeeze(mean(m1.normMUA(sGlobal, :, preIdx), 3, 'omitnan'));
    xAll = double(xAll(:));
    prefVals = xAll(prefTrials);
    nonVals = xAll(nonTrials);

    ex = init_example_struct();
    ex.groupName = exampleRows(i).groupName;
    ex.groupLabel = exampleRows(i).groupLabel;
    ex.siteLocal = sLocal;
    ex.siteGlobal = sGlobal;
    ex.summaryPrestimDiff = exampleRows(i).summaryPrestimDiff;
    ex.prefStim = prefStim(:);
    ex.nonStim = nonStim(:);
    ex.prefTrials = prefTrials(:);
    ex.nonTrials = nonTrials(:);
    ex.prefVals = prefVals(:);
    ex.nonVals = nonVals(:);
    ex.meanPref = mean(prefVals, 'omitnan');
    ex.meanNon = mean(nonVals, 'omitnan');
    ex.medianPref = median(prefVals, 'omitnan');
    ex.medianNon = median(nonVals, 'omitnan');
    ex.trimMeanPref = trimmean_local(prefVals, 20);
    ex.trimMeanNon = trimmean_local(nonVals, 20);
    ex.rankSumP = ranksum_local(prefVals, nonVals);
    ex.nOutPref = tukey_outlier_count(prefVals);
    ex.nOutNon = tukey_outlier_count(nonVals);
    Examples(i) = ex;
end

fprintf('IT prestim trial diagnostic (%s)\n', char(monkeySuffix));
for i = 1:numel(Examples)
    fprintf(['  %s | site %d (global %d): Ntar=%d Ndis=%d | ' ...
        'mean=%.4f vs %.4f | median=%.4f vs %.4f | trim20=%.4f vs %.4f | ' ...
        'outliers=%d vs %d | p=%.3g\n'], ...
        char(Examples(i).groupLabel), Examples(i).siteLocal, Examples(i).siteGlobal, ...
        numel(Examples(i).prefVals), numel(Examples(i).nonVals), ...
        Examples(i).meanPref, Examples(i).meanNon, ...
        Examples(i).medianPref, Examples(i).medianNon, ...
        Examples(i).trimMeanPref, Examples(i).trimMeanNon, ...
        Examples(i).nOutPref, Examples(i).nOutNon, ...
        Examples(i).rankSumP);
end

OUT = struct();
OUT.P = P;
OUT.monkeySuffix = monkeySuffix;
OUT.preWindowMs = P.preWindowMs;
OUT.tb = tb(preIdx);
OUT.exampleRows = exampleRows;
OUT.Examples = Examples;

if P.saveResult
    save(outPath, 'OUT', '-v7.3');
    fprintf('Saved IT prestim trial diagnostic to %s\n', outPath);
end

if P.plotFigure
    make_trial_figure(OUT);
end
end

function defs = make_group_defs(TIME)
defs = struct();
defs.distance_sig.label = 'Distance growth sig (late-pre)';
defs.distance_sig.mask = logical(TIME.isDistSig(:));
defs.distance_only.label = 'Distance growth only';
defs.distance_only.mask = logical(TIME.isDistOnly(:));
defs.direction_sig.label = 'Direction sig (dir|dist)';
defs.direction_sig.mask = logical(TIME.isDirSig(:));
defs.direction_only.label = 'Direction only';
defs.direction_only.mask = logical(TIME.isDirOnly(:));
end

function [label, mask] = resolve_group_mask(defs, name)
field = char(name);
assert(isfield(defs, field), 'Unknown group mode: %s', field);
label = defs.(field).label;
mask = defs.(field).mask;
end

function qVec = choose_quantiles(selectionQuantiles, nExamples)
if nargin < 2 || isempty(nExamples)
    nExamples = numel(selectionQuantiles);
end
if isempty(selectionQuantiles)
    if nExamples == 1
        qVec = 0.5;
    else
        qVec = linspace(0.25, 0.75, nExamples);
    end
else
    qVec = selectionQuantiles(:)';
end
if numel(qVec) ~= nExamples
    qVec = linspace(min(qVec), max(qVec), nExamples);
end
qVec = min(max(qVec, 0), 1);
end

function picked = pick_sites_by_quantile(eligible, metric, qVec)
vals = double(metric(eligible));
targets = prctile(vals, 100 * qVec);
picked = nan(size(targets));
used = false(size(eligible));
for i = 1:numel(targets)
    d = abs(vals - targets(i));
    d(used) = inf;
    [~, ix] = min(d);
    used(ix) = true;
    picked(i) = eligible(ix);
end
picked = picked(:);
end

function [prefStim, nonStim] = build_distance_stim_sets(siteLocal, DIST, pairA, pairB, nQuartets)
targetAdv = double(DIST.targetAdvPx(siteLocal, :));

prefStim = [];
nonStim = [];
for q = 1:nQuartets
    if ~(isfinite(targetAdv(q)) && targetAdv(q) ~= 0)
        continue;
    end
    if targetAdv(q) >= 0
        prefStim = [prefStim, pairA(q, :)]; %#ok<AGROW>
        nonStim = [nonStim, pairB(q, :)]; %#ok<AGROW>
    else
        prefStim = [prefStim, pairB(q, :)]; %#ok<AGROW>
        nonStim = [nonStim, pairA(q, :)]; %#ok<AGROW>
    end
end

prefStim = unique(prefStim(:));
nonStim = unique(nonStim(:));
end

function ex = init_example_struct()
ex = struct( ...
    'groupName', "", ...
    'groupLabel', "", ...
    'siteLocal', NaN, ...
    'siteGlobal', NaN, ...
    'summaryPrestimDiff', NaN, ...
    'prefStim', [], ...
    'nonStim', [], ...
    'prefTrials', [], ...
    'nonTrials', [], ...
    'prefVals', [], ...
    'nonVals', [], ...
    'meanPref', NaN, ...
    'meanNon', NaN, ...
    'medianPref', NaN, ...
    'medianNon', NaN, ...
    'trimMeanPref', NaN, ...
    'trimMeanNon', NaN, ...
    'rankSumP', NaN, ...
    'nOutPref', NaN, ...
    'nOutNon', NaN);
end

function make_trial_figure(OUT)
Examples = OUT.Examples;
nEx = numel(Examples);
cPref = [0.80 0.22 0.16];
cNon = [0.30 0.30 0.30];

fig = figure('Color', 'w', 'Name', sprintf('IT prestim trial distributions (%s)', char(OUT.monkeySuffix)));
useTiled = exist('tiledlayout', 'file') == 2;
if useTiled
    tiledlayout(nEx, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
end

for i = 1:nEx
    ex = Examples(i);
    allVals = [ex.prefVals(:); ex.nonVals(:)];
    yPad = 0.08 * max(std(allVals, 0, 'omitnan'), eps);
    yMin = min(allVals) - yPad;
    yMax = max(allVals) + yPad;

    if useTiled, ax1 = nexttile; else, ax1 = subplot(nEx, 2, 2*i-1); end %#ok<LAXES>
    hold(ax1, 'on');
    xp = 1 + 0.12 * (rand(size(ex.prefVals)) - 0.5);
    xn = 2 + 0.12 * (rand(size(ex.nonVals)) - 0.5);
    scatter(ax1, xp, ex.prefVals, 12, 'MarkerFaceColor', cPref, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.45);
    scatter(ax1, xn, ex.nonVals, 12, 'MarkerFaceColor', cNon, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.35);
    plot(ax1, [0.88 1.12], ex.meanPref * [1 1], '-', 'Color', cPref, 'LineWidth', 2.2);
    plot(ax1, [1.88 2.12], ex.meanNon * [1 1], '-', 'Color', cNon, 'LineWidth', 2.2);
    plot(ax1, [0.90 1.10], ex.medianPref * [1 1], '--', 'Color', [0.35 0 0], 'LineWidth', 1.2);
    plot(ax1, [1.90 2.10], ex.medianNon * [1 1], '--', 'Color', [0 0 0], 'LineWidth', 1.2);
    xlim(ax1, [0.5 2.5]);
    ylim(ax1, [yMin yMax]);
    set(ax1, 'XTick', [1 2], 'XTickLabel', {'tar-close', 'dis-close'});
    ylabel(ax1, 'Prestim trial mean (normMUA)');
    title(ax1, sprintf('%s | site %d (global %d)', char(ex.groupLabel), ex.siteLocal, ex.siteGlobal), ...
        'Interpreter', 'none');
    txt = sprintf('summary diff %.3f | ranksum p=%.3g | outliers %d/%d', ...
        ex.summaryPrestimDiff, ex.rankSumP, ex.nOutPref, ex.nOutNon);
    text(ax1, 0.03, 0.97, txt, 'Units', 'normalized', 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', 'FontSize', 8, 'Color', [0.20 0.20 0.20]);
    grid(ax1, 'on');

    if useTiled, ax2 = nexttile; else, ax2 = subplot(nEx, 2, 2*i); end %#ok<LAXES>
    hold(ax2, 'on');
    if numel(allVals) < 2 || ~all(isfinite(allVals))
        edges = linspace(yMin, yMax, OUT.P.nHistBins + 1);
    else
        edges = linspace(min(allVals), max(allVals), OUT.P.nHistBins + 1);
    end
    histogram(ax2, ex.nonVals, edges, 'Normalization', 'probability', ...
        'FaceColor', cNon, 'FaceAlpha', 0.28, 'EdgeColor', 'none');
    histogram(ax2, ex.prefVals, edges, 'Normalization', 'probability', ...
        'FaceColor', cPref, 'FaceAlpha', 0.32, 'EdgeColor', 'none');
    xline(ax2, ex.meanPref, '-', 'Color', cPref, 'LineWidth', 2.0);
    xline(ax2, ex.meanNon, '-', 'Color', cNon, 'LineWidth', 2.0);
    xlabel(ax2, 'Prestim trial mean (normMUA)');
    ylabel(ax2, 'Trial fraction');
    title(ax2, sprintf('Ntar=%d | Ndis=%d', numel(ex.prefVals), numel(ex.nonVals)));
    if i == 1
        legend(ax2, {'distractor-closer','target-closer'}, 'Location', 'best', 'Box', 'off');
    end
    grid(ax2, 'on');
end

annotation(fig, 'textbox', [0.05 0.965 0.90 0.035], ...
    'String', sprintf(['IT prestim trial distributions (%s). Target-closer and distractor-closer are defined only by RF geometry, ' ...
    'while each dot/histogram sample here is an original single trial from normMUA (correct trials, days 1-2).'], char(OUT.monkeySuffix)), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'FontWeight', 'bold', 'FontSize', 11);
end

function val = trimmean_local(x, pct)
x = double(x(:));
x = x(isfinite(x));
if isempty(x)
    val = NaN;
    return;
end
x = sort(x);
n = numel(x);
k = floor((pct / 100) * n / 2);
if (2 * k) >= n
    val = mean(x, 'omitnan');
else
    val = mean(x((k+1):(n-k)), 'omitnan');
end
end

function p = ranksum_local(x, y)
x = double(x(:));
y = double(y(:));
x = x(isfinite(x));
y = y(isfinite(y));
if numel(x) < 2 || numel(y) < 2
    p = NaN;
    return;
end
if exist('ranksum', 'file') == 2
    p = ranksum(x, y);
else
    [~, p] = ttest2(x, y, 'Vartype', 'unequal');
end
end

function nOut = tukey_outlier_count(x)
x = double(x(:));
x = x(isfinite(x));
if numel(x) < 4
    nOut = 0;
    return;
end
q = prctile(x, [25 75]);
iqrVal = q(2) - q(1);
lo = q(1) - 1.5 * iqrVal;
hi = q(2) + 1.5 * iqrVal;
nOut = nnz(x < lo | x > hi);
end
