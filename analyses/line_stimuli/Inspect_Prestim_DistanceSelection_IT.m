function OUT = Inspect_Prestim_DistanceSelection_IT(Puser)
% INSPECT_PRESTIM_DISTANCESELECTION_IT
% Diagnose the prestimulus offset from the earlier IT distance-based
% target-side timecourse analysis by going back to the original trial data.
%
% The condition split matches Attention_TargetSide_Timecourse_IT exactly:
% - site groups come from the current late-vs-pre distance-growth criterion
% - target-closer / distractor-closer are defined only from the geometric
%   sign of targetAdvPx
% - trial filtering matches the cached response files (correct trials, days
%   1-2)
%
% The figure separates three views for each distance-based site group:
%   1) pooled raw prestim trial distributions
%   2) pooled site-centered prestim trial distributions
%   3) site-level prestim mean difference (target-closer - distractor-closer), which is the
%      quantity most directly linked to the earlier cross-site panel

%% Settings
P = struct();
P.Monkey = 1;                        % 1 = Nilson, 2 = Figaro
P.preWindowMs = [-200 0];
P.onlyCorrect = true;
P.correctCol = 9;
P.correctVal = 1;
P.days = [1 2];
P.dayCol = 11;
P.nHistBins = 32;
P.siteChunk = 8;
P.plotFigure = true;
P.saveResult = true;
P.forceRefit = false;

if nargin >= 1 && ~isempty(Puser)
    userFields = fieldnames(Puser);
    for iF = 1:numel(userFields)
        P.(userFields{iF}) = Puser.(userFields{iF});
    end
end

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
outPath = fullfile(cfg.resultsDir, sprintf('Inspect_Prestim_DistanceSelection_IT_%s.mat', char(monkeySuffix)));
currentExclusions = site_session_exclusions(monkeySuffix);
hasSessionExclusions = ~isempty(currentExclusions);

assert(exist(trialRespPath, 'file') == 2, 'Missing %s.', trialRespPath);
assert(exist(trialInfoPath, 'file') == 2, 'Missing %s.', trialInfoPath);

useCached = exist(outPath, 'file') == 2 && ~P.forceRefit;
if useCached
    S = load(outPath, 'OUT');
    if session_exclusion_cache_matches(S, monkeySuffix)
        assert(isfield(S, 'OUT') && isstruct(S.OUT), '%s must contain OUT.', outPath);
        OUT = S.OUT;
        if P.plotFigure
            make_summary_figure(OUT);
        end
        return;
    end
    if hasSessionExclusions
        fprintf(['Cached IT prestim distance-selection diagnostic does not match the active ' ...
                 'session exclusions for monkey %s; recomputing.\n'], char(monkeySuffix));
    else
        fprintf('Cached IT prestim distance-selection diagnostic is from an exclusion-aware run; recomputing canonical cache.\n');
    end
end

%% Load cached selection logic
if hasSessionExclusions
    fprintf(['Session exclusions are active for monkey %s; refreshing IT target-side ' ...
             'outputs before prestim distance diagnostics.\n'], char(monkeySuffix));
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
    TIME = St.OUT;
    DIST = Sd.OUT;
    BASE = Sb.OUT;
end

RFrange = TIME.RFrange(:);
nIT = numel(RFrange);
assert(numel(DIST.RegressionLate) == nIT, 'Regression/site count mismatch.');

stimRef = double(BASE.QuartetTable.stimRef(:));
nQuartets = height(BASE.QuartetTable);
pairA = [stimRef, stimRef + 4];
pairB = [stimRef + 1, stimRef + 5];

groupList = { ...
    struct('name', 'Distance growth sig (late-pre)', 'mask', logical(TIME.isDistSig(:))), ...
    struct('name', 'Distance growth only', 'mask', logical(TIME.isDistOnly(:))) ...
    };

%% Load raw trial data
m1 = matfile(trialRespPath);
m2 = matfile(trialInfoPath);
ALLMAT = m2.ALLMAT;
tb = double(m2.tb);
tb = tb(:);

stimPerTrial = double(ALLMAT(:,1));
trialInclude = true(size(stimPerTrial));
if P.onlyCorrect
    trialInclude = trialInclude & (double(ALLMAT(:,P.correctCol)) == P.correctVal);
end
if ~isempty(P.days)
    trialInclude = trialInclude & ismember(double(ALLMAT(:,P.dayCol)), P.days(:));
end

preIdx = find(tb >= P.preWindowMs(1) & tb < P.preWindowMs(2));
assert(~isempty(preIdx), 'No prestim raw samples found in [%d %d).', P.preWindowMs(1), P.preWindowMs(2));

%% Collect raw-trial and site-level prestim summaries
Groups = repmat(init_group_struct(), numel(groupList), 1);
for g = 1:numel(groupList)
    mask = groupList{g}.mask(:);
    siteRows = find(mask);

    G = init_group_struct();
    G.name = string(groupList{g}.name);
    G.siteRows = siteRows(:);
    G.siteGlobal = RFrange(siteRows);
    G.nSites = numel(siteRows);
    fprintf('Processing %s: %d sites\n', char(G.name), G.nSites);

    pooledPref = [];
    pooledNon = [];
    pooledPrefCentered = [];
    pooledNonCentered = [];
    siteMeanPref = nan(G.nSites, 1);
    siteMeanNon = nan(G.nSites, 1);
    nTrialPref = zeros(G.nSites, 1);
    nTrialNon = zeros(G.nSites, 1);

    starts = 1:P.siteChunk:G.nSites;
    for k = 1:numel(starts)
        i0 = starts(k);
        i1 = min(G.nSites, i0 + P.siteChunk - 1);
        idxChunk = i0:i1;
        siteLocalChunk = siteRows(idxChunk);
        siteGlobalChunk = RFrange(siteLocalChunk);
        g0 = min(siteGlobalChunk);
        g1 = max(siteGlobalChunk);
        Xfull = m1.normMUA(g0:g1, :, preIdx);
        Xfull = squeeze(mean(Xfull, 3, 'omitnan'));
        rowIdx = siteGlobalChunk - g0 + 1;
        X = Xfull(rowIdx, :);
        if mod(k, 5) == 1 || k == numel(starts)
            fprintf('  chunk %d / %d (%d:%d)\n', k, numel(starts), i0, i1);
        end
        if numel(idxChunk) == 1
            X = reshape(double(X), 1, []);
        else
            X = double(X);
        end

        for jj = 1:numel(idxChunk)
            i = idxChunk(jj);
            sLocal = siteLocalChunk(jj);
            [prefStim, nonStim] = build_distance_stim_sets(sLocal, DIST, pairA, pairB, nQuartets);

            prefTrialMask = trialInclude & ismember(stimPerTrial, prefStim);
            nonTrialMask = trialInclude & ismember(stimPerTrial, nonStim);
            prefTrials = find(prefTrialMask);
            nonTrials = find(nonTrialMask);

            xAll = X(jj, :).';
            xPref = xAll(prefTrials);
            xNon = xAll(nonTrials);

            nTrialPref(i) = numel(xPref);
            nTrialNon(i) = numel(xNon);
            siteMeanPref(i) = mean(xPref, 'omitnan');
            siteMeanNon(i) = mean(xNon, 'omitnan');

            pooledPref = [pooledPref; xPref(:)]; %#ok<AGROW>
            pooledNon = [pooledNon; xNon(:)]; %#ok<AGROW>

            xAllSite = [xPref(:); xNon(:)];
            siteCenter = mean(xAllSite, 'omitnan');
            pooledPrefCentered = [pooledPrefCentered; xPref(:) - siteCenter]; %#ok<AGROW>
            pooledNonCentered = [pooledNonCentered; xNon(:) - siteCenter]; %#ok<AGROW>
        end
    end

    G.pooledPref = pooledPref;
    G.pooledNon = pooledNon;
    G.pooledPrefCentered = pooledPrefCentered;
    G.pooledNonCentered = pooledNonCentered;
    G.siteMeanPref = siteMeanPref;
    G.siteMeanNon = siteMeanNon;
    G.siteDiff = siteMeanPref - siteMeanNon;
    G.nTrialPref = nTrialPref;
    G.nTrialNon = nTrialNon;

    G.meanRawPref = mean(pooledPref, 'omitnan');
    G.meanRawNon = mean(pooledNon, 'omitnan');
    G.meanCenteredPref = mean(pooledPrefCentered, 'omitnan');
    G.meanCenteredNon = mean(pooledNonCentered, 'omitnan');
    G.meanSiteDiff = mean(G.siteDiff, 'omitnan');
    G.medianSiteDiff = median(G.siteDiff, 'omitnan');
    G.trimMeanSiteDiff = trimmean_local(G.siteDiff, 20);
    G.pRaw = pooled_two_sample_p_local(pooledPref, pooledNon);
    G.pCentered = pooled_two_sample_p_local(pooledPrefCentered, pooledNonCentered);
    G.pSite = signrank_local(G.siteDiff);
    G.fracSitePositive = mean(G.siteDiff > 0, 'omitnan');

    Groups(g) = G;
end

fprintf('IT prestim distance-selection diagnostic (%s)\n', char(monkeySuffix));
for g = 1:numel(Groups)
    G = Groups(g);
    fprintf(['  %s: Nsites=%d | trial N=%d vs %d | raw mean=%.4f vs %.4f | ' ...
        'centered mean=%.4f vs %.4f | site diff mean=%.4f median=%.4f trim20=%.4f | ' ...
        'pRaw=%.3g pCentered=%.3g pSite=%.3g fracPos=%.3f\n'], ...
        char(G.name), G.nSites, numel(G.pooledPref), numel(G.pooledNon), ...
        G.meanRawPref, G.meanRawNon, ...
        G.meanCenteredPref, G.meanCenteredNon, ...
        G.meanSiteDiff, G.medianSiteDiff, G.trimMeanSiteDiff, ...
        G.pRaw, G.pCentered, G.pSite, G.fracSitePositive);
end

OUT = struct();
OUT.P = P;
OUT.monkeySuffix = monkeySuffix;
OUT.preWindowMs = P.preWindowMs;
OUT.Groups = Groups;
OUT.siteSessionExclusions = currentExclusions;

if P.saveResult
    save(outPath, 'OUT', '-v7.3');
    fprintf('Saved IT prestim distance-selection diagnostic to %s\n', outPath);
end

if P.plotFigure
    make_summary_figure(OUT);
end
end

function G = init_group_struct()
G = struct( ...
    'name', "", ...
    'siteRows', [], ...
    'siteGlobal', [], ...
    'nSites', 0, ...
    'pooledPref', [], ...
    'pooledNon', [], ...
    'pooledPrefCentered', [], ...
    'pooledNonCentered', [], ...
    'siteMeanPref', [], ...
    'siteMeanNon', [], ...
    'siteDiff', [], ...
    'nTrialPref', [], ...
    'nTrialNon', [], ...
    'meanRawPref', NaN, ...
    'meanRawNon', NaN, ...
    'meanCenteredPref', NaN, ...
    'meanCenteredNon', NaN, ...
    'meanSiteDiff', NaN, ...
    'medianSiteDiff', NaN, ...
    'trimMeanSiteDiff', NaN, ...
    'pRaw', NaN, ...
    'pCentered', NaN, ...
    'pSite', NaN, ...
    'fracSitePositive', NaN);
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
        prefStim = [prefStim, pairA(q,:)]; %#ok<AGROW>
        nonStim = [nonStim, pairB(q,:)]; %#ok<AGROW>
    else
        prefStim = [prefStim, pairB(q,:)]; %#ok<AGROW>
        nonStim = [nonStim, pairA(q,:)]; %#ok<AGROW>
    end
end

prefStim = unique(prefStim(:));
nonStim = unique(nonStim(:));
end

function make_summary_figure(OUT)
Groups = OUT.Groups;
nG = numel(Groups);
cPref = [0.82 0.22 0.16];
cNon = [0.30 0.30 0.30];

fig = figure('Color', 'w', 'Name', sprintf('IT prestim distance-selection (%s)', char(OUT.monkeySuffix)));
useTiled = exist('tiledlayout', 'file') == 2;
if useTiled
    tiledlayout(nG, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
end

for g = 1:nG
    G = Groups(g);

    if useTiled, ax1 = nexttile; else, ax1 = subplot(nG, 3, 3*(g-1)+1); end %#ok<LAXES>
    hold(ax1, 'on');
    edgesRaw = choose_edges([G.pooledPref; G.pooledNon], OUT.P.nHistBins);
    histogram(ax1, G.pooledNon, edgesRaw, 'Normalization', 'probability', ...
        'FaceColor', cNon, 'FaceAlpha', 0.28, 'EdgeColor', 'none');
    histogram(ax1, G.pooledPref, edgesRaw, 'Normalization', 'probability', ...
        'FaceColor', cPref, 'FaceAlpha', 0.34, 'EdgeColor', 'none');
    xline(ax1, G.meanRawNon, '-', 'Color', cNon, 'LineWidth', 2.0);
    xline(ax1, G.meanRawPref, '-', 'Color', cPref, 'LineWidth', 2.0);
    xlabel(ax1, 'Prestim trial mean (raw normMUA)');
    ylabel(ax1, 'Trial fraction');
    title(ax1, sprintf('%s: pooled raw', char(G.name)), 'Interpreter', 'none');
    txtRaw = sprintf('N=%d vs %d | p=%.3g', numel(G.pooledPref), numel(G.pooledNon), G.pRaw);
    text(ax1, 0.03, 0.97, txtRaw, 'Units', 'normalized', 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', 'FontSize', 8, 'Color', [0.2 0.2 0.2]);
    if g == 1
        legend(ax1, {'distractor-closer','target-closer'}, 'Location', 'best', 'Box', 'off');
    end
    grid(ax1, 'on');

    if useTiled, ax2 = nexttile; else, ax2 = subplot(nG, 3, 3*(g-1)+2); end %#ok<LAXES>
    hold(ax2, 'on');
    edgesCtr = choose_edges([G.pooledPrefCentered; G.pooledNonCentered], OUT.P.nHistBins);
    histogram(ax2, G.pooledNonCentered, edgesCtr, 'Normalization', 'probability', ...
        'FaceColor', cNon, 'FaceAlpha', 0.28, 'EdgeColor', 'none');
    histogram(ax2, G.pooledPrefCentered, edgesCtr, 'Normalization', 'probability', ...
        'FaceColor', cPref, 'FaceAlpha', 0.34, 'EdgeColor', 'none');
    xline(ax2, G.meanCenteredNon, '-', 'Color', cNon, 'LineWidth', 2.0);
    xline(ax2, G.meanCenteredPref, '-', 'Color', cPref, 'LineWidth', 2.0);
    xline(ax2, 0, 'k:');
    xlabel(ax2, 'Prestim trial mean (site-centered)');
    ylabel(ax2, 'Trial fraction');
    title(ax2, sprintf('%s: pooled centered', char(G.name)), 'Interpreter', 'none');
    txtCtr = sprintf('p=%.3g', G.pCentered);
    text(ax2, 0.03, 0.97, txtCtr, 'Units', 'normalized', 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', 'FontSize', 8, 'Color', [0.2 0.2 0.2]);
    grid(ax2, 'on');

    if useTiled, ax3 = nexttile; else, ax3 = subplot(nG, 3, 3*(g-1)+3); end %#ok<LAXES>
    hold(ax3, 'on');
    edgesSite = choose_edges(G.siteDiff, max(16, ceil(sqrt(max(G.nSites, 1)))));
    histogram(ax3, G.siteDiff, edgesSite, 'FaceColor', [0.40 0.55 0.80], 'EdgeColor', 'none');
    xline(ax3, 0, 'k:');
    xline(ax3, G.meanSiteDiff, '-', 'Color', [0.10 0.25 0.55], 'LineWidth', 2.0);
    xlabel(ax3, 'Site prestim mean diff (target-closer - distractor-closer)');
    ylabel(ax3, 'N sites');
    title(ax3, sprintf('%s: site means', char(G.name)), 'Interpreter', 'none');
    txtSite = sprintf('N=%d | mean=%.4f | median=%.4f | p=%.3g | frac>0=%.2f', ...
        G.nSites, G.meanSiteDiff, G.medianSiteDiff, G.pSite, G.fracSitePositive);
    text(ax3, 0.03, 0.97, txtSite, 'Units', 'normalized', 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', 'FontSize', 8, 'Color', [0.2 0.2 0.2]);
    grid(ax3, 'on');
end

annotation(fig, 'textbox', [0.05 0.965 0.90 0.035], ...
    'String', sprintf(['IT prestim distributions under the exact target-closer vs distractor-closer geometric split used in the current timecourse analysis (%s). ' ...
    'Left: pooled raw trial values. Middle: pooled values after removing each site''s overall prestim baseline. Right: site-level prestim mean difference, ' ...
    'which most directly corresponds to the earlier suspicious cross-site panel.'], char(OUT.monkeySuffix)), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontWeight', 'bold', 'FontSize', 11);
end

function edges = choose_edges(x, nBins)
x = double(x(:));
x = x(isfinite(x));
if isempty(x)
    edges = linspace(-1, 1, max(nBins, 2) + 1);
    return;
end
if numel(unique(x)) < 2
    dx = max(abs(x(1)), 1) * 0.05 + 1e-3;
    edges = linspace(x(1) - dx, x(1) + dx, max(nBins, 2) + 1);
    return;
end
edges = linspace(min(x), max(x), max(nBins, 2) + 1);
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

function p = pooled_two_sample_p_local(x, y)
x = double(x(:));
y = double(y(:));
x = x(isfinite(x));
y = y(isfinite(y));
if numel(x) < 2 || numel(y) < 2
    p = NaN;
    return;
end
[~, p] = ttest2(x, y, 'Vartype', 'unequal');
end

function p = signrank_local(x)
x = double(x(:));
x = x(isfinite(x));
if numel(x) < 3
    p = NaN;
    return;
end
if exist('signrank', 'file') == 2
    p = signrank(x, 0);
else
    [~, p] = ttest(x, 0);
end
end
