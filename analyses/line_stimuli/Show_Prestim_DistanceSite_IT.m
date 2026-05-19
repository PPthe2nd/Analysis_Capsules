function OUT = Show_Prestim_DistanceSite_IT(Puser)
% SHOW_PRESTIM_DISTANCESITE_IT
% Single-site trial-level diagnostic for a large prestim distance effect in
% IT. The target-closer / distractor-closer classes match the current
% geometry-only distance analysis exactly, and the selected site groups now
% come from the late-vs-pre distance-growth criterion.

%% Settings
P = struct();
P.Monkey = 1;                        % 1 = Nilson, 2 = Figaro
P.preWindowMs = [-200 0];
P.onlyCorrect = true;
P.correctCol = 9;
P.correctVal = 1;
P.days = [1 2];
P.dayCol = 11;
P.plotFigure = true;
P.saveResult = true;
P.nTopOutliers = 5;
P.chooseLargestPositive = true;
P.siteLocal = [];
P.siteGlobal = [];
P.plotSessionFigure = true;

if nargin >= 1 && ~isempty(Puser)
    fn = fieldnames(Puser);
    for i = 1:numel(fn)
        P.(fn{i}) = Puser.(fn{i});
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

distInspectPath = fullfile(cfg.resultsDir, sprintf('Inspect_Prestim_DistanceSelection_IT_%s.mat', char(monkeySuffix)));
distPath = fullfile(cfg.matDir, sprintf('Attention_TargetSide_DistanceControl_IT_%s.mat', char(monkeySuffix)));
basePath = fullfile(cfg.matDir, sprintf('Attention_TargetSide_Tuning_IT_directiondelta_%s.mat', char(monkeySuffix)));
trialRespPath = fullfile(cfg.dataRoot, monkeyFolder, 'ObjAtt_lines_normMUA.mat');
trialInfoPath = fullfile(cfg.dataRoot, monkeyFolder, 'ObjAtt_lines_MUA_trials.mat');
outPath = fullfile(cfg.resultsDir, sprintf('Show_Prestim_DistanceSite_IT_%s.mat', char(monkeySuffix)));
hasSessionExclusions = ~isempty(site_session_exclusions(monkeySuffix));

assert(exist(trialRespPath, 'file') == 2, 'Missing %s.', trialRespPath);
assert(exist(trialInfoPath, 'file') == 2, 'Missing %s.', trialInfoPath);

%% Load site ranking from previous distance-selection diagnostic
if hasSessionExclusions
    SEL = Inspect_Prestim_DistanceSelection_IT(struct('plotFigure', false, 'saveResult', false, 'Monkey', P.Monkey));
elseif exist(distInspectPath, 'file') == 2
    Ssel = load(distInspectPath, 'OUT');
    SEL = Ssel.OUT;
else
    SEL = Inspect_Prestim_DistanceSelection_IT(struct('plotFigure', false, 'saveResult', false, 'Monkey', P.Monkey));
end

if hasSessionExclusions
    fprintf(['Session exclusions are active for monkey %s; refreshing IT target-side ' ...
             'outputs before single-site prestim diagnostics.\n'], char(monkeySuffix));
    DIST = Attention_TargetSide_DistanceControl_IT(struct('makeSummaryFigure', false));
    BASE = Attention_TargetSide_Tuning_IT(struct( ...
        'makeSummaryFigures', false, ...
        'makeExampleFigures', false, ...
        'makeAngleReferenceFigure', false));
else
    assert(exist(distPath, 'file') == 2, 'Missing %s.', distPath);
    assert(exist(basePath, 'file') == 2, 'Missing %s.', basePath);
    Sd = load(distPath, 'OUT');
    Sb = load(basePath, 'OUT');
    DIST = Sd.OUT;
    BASE = Sb.OUT;
end
RFrange = DIST.RFrange(:);

if isempty(P.siteLocal) && isempty(P.siteGlobal)
    [siteLocal, siteGlobal, siteGroup, siteDiff] = choose_site_from_selection(SEL, RFrange, P.chooseLargestPositive);
else
    if ~isempty(P.siteLocal)
        siteLocal = P.siteLocal;
        siteGlobal = RFrange(siteLocal);
    else
        siteGlobal = P.siteGlobal;
        siteLocal = find(RFrange == siteGlobal, 1, 'first');
        assert(~isempty(siteLocal), 'Requested global site %d not found in IT RFrange.', siteGlobal);
    end
    [siteGroup, siteDiff] = lookup_site_in_selection(SEL, siteLocal);
end

stimRef = double(BASE.QuartetTable.stimRef(:));
nQuartets = height(BASE.QuartetTable);
pairA = [stimRef, stimRef + 4];
pairB = [stimRef + 1, stimRef + 5];
[prefStim, nonStim] = build_distance_stim_sets(siteLocal, DIST, pairA, pairB, nQuartets);

%% Load raw trial data for this site
m1 = matfile(trialRespPath);
m2 = matfile(trialInfoPath);
ALLMAT = m2.ALLMAT;
tb = double(m2.tb);
tb = tb(:);
preIdx = find(tb >= P.preWindowMs(1) & tb < P.preWindowMs(2));
assert(~isempty(preIdx), 'No prestim raw samples found in [%d %d).', P.preWindowMs(1), P.preWindowMs(2));

stimPerTrial = double(ALLMAT(:,1));
sessionPerTrial = double(ALLMAT(:,P.dayCol));
trialInclude = true(size(stimPerTrial));
if P.onlyCorrect
    trialInclude = trialInclude & (double(ALLMAT(:,P.correctCol)) == P.correctVal);
end
if ~isempty(P.days)
    trialInclude = trialInclude & ismember(sessionPerTrial, P.days(:));
end

xAll = squeeze(mean(m1.normMUA(siteGlobal, :, preIdx), 3, 'omitnan'));
xAll = double(xAll(:));
tcAll = squeeze(m1.normMUA(siteGlobal, :, :));
if isvector(tcAll)
    tcAll = reshape(double(tcAll), 1, []);
else
    tcAll = double(tcAll);
end

prefMask = trialInclude & ismember(stimPerTrial, prefStim);
nonMask = trialInclude & ismember(stimPerTrial, nonStim);
prefTrials = find(prefMask);
nonTrials = find(nonMask);
allUsed = sort([prefTrials; nonTrials]);
usedVals = xAll(allUsed);
usedStim = stimPerTrial(allUsed);
usedSess = sessionPerTrial(allUsed);
isPrefUsed = prefMask(allUsed);

prefVals = xAll(prefTrials);
nonVals = xAll(nonTrials);
[outPref, prefFence] = tukey_outlier_mask(prefVals);
[outNon, nonFence] = tukey_outlier_mask(nonVals);

meanPref = mean(prefVals, 'omitnan');
meanNon = mean(nonVals, 'omitnan');
medPref = median(prefVals, 'omitnan');
medNon = median(nonVals, 'omitnan');
trimPref = trimmean_local(prefVals, 20);
trimNon = trimmean_local(nonVals, 20);
meanDiff = meanPref - meanNon;
meanDiffNoOut = mean(prefVals(~outPref), 'omitnan') - mean(nonVals(~outNon), 'omitnan');
pRank = ranksum_local(prefVals, nonVals);

sessionIds = unique(sessionPerTrial(trialInclude));
sessionIds = sessionIds(:).';
nSess = numel(sessionIds);
sessionMeanTc = nan(nSess, numel(tb));
sessionSemTc = nan(nSess, numel(tb));
sessionN = zeros(nSess,1);
for iSess = 1:nSess
    idx = trialInclude & (sessionPerTrial == sessionIds(iSess));
    X = tcAll(idx, :);
    if isempty(X)
        continue;
    end
    sessionN(iSess) = size(X, 1);
    sessionMeanTc(iSess, :) = mean(X, 1, 'omitnan');
    sessionSemTc(iSess, :) = std(X, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(X), 1));
end

summaryLines = build_summary_lines(siteGroup, siteLocal, siteGlobal, siteDiff, ...
    prefVals, nonVals, outPref, outNon, meanDiff, meanDiffNoOut, pRank);
outlierLines = build_outlier_lines(prefTrials, prefVals, stimPerTrial(prefTrials), sessionPerTrial(prefTrials), outPref, ...
    nonTrials, nonVals, stimPerTrial(nonTrials), sessionPerTrial(nonTrials), outNon, P.nTopOutliers);

OUT = struct();
OUT.P = P;
OUT.monkeySuffix = monkeySuffix;
OUT.siteLocal = siteLocal;
OUT.siteGlobal = siteGlobal;
OUT.siteGroup = siteGroup;
OUT.summaryPrestimDiff = siteDiff;
OUT.prefStim = prefStim(:);
OUT.nonStim = nonStim(:);
OUT.prefTrials = prefTrials(:);
OUT.nonTrials = nonTrials(:);
OUT.prefVals = prefVals(:);
OUT.nonVals = nonVals(:);
OUT.outPref = outPref(:);
OUT.outNon = outNon(:);
OUT.prefFence = prefFence;
OUT.nonFence = nonFence;
OUT.allUsedTrials = allUsed(:);
OUT.usedVals = usedVals(:);
OUT.usedStim = usedStim(:);
OUT.usedSess = usedSess(:);
OUT.isPrefUsed = isPrefUsed(:);
OUT.tb = tb(:);
OUT.sessionIds = sessionIds(:);
OUT.sessionMeanTc = sessionMeanTc;
OUT.sessionSemTc = sessionSemTc;
OUT.sessionN = sessionN;
OUT.summaryLines = summaryLines;
OUT.outlierLines = outlierLines;

fprintf('Chosen IT site %d (global %d) from %s | summary prestim diff %.4f\n', ...
    siteLocal, siteGlobal, char(siteGroup), siteDiff);
fprintf('  Tar-close vs dis-close trials: %d vs %d | mean diff %.4f | no-outlier diff %.4f | p=%.3g\n', ...
    numel(prefVals), numel(nonVals), meanDiff, meanDiffNoOut, pRank);
fprintf('  Outliers: tar=%d, dis=%d\n', nnz(outPref), nnz(outNon));

if P.saveResult
    save(outPath, 'OUT', '-v7.3');
    fprintf('Saved IT single-site prestim diagnostic to %s\n', outPath);
end

if P.plotFigure
    make_single_site_figure(OUT);
end
if P.plotSessionFigure
    make_session_timecourse_figure(OUT);
end
end

function [siteLocal, siteGlobal, groupName, siteDiff] = choose_site_from_selection(SEL, RFrange, chooseLargestPositive)
allLocal = [];
allGlobal = [];
allDiff = [];
allGroup = strings(0,1);
for g = 1:numel(SEL.Groups)
    G = SEL.Groups(g);
    allLocal = [allLocal; G.siteRows(:)]; %#ok<AGROW>
    allGlobal = [allGlobal; RFrange(G.siteRows(:))]; %#ok<AGROW>
    allDiff = [allDiff; G.siteDiff(:)]; %#ok<AGROW>
    allGroup = [allGroup; repmat(string(G.name), numel(G.siteRows), 1)]; %#ok<AGROW>
end
ok = isfinite(allDiff);
allLocal = allLocal(ok);
allGlobal = allGlobal(ok);
allDiff = allDiff(ok);
allGroup = allGroup(ok);
if chooseLargestPositive
    [siteDiff, ix] = max(allDiff);
else
    [~, ix] = max(abs(allDiff));
    siteDiff = allDiff(ix);
end
siteLocal = allLocal(ix);
siteGlobal = allGlobal(ix);
groupName = allGroup(ix);
end

function [groupName, siteDiff] = lookup_site_in_selection(SEL, siteLocal)
groupName = "";
siteDiff = NaN;
for g = 1:numel(SEL.Groups)
    G = SEL.Groups(g);
    ix = find(G.siteRows(:) == siteLocal, 1, 'first');
    if ~isempty(ix)
        groupName = string(G.name);
        siteDiff = G.siteDiff(ix);
        return;
    end
end
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

function make_single_site_figure(OUT)
cPref = [0.82 0.22 0.16];
cNon = [0.30 0.30 0.30];
fig = figure('Color', 'w', 'Name', sprintf('IT prestim site %d', OUT.siteGlobal));
useTiled = exist('tiledlayout', 'file') == 2;
if useTiled
    tiledlayout(2,2, 'TileSpacing', 'compact', 'Padding', 'compact');
end

% Panel 1: strip plot by condition
if useTiled, ax1 = nexttile; else, ax1 = subplot(2,2,1); end %#ok<LAXES>
hold(ax1, 'on');
xp = 1 + 0.12 * (rand(size(OUT.prefVals)) - 0.5);
xn = 2 + 0.12 * (rand(size(OUT.nonVals)) - 0.5);
scatter(ax1, xp, OUT.prefVals, 14, 'MarkerFaceColor', cPref, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.45);
scatter(ax1, xn, OUT.nonVals, 14, 'MarkerFaceColor', cNon, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.35);
scatter(ax1, xp(OUT.outPref), OUT.prefVals(OUT.outPref), 28, 'MarkerEdgeColor', [0.45 0 0], 'MarkerFaceColor', 'none', 'LineWidth', 1.0);
scatter(ax1, xn(OUT.outNon), OUT.nonVals(OUT.outNon), 28, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', 'none', 'LineWidth', 1.0);
plot(ax1, [0.88 1.12], mean(OUT.prefVals,'omitnan')*[1 1], '-', 'Color', cPref, 'LineWidth', 2.2);
plot(ax1, [1.88 2.12], mean(OUT.nonVals,'omitnan')*[1 1], '-', 'Color', cNon, 'LineWidth', 2.2);
plot(ax1, [0.90 1.10], median(OUT.prefVals,'omitnan')*[1 1], '--', 'Color', [0.45 0 0], 'LineWidth', 1.2);
plot(ax1, [1.90 2.10], median(OUT.nonVals,'omitnan')*[1 1], '--', 'Color', [0 0 0], 'LineWidth', 1.2);
set(ax1, 'XTick', [1 2], 'XTickLabel', {'tar-close','dis-close'});
xlim(ax1, [0.5 2.5]);
ylabel(ax1, 'Prestim trial mean (normMUA)');
title(ax1, sprintf('%s | site %d (global %d)', char(OUT.siteGroup), OUT.siteLocal, OUT.siteGlobal), 'Interpreter', 'none');
grid(ax1, 'on');

% Panel 2: histogram overlay
if useTiled, ax2 = nexttile; else, ax2 = subplot(2,2,2); end %#ok<LAXES>
hold(ax2, 'on');
edges = choose_edges([OUT.prefVals; OUT.nonVals], 30);
histogram(ax2, OUT.nonVals, edges, 'Normalization', 'probability', ...
    'FaceColor', cNon, 'FaceAlpha', 0.28, 'EdgeColor', 'none');
histogram(ax2, OUT.prefVals, edges, 'Normalization', 'probability', ...
    'FaceColor', cPref, 'FaceAlpha', 0.34, 'EdgeColor', 'none');
xline(ax2, mean(OUT.prefVals,'omitnan'), '-', 'Color', cPref, 'LineWidth', 2.0);
xline(ax2, mean(OUT.nonVals,'omitnan'), '-', 'Color', cNon, 'LineWidth', 2.0);
xlabel(ax2, 'Prestim trial mean (normMUA)');
ylabel(ax2, 'Trial fraction');
title(ax2, sprintf('N=%d vs %d', numel(OUT.prefVals), numel(OUT.nonVals)));
legend(ax2, {'distractor-closer','target-closer'}, 'Location', 'best', 'Box', 'off');
grid(ax2, 'on');

% Panel 3: trial order
if useTiled, ax3 = nexttile; else, ax3 = subplot(2,2,3); end %#ok<LAXES>
hold(ax3, 'on');
trialAxis = 1:numel(OUT.allUsedTrials);
scatter(ax3, trialAxis(~OUT.isPrefUsed), OUT.usedVals(~OUT.isPrefUsed), 12, ...
    'MarkerFaceColor', cNon, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.35);
scatter(ax3, trialAxis(OUT.isPrefUsed), OUT.usedVals(OUT.isPrefUsed), 12, ...
    'MarkerFaceColor', cPref, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.45);
outMaskUsed = false(size(OUT.allUsedTrials));
outMaskUsed(OUT.isPrefUsed) = OUT.outPref;
outMaskUsed(~OUT.isPrefUsed) = OUT.outNon;
scatter(ax3, trialAxis(outMaskUsed), OUT.usedVals(outMaskUsed), 30, ...
    'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', 'none', 'LineWidth', 1.0);
sessChange = find(diff(OUT.usedSess(:)) ~= 0);
for i = 1:numel(sessChange)
    xline(ax3, sessChange(i) + 0.5, ':', 'Color', [0.75 0.75 0.75]);
end
xlabel(ax3, 'Included trial order');
ylabel(ax3, 'Prestim trial mean (normMUA)');
title(ax3, 'Trial order (dotted lines = session changes)');
grid(ax3, 'on');

% Panel 4: text summary
if useTiled, ax4 = nexttile; else, ax4 = subplot(2,2,4); end %#ok<LAXES>
axis(ax4, 'off');
lines = [OUT.summaryLines(:); " "; OUT.outlierLines(:)];
text(ax4, 0.02, 0.98, strjoin(cellstr(lines), newline), ...
    'Units', 'normalized', 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top', 'FontName', 'Courier', 'FontSize', 8);

annotation(fig, 'textbox', [0.05 0.965 0.90 0.035], ...
    'String', sprintf(['Single-site IT prestim trial diagnostic (%s). Target-closer and distractor-closer are defined only by RF geometry; ' ...
    'open circles mark Tukey outliers within each condition.'], char(OUT.monkeySuffix)), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontWeight', 'bold', 'FontSize', 11);
end

function make_session_timecourse_figure(OUT)
tb = OUT.tb(:);
fig = figure('Color', 'w', 'Name', sprintf('IT session timecourse site %d', OUT.siteGlobal));
ax = axes('Parent', fig); hold(ax, 'on');

cols = [0.15 0.35 0.75; 0.80 0.25 0.15; 0.25 0.55 0.25; 0.55 0.25 0.55];
nSess = numel(OUT.sessionIds);
lineH = gobjects(nSess, 1);
for iSess = 1:nSess
    c = cols(1 + mod(iSess-1, size(cols,1)), :);
    m = OUT.sessionMeanTc(iSess, :).';
    s = OUT.sessionSemTc(iSess, :).';
    add_sem_patch_local(ax, tb, m, s, c, 0.16);
    lineH(iSess) = plot(ax, tb, m, '-', 'Color', c, 'LineWidth', 2.1);
end
line(ax, [0 0], ylim(ax), 'Color', [0 0 0], 'LineStyle', '-', 'LineWidth', 1.0);
xlabel(ax, 'Time from stimulus onset (ms)');
ylabel(ax, 'Mean response (normMUA)');
ttl = sprintf('Site %d (global %d): all-stimuli mean timecourse by session', OUT.siteLocal, OUT.siteGlobal);
title(ax, ttl, 'Interpreter', 'none');
leg = strings(nSess,1);
for iSess = 1:nSess
    leg(iSess) = sprintf('session %d (N=%d)', OUT.sessionIds(iSess), OUT.sessionN(iSess));
end
legend(ax, lineH, cellstr(leg), 'Location', 'best', 'Box', 'off');
grid(ax, 'on');

preMask = (tb >= OUT.P.preWindowMs(1)) & (tb < OUT.P.preWindowMs(2));
txt = strings(0,1);
for iSess = 1:nSess
    txt(end+1) = sprintf('sess %d pre mean %.4f', OUT.sessionIds(iSess), mean(OUT.sessionMeanTc(iSess, preMask), 'omitnan')); %#ok<AGROW>
end
annotation(fig, 'textbox', [0.12 0.92 0.76 0.05], ...
    'String', strjoin(cellstr(txt), '   |   '), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'FontSize', 10);
end

function lines = build_summary_lines(groupName, siteLocal, siteGlobal, siteDiff, prefVals, nonVals, outPref, outNon, meanDiff, meanDiffNoOut, pRank)
lines = strings(0,1);
lines(end+1) = sprintf('Group: %s', char(groupName));
lines(end+1) = sprintf('Site local/global: %d / %d', siteLocal, siteGlobal);
lines(end+1) = sprintf('Summary-site diff: %.4f', siteDiff);
lines(end+1) = sprintf('Trials tar/dis: %d / %d', numel(prefVals), numel(nonVals));
lines(end+1) = sprintf('Mean tar/dis: %.4f / %.4f', mean(prefVals,'omitnan'), mean(nonVals,'omitnan'));
lines(end+1) = sprintf('Median tar/dis: %.4f / %.4f', median(prefVals,'omitnan'), median(nonVals,'omitnan'));
lines(end+1) = sprintf('Trim20 tar/dis: %.4f / %.4f', trimmean_local(prefVals,20), trimmean_local(nonVals,20));
lines(end+1) = sprintf('Mean diff: %.4f', meanDiff);
lines(end+1) = sprintf('No-outlier diff: %.4f', meanDiffNoOut);
lines(end+1) = sprintf('Ranksum p: %.3g', pRank);
lines(end+1) = sprintf('Outliers tar/dis: %d / %d', nnz(outPref), nnz(outNon));
end

function lines = build_outlier_lines(prefTrials, prefVals, prefStim, prefSess, outPref, nonTrials, nonVals, nonStim, nonSess, outNon, nTop)
lines = strings(0,1);
lines(end+1) = "Top outliers:";
tmp = format_one_condition_lines('tar-close', prefTrials, prefVals, prefStim, prefSess, outPref, nTop);
lines = [lines; tmp(:)]; %#ok<AGROW>
tmp = format_one_condition_lines('dis-close', nonTrials, nonVals, nonStim, nonSess, outNon, nTop);
lines = [lines; tmp(:)]; %#ok<AGROW>
end

function lines = format_one_condition_lines(label, trialIdx, vals, stim, sess, outMask, nTop)
lines = strings(0,1);
vals = vals(:);
trialIdx = trialIdx(:);
stim = stim(:);
sess = sess(:);
med = median(vals, 'omitnan');
score = abs(vals - med);
score(~outMask(:)) = -inf;
[~, ord] = sort(score, 'descend');
ord = ord(1:min(nTop, nnz(isfinite(score) & score > -inf)));
if isempty(ord)
    lines(end+1) = sprintf('  %s: no flagged outliers', label);
    return;
end
lines(end+1) = sprintf('  %s:', label);
for i = 1:numel(ord)
    k = ord(i);
    lines(end+1) = sprintf('    trial %d | stim %d | sess %d | %.4f', ...
        trialIdx(k), stim(k), sess(k), vals(k));
end
end

function edges = choose_edges(x, nBins)
x = double(x(:));
x = x(isfinite(x));
if isempty(x)
    edges = linspace(-1, 1, nBins + 1);
    return;
end
if numel(unique(x)) < 2
    dx = max(abs(x(1)),1)*0.05 + 1e-3;
    edges = linspace(x(1)-dx, x(1)+dx, nBins + 1);
else
    edges = linspace(min(x), max(x), nBins + 1);
end
end

function [mask, fence] = tukey_outlier_mask(x)
x = double(x(:));
x = x(isfinite(x));
mask = false(size(x));
fence = [NaN NaN];
if numel(x) < 4
    return;
end
q = prctile(x, [25 75]);
iqrVal = q(2) - q(1);
lo = q(1) - 1.5 * iqrVal;
hi = q(2) + 1.5 * iqrVal;
mask = (x < lo) | (x > hi);
fence = [lo hi];
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
k = floor((pct/100) * n / 2);
if (2*k) >= n
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

function add_sem_patch_local(ax, t, m, s, col, faceAlpha)
t = t(:);
m = m(:);
s = s(:);
ok = isfinite(t) & isfinite(m) & isfinite(s);
if nnz(ok) < 3
    return;
end
tt = t(ok);
lo = m(ok) - s(ok);
hi = m(ok) + s(ok);
patch(ax, [tt; flipud(tt)], [lo; flipud(hi)], col, ...
    'FaceAlpha', faceAlpha, 'EdgeColor', 'none', ...
    'HandleVisibility', 'off');
end
