function OUT = NoiseCorr_V1V4_PrecedingTrialType(Monkey, Puser)
% NOISECORR_V1V4_PRECEDINGTRIALTYPE
% Diagnose whether current TT and DD trials occur in different lag-1 trial
% contexts for the V1-V4 pair set used in the RF-distance analysis.

if nargin < 1 || isempty(Monkey)
    Monkey = 1;
end
if nargin < 2 || isempty(Puser)
    Puser = struct();
end

P = struct();
P.cacheTag = "all";
P.rfMinPx = 50;
P.minTrialsPerCond = 15;
P.sessions = [1 2];
P.onlyCorrect = true;
P.correctCol = 9;
P.correctVal = 1;
P.stimCol = 1;
P.dayCol = 11;
P.useCache = true;
P.saveResult = true;
P.plotFigure = true;
P.verbose = true;

if ~isempty(Puser)
    fn = fieldnames(Puser);
    for i = 1:numel(fn)
        P.(fn{i}) = Puser.(fn{i});
    end
end

cfg = config();
[monkeySuffix, monkeyFolder] = local_monkey_info(Monkey);
outPath = fullfile(cfg.resultsDir, sprintf('NoiseCorr_V1V4_PrecedingTrialType_%s_%s.mat', char(monkeySuffix), char(P.cacheTag)));
cacheParams = local_cache_params(Monkey, P);

if P.useCache && exist(outPath, 'file') == 2
    S = load(outPath, 'OUT');
    if isfield(S, 'OUT') && local_cache_matches(S.OUT, cacheParams) && ...
            session_exclusion_cache_matches(S.OUT, monkeySuffix)
        OUT = S.OUT;
        OUT.P = P;
        if P.verbose
            fprintf('Loaded cached V1-V4 preceding-trial diagnostic from %s\n', outPath);
        end
        if P.plotFigure
            local_plot_summary(OUT);
        end
        return;
    end
end

Main = NoiseCorr_V1V4_RFDistance_TargetDistractor(Monkey, struct('windowIdx', [1 2 3], ...
    'plotFigure', false, 'verbose', false, 'useCache', true));

trialInfoPath = fullfile(cfg.dataRoot, monkeyFolder, 'ObjAtt_lines_MUA_trials.mat');
assert(exist(trialInfoPath, 'file') == 2, 'Missing %s.', trialInfoPath);
m2 = matfile(trialInfoPath);
ALLMAT = double(m2.ALLMAT);

stimAll = double(ALLMAT(:, P.stimCol));
dayAll = double(ALLMAT(:, P.dayCol));
validStim = isfinite(stimAll) & stimAll >= 1 & stimAll <= 384 & (floor(stimAll) == stimAll);
validSess = ismember(dayAll, P.sessions(:));
validRaw = validStim & validSess;
currMask = validRaw;
if P.onlyCorrect
    currMask = currMask & (double(ALLMAT(:, P.correctCol)) == P.correctVal);
end

prevPresented = nan(size(stimAll));
prevIncluded = nan(size(stimAll));
for k = 1:numel(P.sessions)
    sess = P.sessions(k);
    rowsRaw = find(dayAll == sess & validRaw);
    if numel(rowsRaw) >= 2
        prevPresented(rowsRaw(2:end)) = rowsRaw(1:end-1);
    end
    rowsIncl = find(dayAll == sess & currMask);
    if numel(rowsIncl) >= 2
        prevIncluded(rowsIncl(2:end)) = rowsIncl(1:end-1);
    end
end

trialRowsCurr = find(currMask);
stimCurr = stimAll(trialRowsCurr);
prevPresentedRowsCurr = prevPresented(trialRowsCurr);
prevIncludedRowsCurr = prevIncluded(trialRowsCurr);

siteMapV1 = nan(max(double(Main.SiteTableV1.siteIdx)), 1);
siteMapV1(double(Main.SiteTableV1.siteIdx)) = 1:height(Main.SiteTableV1);
siteMapV4 = nan(max(double(Main.SiteTableV4.siteIdx)), 1);
siteMapV4(double(Main.SiteTableV4.siteIdx)) = 1:height(Main.SiteTableV4);

pairMask = isfinite(Main.rfDistancePx) & (Main.rfDistancePx > P.rfMinPx);
pairIdxUse = find(pairMask);
nPairsUse = numel(pairIdxUse);

currCondNames = ["T-T"; "D-D"];
prevCondNames = ["T-T"; "D-D"; "T-D"; "B-B"; "B-T"; "B-D"; "other"; "start"];
nPrev = numel(prevCondNames);

pairCurrentCount = zeros(nPairsUse, 2, 'uint16');
pairPrevFracPresented = nan(nPairsUse, 2, nPrev);
pairPrevFracIncluded = nan(nPairsUse, 2, nPrev);

for kp = 1:nPairsUse
    p = pairIdxUse(kp);
    rowV1 = siteMapV1(double(Main.PairTable.site1Idx(p)));
    rowV4 = siteMapV4(double(Main.PairTable.site2Idx(p)));
    assert(isfinite(rowV1) && isfinite(rowV4), 'Could not map pair %d to V1/V4 rows.', p);

    condByStim = local_pair_cond_by_stim(Main.assignCodeV1(rowV1, :), Main.assignCodeV4(rowV4, :));
    currCond = condByStim(stimCurr);

    prevCodePresented = local_prev_codes(condByStim, stimAll, prevPresentedRowsCurr);
    prevCodeIncluded = local_prev_codes(condByStim, stimAll, prevIncludedRowsCurr);

    for c = 1:2
        sel = (currCond == c);
        nNow = nnz(sel);
        pairCurrentCount(kp, c) = uint16(nNow);
        if nNow < P.minTrialsPerCond
            continue;
        end
        pairPrevFracPresented(kp, c, :) = local_fraction_by_code(prevCodePresented(sel), nPrev);
        pairPrevFracIncluded(kp, c, :) = local_fraction_by_code(prevCodeIncluded(sel), nPrev);
    end

    if P.verbose && (mod(kp, 2000) == 0 || kp == nPairsUse)
        fprintf('  Preceding-trial diagnostic pairs: %d / %d\n', kp, nPairsUse);
    end
end

validPair = squeeze(pairCurrentCount(:, 1) >= P.minTrialsPerCond & pairCurrentCount(:, 2) >= P.minTrialsPerCond);
nValidPairs = nnz(validPair);

[meanPresented, semPresented] = local_mean_sem(pairPrevFracPresented(validPair, :, :));
[meanIncluded, semIncluded] = local_mean_sem(pairPrevFracIncluded(validPair, :, :));
[meanDiffPresented, semDiffPresented] = local_mean_sem_diff(pairPrevFracPresented(validPair, :, :));
[meanDiffIncluded, semDiffIncluded] = local_mean_sem_diff(pairPrevFracIncluded(validPair, :, :));

OUT = struct();
OUT.P = P;
OUT.cacheParams = cacheParams;
OUT.siteSessionExclusions = site_session_exclusions(monkeySuffix);
OUT.monkeySuffix = monkeySuffix;
OUT.mainPath = fullfile(cfg.resultsDir, sprintf('NoiseCorr_V1V4_RFDistance_TargetDistractor_%s_%s.mat', char(monkeySuffix), char(P.cacheTag)));
OUT.rfMinPx = P.rfMinPx;
OUT.pairIdxUse = pairIdxUse;
OUT.validPair = validPair;
OUT.nValidPairs = nValidPairs;
OUT.currCondNames = currCondNames;
OUT.prevCondNames = prevCondNames;
OUT.pairCurrentCount = pairCurrentCount;
OUT.pairPrevFracPresented = pairPrevFracPresented;
OUT.pairPrevFracIncluded = pairPrevFracIncluded;
OUT.meanPresented = meanPresented;
OUT.semPresented = semPresented;
OUT.meanIncluded = meanIncluded;
OUT.semIncluded = semIncluded;
OUT.meanDiffPresented = meanDiffPresented;
OUT.semDiffPresented = semDiffPresented;
OUT.meanDiffIncluded = meanDiffIncluded;
OUT.semDiffIncluded = semDiffIncluded;

if P.saveResult
    save(outPath, 'OUT', '-v7.3');
    if P.verbose
        fprintf('Saved V1-V4 preceding-trial diagnostic to %s\n', outPath);
    end
end

if P.plotFigure
    local_plot_summary(OUT);
end
end

function condByStim = local_pair_cond_by_stim(code1, code2)
condByStim = zeros(384, 1, 'uint8');
tt = (code1 == 1) & (code2 == 1);
dd = (code1 == 2) & (code2 == 2);
td = ((code1 == 1) & (code2 == 2)) | ((code1 == 2) & (code2 == 1));
bb = (code1 == 3) & (code2 == 3);
bt = ((code1 == 3) & (code2 == 1)) | ((code1 == 1) & (code2 == 3));
bd = ((code1 == 3) & (code2 == 2)) | ((code1 == 2) & (code2 == 3));
condByStim(tt) = 1;
condByStim(dd) = 2;
condByStim(td) = 3;
condByStim(bb) = 4;
condByStim(bt) = 5;
condByStim(bd) = 6;
end

function prevCode = local_prev_codes(condByStim, stimAll, prevRows)
prevCode = uint8(8 * ones(size(prevRows))); % 8 = start
hasPrev = isfinite(prevRows) & prevRows >= 1;
if ~any(hasPrev)
    return;
end
prevStim = stimAll(prevRows(hasPrev));
validPrevStim = isfinite(prevStim) & prevStim >= 1 & prevStim <= 384 & (floor(prevStim) == prevStim);
codes = uint8(7 * ones(size(prevStim))); % 7 = other
codes(validPrevStim) = condByStim(prevStim(validPrevStim));
codes(codes == 0) = 7;
prevCode(hasPrev) = codes;
end

function frac = local_fraction_by_code(codeVec, nPrev)
frac = nan(1, 1, nPrev);
codeVec = double(codeVec(:));
if isempty(codeVec)
    return;
end
counts = histcounts(codeVec, 0.5:1:(nPrev + 0.5));
frac(1, 1, :) = counts ./ sum(counts);
end

function [mu, sem] = local_mean_sem(X)
mu = squeeze(mean(X, 1, 'omitnan'));
sem = squeeze(std(X, 0, 1, 'omitnan')) ./ sqrt(size(X, 1));
end

function [mu, sem] = local_mean_sem_diff(X)
delta = squeeze(X(:, 1, :) - X(:, 2, :));
mu = mean(delta, 1, 'omitnan')';
sem = std(delta, 0, 1, 'omitnan')' ./ sqrt(size(delta, 1));
end

function [monkeySuffix, monkeyFolder] = local_monkey_info(monkeyId)
switch monkeyId
    case 1
        monkeySuffix = "N";
        monkeyFolder = 'Mr Nilson';
    case 2
        monkeySuffix = "F";
        monkeyFolder = 'Figaro';
    otherwise
        error('Monkey must be 1 (Nilson) or 2 (Figaro).');
end
end

function S = local_cache_params(Monkey, P)
S = struct();
S.cacheVersion = 1;
S.Monkey = double(Monkey);
S.cacheTag = string(P.cacheTag);
S.rfMinPx = double(P.rfMinPx);
S.minTrialsPerCond = double(P.minTrialsPerCond);
S.sessions = double(P.sessions(:)');
S.onlyCorrect = logical(P.onlyCorrect);
S.correctCol = double(P.correctCol);
S.correctVal = double(P.correctVal);
S.stimCol = double(P.stimCol);
S.dayCol = double(P.dayCol);
end

function tf = local_cache_matches(OUT, cacheParams)
tf = isstruct(OUT) && isfield(OUT, 'cacheParams') && isequaln(OUT.cacheParams, cacheParams);
end

function local_plot_summary(OUT)
figure('Color', 'w', 'Name', sprintf('V1-V4 preceding-trial diagnostic (%s)', char(OUT.monkeySuffix)));
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

ax1 = nexttile;
imagesc(ax1, OUT.meanPresented);
colorbar(ax1);
set(ax1, 'XTick', 1:numel(OUT.prevCondNames), 'XTickLabel', cellstr(OUT.prevCondNames), ...
    'YTick', 1:numel(OUT.currCondNames), 'YTickLabel', cellstr(OUT.currCondNames));
xtickangle(ax1, 30);
title(ax1, sprintf('Prev presented trial | RF > %d px | n=%d pairs', OUT.rfMinPx, OUT.nValidPairs));
ylabel(ax1, 'Current trial type');

ax2 = nexttile;
hold(ax2, 'on');
x = 1:numel(OUT.prevCondNames);
bar(ax2, x, OUT.meanDiffPresented, 'FaceColor', [0.25 0.45 0.80], 'EdgeColor', 'none');
errorbar(ax2, x, OUT.meanDiffPresented, OUT.semDiffPresented, 'k.', 'LineWidth', 1.0, 'HandleVisibility', 'off');
plot(ax2, x, zeros(size(x)), 'k:', 'HandleVisibility', 'off');
set(ax2, 'XTick', x, 'XTickLabel', cellstr(OUT.prevCondNames));
xtickangle(ax2, 30);
ylabel(ax2, 'TT minus DD fraction');
title(ax2, 'Prev presented trial difference');
grid(ax2, 'on');

ax3 = nexttile;
imagesc(ax3, OUT.meanIncluded);
colorbar(ax3);
set(ax3, 'XTick', 1:numel(OUT.prevCondNames), 'XTickLabel', cellstr(OUT.prevCondNames), ...
    'YTick', 1:numel(OUT.currCondNames), 'YTickLabel', cellstr(OUT.currCondNames));
xtickangle(ax3, 30);
title(ax3, 'Prev included-correct trial');
ylabel(ax3, 'Current trial type');

ax4 = nexttile;
hold(ax4, 'on');
x = 1:numel(OUT.prevCondNames);
bar(ax4, x, OUT.meanDiffIncluded, 'FaceColor', [0.80 0.45 0.20], 'EdgeColor', 'none');
errorbar(ax4, x, OUT.meanDiffIncluded, OUT.semDiffIncluded, 'k.', 'LineWidth', 1.0, 'HandleVisibility', 'off');
plot(ax4, x, zeros(size(x)), 'k:', 'HandleVisibility', 'off');
set(ax4, 'XTick', x, 'XTickLabel', cellstr(OUT.prevCondNames));
xtickangle(ax4, 30);
ylabel(ax4, 'TT minus DD fraction');
title(ax4, 'Prev included-correct difference');
grid(ax4, 'on');
end
