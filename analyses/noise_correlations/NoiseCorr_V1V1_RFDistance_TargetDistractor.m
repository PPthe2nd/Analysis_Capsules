function OUT = NoiseCorr_V1V1_RFDistance_TargetDistractor(Monkey, Puser)
% NOISECORR_V1V1_RFDISTANCE_TARGETDISTRACTOR
% Plot V1-V1 noise correlation versus RF distance for target/distractor/
% background pair conditions. RF distance uses the V1 RF center
% coordinates in pixels from Tall_V1; condition labels also use Tall_V1.

if nargin < 1 || isempty(Monkey)
    Monkey = 1;
end
if nargin < 2 || isempty(Puser)
    Puser = struct();
end

P = struct();
P.cacheTag = "all";
P.windowIdx = 3;              % default = 300-500 ms; can also be [1 2 3]
P.binWidthPx = 10;
P.minTrials = 15;
P.summaryMinRfPx = 50;
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
[monkeySuffix, monkeyFolder] = local_monkey_info(Monkey); %#ok<ASGLU>
outPath = fullfile(cfg.resultsDir, sprintf('NoiseCorr_V1V1_RFDistance_TargetDistractor_%s_%s.mat', char(monkeySuffix), char(P.cacheTag)));
cacheParams = local_cache_params(Monkey, P);

if P.useCache && exist(outPath, 'file') == 2
    S = load(outPath, 'OUT');
    if isfield(S, 'OUT') && local_cache_matches(S.OUT, cacheParams) && ...
            session_exclusion_cache_matches(S.OUT, monkeySuffix)
        OUT = S.OUT;
        OUT.P = P;
        if P.verbose
            fprintf('Loaded cached V1-V1 RF-distance noise-correlation analysis from %s\n', outPath);
        end
        if P.plotFigure
            local_plot_summary(OUT);
            local_plot_difference_summary(OUT);
        end
        return;
    end
end

masterPath = fullfile(cfg.resultsDir, sprintf('NoiseCorr_Cache_%s_%s.mat', char(monkeySuffix), char(P.cacheTag)));
assert(exist(masterPath, 'file') == 2, ...
    'Missing %s. Run Build_NoiseCorr_Cache first.', masterPath);
S = load(masterPath, 'OUT');
assert(isfield(S, 'OUT') && isstruct(S.OUT), '%s must contain OUT.', masterPath);
Master = S.OUT;
assert(isfield(Master, 'SiteTable') && isfield(Master, 'PairClassTable'), ...
    'Master noise-correlation cache is missing required fields.');
windowIdx = unique(double(P.windowIdx(:)'));
assert(~isempty(windowIdx) && all(isfinite(windowIdx)) && ...
    all(windowIdx >= 1) && all(windowIdx <= size(Master.timeWindows, 1)), ...
    'windowIdx must contain values between 1 and %d.', size(Master.timeWindows, 1));
nWinSel = numel(windowIdx);

if P.verbose
    fprintf('Loaded master noise-correlation cache from %s\n', masterPath);
end

pairBlock = local_load_v1v1_pair_block(Master, cfg, monkeySuffix, P.cacheTag);
PairTable = pairBlock.PairTable;
sumXY = pairBlock.sumXY;

siteRowsV1 = find(Master.SiteTable.area == "V1");
assert(~isempty(siteRowsV1), 'No V1 sites found in the master cache.');
SiteTableV1 = Master.SiteTable(siteRowsV1, :);
nV1 = height(SiteTableV1);
nStim = size(Master.nTrials, 1);
nSess = size(Master.nTrials, 2);

tallPath = fullfile(cfg.matDir, sprintf('Tall_V1_lines_%s.mat', char(monkeySuffix)));
assert(exist(tallPath, 'file') == 2, 'Missing %s.', tallPath);
Sgeo = load(tallPath, 'Tall_V1');
Tall_V1 = local_sort_tall(Sgeo.Tall_V1, 'Tall_V1');
assignCode = local_build_assign_code(Tall_V1, SiteTableV1.areaLocalSite);
xPxAll = double(Tall_V1(1).T.x_px(1:512));
yPxAll = double(Tall_V1(1).T.y_px(1:512));
xPx = xPxAll(double(SiteTableV1.areaLocalSite));
yPx = yPxAll(double(SiteTableV1.areaLocalSite));

siteIdxToV1Row = nan(height(Master.SiteTable), 1);
siteIdxToV1Row(siteRowsV1) = 1:nV1;
pairSite1V1 = siteIdxToV1Row(double(PairTable.site1Idx));
pairSite2V1 = siteIdxToV1Row(double(PairTable.site2Idx));
assert(all(isfinite(pairSite1V1)) && all(isfinite(pairSite2V1)), ...
    'V1-V1 pair block contains non-V1 sites.');

nPairs = height(PairTable);
rfDistancePx = hypot(xPx(pairSite1V1) - xPx(pairSite2V1), ...
                     yPx(pairSite1V1) - yPx(pairSite2V1));

condNames = ["T-T"; "D-D"; "T-D"; "B-B"; "B-T"; "B-D"];
nCond = numel(condNames);
pairCorrCond = nan(nPairs, nCond, nWinSel);
pairTrialCond = zeros(nPairs, nCond, nWinSel);
pairStimCond = zeros(nPairs, nCond);

sumX2Master = double(Master.sumX2(:,:, :, windowIdx));
nTrialsSiteStimSess = double(Master.nTrialsSiteStimSess);

for p = 1:nPairs
    i1 = pairSite1V1(p);
    i2 = pairSite2V1(p);
    code1 = assignCode(i1, :);
    code2 = assignCode(i2, :);

    condMasks = false(nCond, nStim);
    condMasks(1, :) = (code1 == 1) & (code2 == 1); % T-T
    condMasks(2, :) = (code1 == 2) & (code2 == 2); % D-D
    condMasks(3, :) = ((code1 == 1) & (code2 == 2)) | ((code1 == 2) & (code2 == 1)); % T-D pooled
    condMasks(4, :) = (code1 == 3) & (code2 == 3); % B-B
    condMasks(5, :) = ((code1 == 3) & (code2 == 1)) | ((code1 == 1) & (code2 == 3)); % B-T pooled
    condMasks(6, :) = ((code1 == 3) & (code2 == 2)) | ((code1 == 2) & (code2 == 3)); % B-D pooled

    site1MasterIdx = double(PairTable.site1Idx(p));
    site2MasterIdx = double(PairTable.site2Idx(p));
    n1 = squeeze(nTrialsSiteStimSess(site1MasterIdx, :, :));
    n2 = squeeze(nTrialsSiteStimSess(site2MasterIdx, :, :));
    nPairStimSess = min(n1, n2);

    for c = 1:nCond
        stimMask = condMasks(c, :);
        pairStimCond(p, c) = nnz(stimMask);
        if ~any(stimMask)
            continue;
        end
        nUse = sum(nPairStimSess(stimMask, :), 'all', 'omitnan');
        if nUse < P.minTrials
            continue;
        end

        for wSel = 1:nWinSel
            wIdx = windowIdx(wSel);
            cxyStimSess = squeeze(double(sumXY(p, :, :, wIdx)));
            sx2_1 = squeeze(sumX2Master(site1MasterIdx, :, :, wSel));
            sx2_2 = squeeze(sumX2Master(site2MasterIdx, :, :, wSel));
            pairTrialCond(p, c, wSel) = nUse;
            num = sum(cxyStimSess(stimMask, :), 'all', 'omitnan');
            den1 = sum(sx2_1(stimMask, :), 'all', 'omitnan');
            den2 = sum(sx2_2(stimMask, :), 'all', 'omitnan');
            denom = sqrt(den1 * den2);
            if isfinite(denom) && denom > 0
                pairCorrCond(p, c, wSel) = num / denom;
            end
        end
    end
end

maxDist = max(rfDistancePx(isfinite(rfDistancePx)));
binEdges = 0:P.binWidthPx:(ceil(maxDist / P.binWidthPx) * P.binWidthPx + P.binWidthPx);
if numel(binEdges) < 2
    binEdges = [0 P.binWidthPx];
end
binCenters = 0.5 * (binEdges(1:end-1) + binEdges(2:end));
nBins = numel(binCenters);

meanByBin = nan(nBins, nCond, nWinSel);
semByBin = nan(nBins, nCond, nWinSel);
nPairsByBin = zeros(nBins, nCond, nWinSel);
meanAboveMinRf = nan(nCond, nWinSel);
semAboveMinRf = nan(nCond, nWinSel);
nPairsAboveMinRf = zeros(nCond, nWinSel);
for wSel = 1:nWinSel
    for c = 1:nCond
        valid = isfinite(pairCorrCond(:, c, wSel)) & isfinite(rfDistancePx);
        binIdx = discretize(rfDistancePx(valid), binEdges);
        vals = pairCorrCond(valid, c, wSel);
        for b = 1:nBins
            xb = vals(binIdx == b);
            nPairsByBin(b, c, wSel) = numel(xb);
            if isempty(xb)
                continue;
            end
            meanByBin(b, c, wSel) = mean(xb, 'omitnan');
            if numel(xb) >= 2
                semByBin(b, c, wSel) = std(xb, 0, 'omitnan') / sqrt(numel(xb));
            else
                semByBin(b, c, wSel) = 0;
            end
        end

        validFar = valid & (rfDistancePx > P.summaryMinRfPx);
        valsFar = pairCorrCond(validFar, c, wSel);
        nPairsAboveMinRf(c, wSel) = numel(valsFar);
        if isempty(valsFar)
            continue;
        end
        meanAboveMinRf(c, wSel) = mean(valsFar, 'omitnan');
        if numel(valsFar) >= 2
            semAboveMinRf(c, wSel) = std(valsFar, 0, 'omitnan') / sqrt(numel(valsFar));
        else
            semAboveMinRf(c, wSel) = 0;
        end
    end
end

prePos = find(windowIdx == 1, 1, 'first');
latePos = find(windowIdx == 3, 1, 'first');
meanLateMinusPre = nan(nCond, 1);
semLateMinusPre = nan(nCond, 1);
nPairsLateMinusPre = zeros(nCond, 1);
if ~isempty(prePos) && ~isempty(latePos)
    for c = 1:nCond
        validDiff = isfinite(rfDistancePx) & (rfDistancePx > P.summaryMinRfPx) & ...
            isfinite(pairCorrCond(:, c, prePos)) & isfinite(pairCorrCond(:, c, latePos));
        delta = pairCorrCond(validDiff, c, latePos) - pairCorrCond(validDiff, c, prePos);
        nPairsLateMinusPre(c) = numel(delta);
        if isempty(delta)
            continue;
        end
        meanLateMinusPre(c) = mean(delta, 'omitnan');
        if numel(delta) >= 2
            semLateMinusPre(c) = std(delta, 0, 'omitnan') / sqrt(numel(delta));
        else
            semLateMinusPre(c) = 0;
        end
    end
end

OUT = struct();
OUT.P = P;
OUT.cacheParams = cacheParams;
OUT.siteSessionExclusions = site_session_exclusions(monkeySuffix);
OUT.masterPath = masterPath;
OUT.pairFile = pairBlock.pairFile;
OUT.monkeySuffix = monkeySuffix;
OUT.windowIdx = windowIdx(:);
OUT.timeWindowsSelected = Master.timeWindows(windowIdx, :);
if numel(windowIdx) == 1
    OUT.timeWindow = Master.timeWindows(windowIdx, :);
end
OUT.condNames = condNames;
OUT.SiteTableV1 = SiteTableV1;
OUT.PairTable = PairTable;
OUT.rfDistancePx = rfDistancePx;
OUT.pairCorrCond = pairCorrCond;
OUT.pairTrialCond = pairTrialCond;
OUT.pairStimCond = pairStimCond;
OUT.assignCode = assignCode;
OUT.binEdges = binEdges;
OUT.binCenters = binCenters;
OUT.meanByBin = meanByBin;
OUT.semByBin = semByBin;
OUT.nPairsByBin = nPairsByBin;
OUT.meanAboveMinRf = meanAboveMinRf;
OUT.semAboveMinRf = semAboveMinRf;
OUT.nPairsAboveMinRf = nPairsAboveMinRf;
OUT.meanLateMinusPre = meanLateMinusPre;
OUT.semLateMinusPre = semLateMinusPre;
OUT.nPairsLateMinusPre = nPairsLateMinusPre;
OUT.preWindowPos = prePos;
OUT.lateWindowPos = latePos;

if P.saveResult
    save(outPath, 'OUT', '-v7.3');
    if P.verbose
        fprintf('Saved V1-V1 RF-distance analysis to %s\n', outPath);
    end
end

if P.plotFigure
    local_plot_summary(OUT);
    local_plot_difference_summary(OUT);
end
end

function pairBlock = local_load_v1v1_pair_block(Master, cfg, monkeySuffix, cacheTag)
pairBlock = struct();
pairBlock.pairFile = "";

if isfield(Master, 'PairClassTable') && istable(Master.PairClassTable)
    idx = find(Master.PairClassTable.pairClass == "V1-V1", 1, 'first');
    assert(~isempty(idx), 'Master cache does not contain V1-V1 pair metadata.');
    if ismember('pairFile', Master.PairClassTable.Properties.VariableNames)
        pairFile = string(Master.PairClassTable.pairFile(idx));
        assert(exist(char(pairFile), 'file') == 2, ...
            'Missing %s. Rebuild the noise-correlation cache.', char(pairFile));
        S = load(char(pairFile), 'PairOut');
        assert(isfield(S, 'PairOut') && isstruct(S.PairOut), '%s must contain PairOut.', char(pairFile));
        pairBlock.PairTable = S.PairOut.PairTable;
        pairBlock.sumXY = S.PairOut.sumXY;
        pairBlock.pairFile = pairFile;
        return;
    end
end

% Backward compatibility with older monolithic master cache.
assert(isfield(Master, 'sumXY') && isfield(Master, 'PairTable') && isfield(Master, 'PairClassTable'), ...
    'Master cache does not contain V1-V1 pair data and no pair file is available.');
idx = find(Master.PairClassTable.pairClass == "V1-V1", 1, 'first');
assert(~isempty(idx), 'Master cache does not contain V1-V1 pair metadata.');
s = Master.PairClassTable.startIdx(idx);
e = Master.PairClassTable.endIdx(idx);
assert(isfinite(s) && isfinite(e), 'Invalid V1-V1 pair range in master cache.');
pairBlock.PairTable = Master.PairTable(s:e, :);
pairBlock.sumXY = Master.sumXY(s:e, :, :, :);
pairBlock.pairFile = fullfile(cfg.resultsDir, sprintf('NoiseCorr_Cache_%s_%s.mat', char(monkeySuffix), char(cacheTag)));
end

function assignCode = local_build_assign_code(Tall_V1, areaLocalSite)
nSites = numel(areaLocalSite);
nStim = numel(Tall_V1);
assignCode = zeros(nSites, nStim, 'uint8');
for stim = 1:nStim
    T = Tall_V1(stim).T;
    assign = string(T.assignment(areaLocalSite));
    overlap = false(nSites, 1);
    if ismember('overlap', T.Properties.VariableNames)
        overlap = logical(T.overlap(areaLocalSite));
    end
    code = zeros(nSites, 1, 'uint8');
    code(assign == "target" & ~overlap) = 1;
    code(assign == "distractor" & ~overlap) = 2;
    code(assign == "background" & ~overlap) = 3;
    assignCode(:, stim) = code;
end
end

function TallSorted = local_sort_tall(Tall, tallFieldName)
stimNums = arrayfun(@(x) x.stimNum, Tall(:));
[stimNumsSorted, ordStim] = sort(stimNums(:));
assert(all(stimNumsSorted(:).' == 1:numel(Tall)), ...
    '%s.stimNum must cover 1..%d exactly.', tallFieldName, numel(Tall));
TallSorted = Tall(ordStim);
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
S.cacheVersion = 6;
S.Monkey = double(Monkey);
S.cacheTag = string(P.cacheTag);
S.windowIdx = double(unique(P.windowIdx(:)'));
S.binWidthPx = double(P.binWidthPx);
S.minTrials = double(P.minTrials);
S.summaryMinRfPx = double(P.summaryMinRfPx);
end

function local_plot_difference_summary(OUT)
if isempty(OUT.preWindowPos) || isempty(OUT.lateWindowPos)
    return;
end

colors = [ ...
    0.15 0.45 0.85; ... % T-T
    0.85 0.40 0.15; ... % D-D
    0.20 0.60 0.30; ... % T-D
    0.15 0.15 0.15; ... % B-B
    0.45 0.45 0.45; ... % B-T
    0.72 0.72 0.72];    % B-D

twPre = OUT.timeWindowsSelected(OUT.preWindowPos, :);
twLate = OUT.timeWindowsSelected(OUT.lateWindowPos, :);

figure('Color', 'w', 'Name', sprintf('V1-V1 late-minus-pre summary (%s)', char(OUT.monkeySuffix)));
ax = axes();
hold(ax, 'on');
xBar = 1:numel(OUT.condNames);
yBar = OUT.meanLateMinusPre(:);
eBar = OUT.semLateMinusPre(:);
b = bar(ax, xBar, yBar, 'FaceColor', 'flat', 'EdgeColor', 'none');
b.CData = colors;
good = isfinite(yBar);
if any(good)
    errorbar(ax, xBar(good), yBar(good), eBar(good), 'k.', 'LineWidth', 1.0, ...
        'HandleVisibility', 'off');
end
plot(ax, xBar, zeros(size(xBar)), 'k:', 'HandleVisibility', 'off');
for c = 1:numel(OUT.condNames)
    if isfinite(yBar(c))
        text(ax, xBar(c), yBar(c), sprintf('  n=%d', OUT.nPairsLateMinusPre(c)), ...
            'Rotation', 90, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', ...
            'FontSize', 8, 'Color', [0.2 0.2 0.2]);
    end
end
set(ax, 'XTick', xBar, 'XTickLabel', cellstr(OUT.condNames));
xtickangle(ax, 30);
ylabel(ax, '\Delta noise correlation (late - pre)');
title(ax, sprintf('RF distance > %d px | %d-%d minus %d-%d ms', ...
    OUT.P.summaryMinRfPx, twLate(1), twLate(2), twPre(1), twPre(2)));
grid(ax, 'on');
end

function tf = local_cache_matches(OUT, cacheParams)
tf = isstruct(OUT) && isfield(OUT, 'cacheParams') && isequaln(OUT.cacheParams, cacheParams);
end

function local_plot_summary(OUT)
colors = [ ...
    0.15 0.45 0.85; ... % T-T
    0.85 0.40 0.15; ... % D-D
    0.20 0.60 0.30; ... % T-D
    0.15 0.15 0.15; ... % B-B
    0.45 0.45 0.45; ... % B-T
    0.72 0.72 0.72];    % B-D
figure('Color', 'w', 'Name', sprintf('V1-V1 noise corr vs RF distance (%s)', char(OUT.monkeySuffix)));
tiledlayout(numel(OUT.windowIdx),2,'Padding','compact','TileSpacing','compact');
lineCondIdx = [1 3 2];

for wSel = 1:numel(OUT.windowIdx)
    tw = OUT.timeWindowsSelected(wSel, :);

    ax1 = nexttile;
    hold(ax1, 'on');
    for c = lineCondIdx
        x = OUT.binCenters;
        y = OUT.meanByBin(:, c, wSel);
        e = OUT.semByBin(:, c, wSel);
        good = isfinite(y);
        if any(good)
            patch(ax1, [x(good) fliplr(x(good))], ...
                [transpose(y(good)-e(good)) fliplr(transpose(y(good)+e(good)))], ...
                colors(c,:), 'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            plot(ax1, x(good), y(good), '-', 'Color', colors(c,:), 'LineWidth', 1.8, ...
                'DisplayName', sprintf('%s', OUT.condNames(c)));
        end
    end
    xlabel(ax1, 'RF distance (px)');
    ylabel(ax1, 'Mean noise correlation');
    title(ax1, sprintf('V1-V1 noise correlation | %d-%d ms | min trials = %d', ...
        tw(1), tw(2), OUT.P.minTrials));
    grid(ax1, 'on');
    legend(ax1, 'Location', 'best', 'Box', 'off');

    ax2 = nexttile;
    hold(ax2, 'on');
    xBar = 1:numel(OUT.condNames);
    yBar = OUT.meanAboveMinRf(:, wSel);
    eBar = OUT.semAboveMinRf(:, wSel);
    b = bar(ax2, xBar, yBar, 'FaceColor', 'flat', 'EdgeColor', 'none');
    b.CData = colors;
    good = isfinite(yBar);
    if any(good)
        errorbar(ax2, xBar(good), yBar(good), eBar(good), 'k.', 'LineWidth', 1.0, ...
            'HandleVisibility', 'off');
    end
    for c = 1:numel(OUT.condNames)
        if isfinite(yBar(c))
            text(ax2, xBar(c), yBar(c), sprintf('  n=%d', OUT.nPairsAboveMinRf(c, wSel)), ...
                'Rotation', 90, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', ...
                'FontSize', 8, 'Color', [0.2 0.2 0.2]);
        end
    end
    set(ax2, 'XTick', xBar, 'XTickLabel', cellstr(OUT.condNames));
    xtickangle(ax2, 30);
    ylabel(ax2, 'Mean noise correlation');
    title(ax2, sprintf('RF distance > %d px | %d-%d ms', OUT.P.summaryMinRfPx, tw(1), tw(2)));
    grid(ax2, 'on');
end
end
