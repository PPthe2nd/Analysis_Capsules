function R = analyze_anchor_knn_timeseries(D, Tall_V1, OUT3, stimID_example, varargin)
% ANALYZE_ANCHOR_KNN_TIMESERIES
% Select valid anchor combinations and compare raw vs global pixel-space KNN traces.
%
% Inputs:
%   D               : struct from post_affine_delta_points_allbins*.mat (field D.bins)
%   Tall_V1         : geometry table struct array
%   OUT3            : attention result struct (for significance mask)
%   stimID_example  : stimulus index used for output naming/title context
%
% Name/value options:
%   'siteRange'   (default 1:512)
%   'siteMask'    (default []) logical mask over absolute site index; if
%                 provided, this overrides pThresh-based selection
%   'pThresh'     (default 0.05)
%   'K'           (default 30)
%   'nPick'       (default 5)
%   'makePlot'    (default true)
%   'verbose'     (default true)
%   'saveFile'    (default '')  % optional .mat output
%
% Output R fields:
%   R.AnchorTable, R.tcRaw, R.tcKNN, R.tCenter, R.minAlong, R.maxAlong,
%   R.KNNContrib (cell of tables), R.TAlong

p = inputParser;
p.addParameter('siteRange', 1:512, @(x) isnumeric(x) && isvector(x));
p.addParameter('siteMask', [], @(x) isempty(x) || (islogical(x) && isvector(x)));
p.addParameter('pThresh', 0.05, @(x) isnumeric(x) && isscalar(x));
p.addParameter('K', 30, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('nPick', 5, @(x) isnumeric(x) && isscalar(x) && x >= 2);
p.addParameter('makePlot', true, @(x) islogical(x) && isscalar(x));
p.addParameter('verbose', true, @(x) islogical(x) && isscalar(x));
p.addParameter('saveFile', '', @(x) ischar(x) || isstring(x));
p.parse(varargin{:});
opt = p.Results;

assert(isstruct(D) && isfield(D, 'bins'), 'D must contain field bins.');

B = D.bins;
nBins = numel(B);
siteRange = opt.siteRange(:);

if isempty(siteRange)
    error('siteRange is empty.');
end
if ~isfield(OUT3, 'pValueTD') || numel(OUT3.pValueTD) < max(siteRange)
    error('OUT3.pValueTD is missing or smaller than siteRange.');
end

if isempty(opt.siteMask)
    sigMaskGlobal = false(max(siteRange),1);
    sigMaskGlobal(siteRange) = isfinite(OUT3.pValueTD(siteRange)) & ...
        (OUT3.pValueTD(siteRange) < opt.pThresh);
else
    sigMaskIn = opt.siteMask(:);
    if numel(sigMaskIn) < max(siteRange)
        error('siteMask has only %d elements; need at least %d.', numel(sigMaskIn), max(siteRange));
    end
    sigMaskGlobal = false(numel(sigMaskIn),1);
    sigMaskGlobal(1:numel(sigMaskIn)) = sigMaskIn;
end

% 1) Build set of fully valid (site, quartet, stimSource) in each bin, then intersect.
pairSets = cell(nBins,1);
for tb = 1:nBins
    xTb = double(B(tb).x_px);
    yTb = double(B(tb).y_px);
    sTb = double(B(tb).siteIdx);
    qTb = double(B(tb).quartetIdx);
    stTb = double(B(tb).stimIdxSource);
    dTb = double(B(tb).delta);

    ok = isfinite(xTb) & isfinite(yTb) & isfinite(sTb) & isfinite(qTb) & ...
         isfinite(stTb) & isfinite(dTb);
    pairSets{tb} = unique([sTb(ok), qTb(ok), stTb(ok)], 'rows');
end

Pall = pairSets{1};
for tb = 2:nBins
    Pall = intersect(Pall, pairSets{tb}, 'rows');
end
if isempty(Pall)
    error('No anchor combinations are present in all bins (finite x,y,delta).');
end

% Keep only requested sites.
keepSig = false(size(Pall,1),1);
for i = 1:size(Pall,1)
    s = Pall(i,1);
    if s >= 1 && s <= numel(sigMaskGlobal)
        keepSig(i) = sigMaskGlobal(s);
    end
end
Pall = Pall(keepSig,:);
if isempty(Pall)
    error('No all-bin anchor combinations remain after site filtering.');
end

% 2) Attach along_GC from source stimulus and select spaced anchors.
along = nan(size(Pall,1),1);
for i = 1:size(Pall,1)
    s0 = Pall(i,1);
    st0 = Pall(i,3);

    if st0 < 1 || st0 > numel(Tall_V1) || ~isfield(Tall_V1(st0), 'T')
        continue;
    end
    T = Tall_V1(st0).T;
    if s0 < 1 || s0 > height(T)
        continue;
    end
    along(i) = double(T.along_GC(s0));
end

okAlong = isfinite(along);
Pall = Pall(okAlong,:);
along = along(okAlong);
if numel(along) < opt.nPick
    error('Only %d valid anchors with finite along_GC; need %d.', ...
        numel(along), opt.nPick);
end

qTargets = quantile(along, linspace(0.1, 0.9, opt.nPick));
picked = false(size(along));
idxSel = nan(opt.nPick,1);
for k = 1:opt.nPick
    [~, ord] = sort(abs(along - qTargets(k)), 'ascend');
    j = ord(find(~picked(ord), 1, 'first'));
    if isempty(j)
        break;
    end
    idxSel(k) = j;
    picked(j) = true;
end
idxSel = idxSel(isfinite(idxSel));
if numel(idxSel) < opt.nPick
    rem = find(~picked);
    need = opt.nPick - numel(idxSel);
    idxSel = [idxSel; rem(1:need)]; %#ok<AGROW>
end

A = Pall(idxSel,:);
aSel = along(idxSel);
AnchorTable = table(A(:,1), A(:,3), A(:,2), aSel, ...
    'VariableNames', {'siteIdx','stimIdx','quartetIdx','along_GC'});
AnchorTable = sortrows(AnchorTable, 'along_GC');

if opt.verbose
    fprintf('Selected valid anchors (N=%d):\n', height(AnchorTable));
    disp(AnchorTable);
end

% 3) Extract raw + global KNN traces and along_GC span of KNN groups.
nC = height(AnchorTable);
tcRaw = nan(nC, nBins);
tcKNN = nan(nC, nBins);
tCenter = nan(nBins,1);
minAlong = nan(nC, nBins);
maxAlong = nan(nC, nBins);
knnSiteStim = cell(nC,1);
for i = 1:nC
    knnSiteStim{i} = zeros(0,2);
end

for tb = 1:nBins
    tCenter(tb) = mean(double(B(tb).timeWindow));

    xTb = double(B(tb).x_px);
    yTb = double(B(tb).y_px);
    dTb = double(B(tb).delta);
    sTb = double(B(tb).siteIdx);
    qTb = double(B(tb).quartetIdx);
    stTb = double(B(tb).stimIdxSource);

    finiteAll = isfinite(xTb) & isfinite(yTb) & isfinite(dTb) & ...
                isfinite(sTb) & isfinite(qTb) & isfinite(stTb);
    if ~any(finiteAll)
        continue;
    end

    xAll = xTb(finiteAll);
    yAll = yTb(finiteAll);
    dAll = dTb(finiteAll);
    sAll = sTb(finiteAll);
    stAll = stTb(finiteAll);

    alongAll = nan(size(sAll));
    for j = 1:numel(sAll)
        Tj = Tall_V1(stAll(j)).T;
        alongAll(j) = double(Tj.along_GC(sAll(j)));
    end

    for i = 1:nC
        s0 = double(AnchorTable.siteIdx(i));
        q0 = double(AnchorTable.quartetIdx(i));
        st0 = double(AnchorTable.stimIdx(i));

        mA = finiteAll & (sTb == s0) & (qTb == q0) & (stTb == st0);
        if ~any(mA)
            continue;
        end

        x0 = mean(xTb(mA));
        y0 = mean(yTb(mA));
        tcRaw(i,tb) = mean(dTb(mA), 'omitnan');

        dist2 = (xAll - x0).^2 + (yAll - y0).^2;
        [~, ord] = sort(dist2, 'ascend');
        if numel(ord) < opt.K
            error(['Bin %d has only %d projected points after filtering; ' ...
                   'cannot enforce K=%d global neighbors.'], tb, numel(ord), opt.K);
        end
        idxK = ord(1:opt.K);

        tcKNN(i,tb) = mean(dAll(idxK), 'omitnan');
        knnSiteStim{i} = [knnSiteStim{i}; [sAll(idxK), stAll(idxK)]]; %#ok<AGROW>

        aK = alongAll(idxK);
        aK = aK(isfinite(aK));
        if ~isempty(aK)
            minAlong(i,tb) = min(aK);
            maxAlong(i,tb) = max(aK);
        end
    end
end

KNNContrib = cell(nC,1);
for i = 1:nC
    SS = knnSiteStim{i};
    if isempty(SS)
        KNNContrib{i} = table([], [], [], ...
            'VariableNames', {'siteIdx','stimIdxSource','countUsed'});
    else
        [uSS, ~, ic] = unique(SS, 'rows');
        counts = accumarray(ic, 1);
        Tcon = table(uSS(:,1), uSS(:,2), counts, ...
            'VariableNames', {'siteIdx','stimIdxSource','countUsed'});
        KNNContrib{i} = sortrows(Tcon, {'countUsed','siteIdx'}, {'descend','ascend'});
    end
end

TAlong = table(AnchorTable.siteIdx, AnchorTable.stimIdx, AnchorTable.quartetIdx, ...
    AnchorTable.along_GC, min(minAlong,[],2,'omitnan'), max(maxAlong,[],2,'omitnan'), ...
    'VariableNames', {'siteIdx','stimIdx','quartetIdx','anchorAlong_GC','minAlong_KNN','maxAlong_KNN'});

if opt.verbose
    fprintf('along_GC span per anchor (across bins, K=%d):\n', opt.K);
    disp(TAlong);
    fprintf('Finite points per raw trace:\n');
    disp(sum(isfinite(tcRaw),2));
    fprintf('Finite points per KNN trace:\n');
    disp(sum(isfinite(tcKNN),2));

    fprintf('\n==== KNN contributor summary (siteIdx, stimIdxSource) ====\n');
    for i = 1:nC
        fprintf('\nAnchor %d: site %d | stim %d | q %d | along %.4f\n', ...
            i, AnchorTable.siteIdx(i), AnchorTable.stimIdx(i), ...
            AnchorTable.quartetIdx(i), AnchorTable.along_GC(i));
        if isempty(KNNContrib{i}) || height(KNNContrib{i}) == 0
            fprintf('    no contributors found\n');
        else
            disp(KNNContrib{i});
            fprintf('    along_GC span of KNN groups: [%.4f, %.4f]\n', ...
                TAlong.minAlong_KNN(i), TAlong.maxAlong_KNN(i));
        end
    end
end

if opt.makePlot
    figure('Color','w');
    tiledlayout(2,1,'Padding','compact','TileSpacing','compact');
    cols = lines(nC);

    ax1 = nexttile; hold on;
    for i = 1:nC
        plot(tCenter, tcRaw(i,:), 'LineWidth', 1.5, 'Color', cols(i,:));
    end
    xline(0,'--k'); yline(0,'--k');
    title('Raw anchor traces (valid anchors only)');
    ylabel('\Delta (T-D)');
    grid on;

    ax2 = nexttile; hold on;
    for i = 1:nC
        plot(tCenter, tcKNN(i,:), 'LineWidth', 1.8, 'Color', cols(i,:));
    end
    xline(0,'--k'); yline(0,'--k');
    title(sprintf('Global pixel-space KNN averaged traces (K=%d)', opt.K));
    xlabel('Time (ms)');
    ylabel('\Delta (T-D)');
    grid on;
    linkaxes([ax1 ax2], 'x');

    legend(ax2, compose('site %d | stim %d | q %d | along %.3f', ...
        AnchorTable.siteIdx, AnchorTable.stimIdx, AnchorTable.quartetIdx, ...
        AnchorTable.along_GC), 'Location','eastoutside');
end

R = struct();
R.meta = struct('stimID_example', stimID_example, 'K', opt.K, 'nPick', opt.nPick, ...
    'pThresh', opt.pThresh, 'siteRange', siteRange(:)');
R.AnchorTable = AnchorTable;
R.Comb5 = AnchorTable;
R.tcRaw = tcRaw;
R.tcKNN = tcKNN;
R.tCenter = tCenter;
R.minAlong = minAlong;
R.maxAlong = maxAlong;
R.TAlong = TAlong;
R.KNNContrib = KNNContrib;

saveFile = char(opt.saveFile);
if ~isempty(saveFile)
    save(saveFile, 'R', '-v7.3');
    if opt.verbose
        fprintf('Saved anchor/KNN diagnostics to: %s\n', saveFile);
    end
end

end
