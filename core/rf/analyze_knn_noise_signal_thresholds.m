function R = analyze_knn_noise_signal_thresholds(D, Tall_V1, varargin)
% ANALYZE_KNN_NOISE_SIGNAL_THRESHOLDS
% Quantify pre-stim noise and late signal after post-affine pixel-space KNN pooling.
%
% Input:
%   D : struct with field D.bins (from compute_projected_delta_points_allbins)
%   Tall_V1 : struct array with Tall_V1(stim).T.along_GC
%
% Name/value options:
%   'KList'           (default [1 10 20 30 50 80])
%   'preEndMs'        (default 0)      % pre bins: window end <= preEndMs
%   'postStartMs'     (default 300)    % post bins: window start >= postStartMs
%   'preQuantilePct'  (default 95)     % threshold from pre-stim |value| percentile
%   'alphaFloorPct'   (default 80)     % suggested alpha floor percentile on pre-stim |value|
%   'cMaxPostPct'     (default 95)     % suggested cMax percentile on post-stim |value|
%   'kRef'            (default 30)     % K used for default movie suggestions
%   'enforceK'        (default true)   % error if a bin has fewer than max(KList) points
%   'siteWeightsByIndex' (default [])  % fixed per-site weights for KNN averaging
%   'makePlot'        (default true)
%   'verbose'         (default true)
%   'saveFile'        (default '')
%
% Output R:
%   R.summary         table with threshold and pre/post exceed metrics per K
%   R.tCenter         [nBins x 1] bin centers
%   R.preMask         [nBins x 1]
%   R.postMask        [nBins x 1]
%   R.suggested       struct with suggested alpha/color settings at kRef

p = inputParser;
p.addParameter('KList', [1 10 20 30 50 80], @(x) isnumeric(x) && isvector(x) && all(x >= 1));
p.addParameter('preEndMs', 0, @(x) isnumeric(x) && isscalar(x));
p.addParameter('postStartMs', 300, @(x) isnumeric(x) && isscalar(x));
p.addParameter('preQuantilePct', 95, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 100);
p.addParameter('alphaFloorPct', 80, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 100);
p.addParameter('cMaxPostPct', 95, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 100);
p.addParameter('kRef', 30, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('enforceK', true, @(x) islogical(x) && isscalar(x));
p.addParameter('siteWeightsByIndex', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.addParameter('makePlot', true, @(x) islogical(x) && isscalar(x));
p.addParameter('verbose', true, @(x) islogical(x) && isscalar(x));
p.addParameter('saveFile', '', @(x) ischar(x) || isstring(x));
p.parse(varargin{:});
opt = p.Results;

assert(isstruct(D) && isfield(D,'bins') && ~isempty(D.bins), 'D must contain non-empty field bins.');
assert(isstruct(Tall_V1) && ~isempty(Tall_V1), 'Tall_V1 must be a non-empty struct array.');
B = D.bins;
nBins = numel(B);

KList = unique(round(opt.KList(:)'));
if any(~isfinite(KList)) || isempty(KList)
    error('KList must contain finite positive integers.');
end
KMax = max(KList);
nK = numel(KList);

% Time masks
win = nan(nBins,2);
tCenter = nan(nBins,1);
for tb = 1:nBins
    win(tb,:) = double(B(tb).timeWindow);
    tCenter(tb) = mean(win(tb,:));
end
preMask = win(:,2) <= opt.preEndMs;
postMask = win(:,1) >= opt.postStartMs;
if ~any(preMask)
    error('No pre-stim bins found with window end <= %.1f ms.', opt.preEndMs);
end
if ~any(postMask)
    error('No post bins found with window start >= %.1f ms.', opt.postStartMs);
end

% Build global point index (site, quartet, stream) so target and distractor
% routes are tracked separately, matching the movie renderer.
pairAll = zeros(0,3);
for tb = 1:nBins
    s = double(B(tb).siteIdx(:));
    q = double(B(tb).quartetIdx(:));
    d = double(B(tb).delta(:));
    x = double(B(tb).x_px(:));
    y = double(B(tb).y_px(:));
    if isfield(B(tb), 'stream')
        str = double(B(tb).stream(:));
    else
        str = ones(size(d));
    end
    ok = isfinite(s) & isfinite(q) & isfinite(d) & isfinite(x) & isfinite(y) & isfinite(str);
    if any(ok)
        pairAll = [pairAll; [s(ok), q(ok), str(ok)]]; %#ok<AGROW>
    end
end
pairAll = unique(pairAll, 'rows');
if isempty(pairAll)
    error('No finite projected points found in D.bins.');
end
nComb = size(pairAll,1);
if opt.verbose
    fprintf('KNN noise/signal prep: %d combos x %d bins\n', nComb, nBins);
    fprintf('Pre bins: %d | Post bins: %d\n', nnz(preMask), nnz(postMask));
end

% Accumulate values by K
V = nan(nComb, nBins, nK, 'single');
alongSpanSum = zeros(nK,1);
alongSpanN = zeros(nK,1);
alongSpanSumPre = zeros(nK,1);
alongSpanNPre = zeros(nK,1);
alongSpanSumPost = zeros(nK,1);
alongSpanNPost = zeros(nK,1);

for tb = 1:nBins
    s = double(B(tb).siteIdx(:));
    q = double(B(tb).quartetIdx(:));
    d = double(B(tb).delta(:));
    x = double(B(tb).x_px(:));
    y = double(B(tb).y_px(:));
    if isfield(B(tb), 'stream')
        str = double(B(tb).stream(:));
    else
        str = ones(size(d));
    end

    ok = isfinite(s) & isfinite(q) & isfinite(d) & isfinite(x) & isfinite(y) & isfinite(str);
    if ~any(ok)
        continue;
    end

    s = s(ok); q = q(ok); d = d(ok); x = x(ok); y = y(ok); str = str(ok);
    if isempty(opt.siteWeightsByIndex)
        wSite = ones(size(s));
    else
        wBySite = double(opt.siteWeightsByIndex(:));
        wSite = zeros(size(s));
        sRound = round(s);
        okW = isfinite(sRound) & sRound >= 1 & sRound <= numel(wBySite);
        wSite(okW) = wBySite(sRound(okW));
        wSite(~isfinite(wSite) | wSite < 0) = 0;
    end
    if ~isfield(B(tb), 'stimIdxSource')
        error('D.bins(%d) is missing stimIdxSource; re-run post-affine export.', tb);
    end
    st = double(B(tb).stimIdxSource(:));
    st = st(ok);
    along = nan(size(s));
    for ii = 1:numel(s)
        stimNum = st(ii);
        siteNum = s(ii);
        if stimNum >= 1 && stimNum <= numel(Tall_V1) && ...
                isfield(Tall_V1(stimNum), 'T') && istable(Tall_V1(stimNum).T)
            Ttab = Tall_V1(stimNum).T;
            if siteNum >= 1 && siteNum <= height(Ttab) && ismember('along_GC', Ttab.Properties.VariableNames)
                along(ii) = double(Ttab.along_GC(siteNum));
            end
        end
    end
    nP = numel(d);
    if nP < KMax
        if opt.enforceK
            error('Bin %d has only %d points; cannot enforce K=%d.', tb, nP, KMax);
        end
        kEff = KList(KList <= nP);
        if isempty(kEff)
            continue;
        end
    else
        kEff = KList;
    end

    % all K at once via cumulative neighbor sums
    vLocal = nan(nP, nK);
    spanLocal = nan(nP, nK);
    for i = 1:nP
        sameObj = (str == str(i));
        idxPool = find(sameObj);
        xPool = x(idxPool);
        yPool = y(idxPool);
        dPool = d(idxPool);
        wPool = wSite(idxPool);
        aPool = along(idxPool);

        dist2 = (xPool - x(i)).^2 + (yPool - y(i)).^2;
        [~, ord] = sort(dist2, 'ascend');
        kUseMax = min(KMax, numel(ord));
        if kUseMax < 1
            continue;
        end
        idx = ord(1:kUseMax);
        csum = cumsum(dPool(idx));
        csumW = cumsum(wPool(idx));
        csumWD = cumsum(wPool(idx) .* dPool(idx));
        for kk = 1:nK
            k = KList(kk);
            if k <= kUseMax
                if csumW(k) > 0
                    vLocal(i,kk) = csumWD(k) / csumW(k);
                else
                    vLocal(i,kk) = csum(k) / k;
                end
                a = aPool(idx(1:k));
                a = a(isfinite(a));
                if ~isempty(a)
                    spanLocal(i,kk) = max(a) - min(a);
                end
            end
        end
    end

    thisPairs = [s q str];
    [tf, loc] = ismember(thisPairs, pairAll, 'rows');
    if any(~tf)
        warning('Some pairs in bin %d were not found in global pairAll; skipping them.', tb);
    end
    loc = loc(tf);
    vals = vLocal(tf,:);
    spans = spanLocal(tf,:);
    strLoc = str(tf);

    % Apply target/distractor routing to match plotting semantics.
    for kk = 1:nK
        vals(:,kk) = route_by_stream(vals(:,kk), strLoc);
    end

    for kk = 1:nK
        V(loc, tb, kk) = single(vals(:,kk));
        sp = spans(:,kk);
        good = isfinite(sp);
        if any(good)
            alongSpanSum(kk) = alongSpanSum(kk) + sum(sp(good));
            alongSpanN(kk) = alongSpanN(kk) + nnz(good);
            if preMask(tb)
                alongSpanSumPre(kk) = alongSpanSumPre(kk) + sum(sp(good));
                alongSpanNPre(kk) = alongSpanNPre(kk) + nnz(good);
            elseif postMask(tb)
                alongSpanSumPost(kk) = alongSpanSumPost(kk) + sum(sp(good));
                alongSpanNPost(kk) = alongSpanNPost(kk) + nnz(good);
            end
        end
    end

    if opt.verbose && (tb == 1 || tb == nBins || mod(tb,10) == 0)
        fprintf('  processed bin %2d/%2d (%d points)\n', tb, nBins, nP);
    end
end

% Summary metrics per K
Kcol = zeros(nK,1);
thr = nan(nK,1);
preExceed = nan(nK,1);
postExceed = nan(nK,1);
preMean = nan(nK,1);
postMean = nan(nK,1);
preMedian = nan(nK,1);
postMedian = nan(nK,1);
dPrime = nan(nK,1);
alphaFloor = nan(nK,1);
cMaxPost = nan(nK,1);
meanAlongSpan = nan(nK,1);
meanAlongSpanPre = nan(nK,1);
meanAlongSpanPost = nan(nK,1);

for kk = 1:nK
    Kcol(kk) = KList(kk);
    A = double(V(:,:,kk));

    preVals = A(:, preMask);
    preVals = preVals(isfinite(preVals));
    postVals = A(:, postMask);
    postVals = postVals(isfinite(postVals));

    if isempty(preVals) || isempty(postVals)
        continue;
    end

    thr(kk) = prctile(preVals, opt.preQuantilePct);
    preExceed(kk) = mean(preVals > thr(kk));
    postExceed(kk) = mean(postVals > thr(kk));

    preMean(kk) = mean(preVals);
    postMean(kk) = mean(postVals);
    preMedian(kk) = median(preVals);
    postMedian(kk) = median(postVals);

    vpre = var(preVals, 1);
    vpost = var(postVals, 1);
    denom = sqrt(0.5*(vpre + vpost));
    if isfinite(denom) && denom > 0
        dPrime(kk) = (postMean(kk) - preMean(kk)) / denom;
    end

    alphaFloor(kk) = prctile(preVals, opt.alphaFloorPct);
    cMaxPost(kk) = prctile(postVals, opt.cMaxPostPct);
    if alongSpanN(kk) > 0
        meanAlongSpan(kk) = alongSpanSum(kk) / alongSpanN(kk);
    end
    if alongSpanNPre(kk) > 0
        meanAlongSpanPre(kk) = alongSpanSumPre(kk) / alongSpanNPre(kk);
    end
    if alongSpanNPost(kk) > 0
        meanAlongSpanPost(kk) = alongSpanSumPost(kk) / alongSpanNPost(kk);
    end
end

summary = table(Kcol, thr, preExceed, postExceed, preMean, postMean, preMedian, postMedian, dPrime, alphaFloor, cMaxPost, ...
    meanAlongSpan, meanAlongSpanPre, meanAlongSpanPost, ...
    'VariableNames', {'K','thresholdPreQ','preExceedFrac','postExceedFrac','preMeanAbs','postMeanAbs', ...
                      'preMedianAbs','postMedianAbs','dPrimeAbs','alphaFloorSuggest','cMaxSuggest', ...
                      'meanAlongSpanKNN','meanAlongSpanKNN_Pre','meanAlongSpanKNN_Post'});

% Suggested settings at reference K
if ismember(opt.kRef, KList)
    kRefUse = opt.kRef;
else
    [~,iNear] = min(abs(KList - opt.kRef));
    kRefUse = KList(iNear);
end
rowRef = find(summary.K == kRefUse, 1, 'first');

suggested = struct();
suggested.kRef = kRefUse;
suggested.alphaFullAt = summary.thresholdPreQ(rowRef);
suggested.alphaValueFloor = summary.alphaFloorSuggest(rowRef);
suggested.colorRedAt = summary.thresholdPreQ(rowRef);
suggested.cMaxFixed = summary.cMaxSuggest(rowRef);
suggested.preQuantilePct = opt.preQuantilePct;

if opt.verbose
    fprintf('\nK sweep summary (threshold from pre-stim |value| p%.1f):\n', opt.preQuantilePct);
    disp(summary);
    fprintf('Suggested movie settings at K=%d:\n', kRefUse);
    fprintf('  alphaFullAt      = %.6g\n', suggested.alphaFullAt);
    fprintf('  alphaValueFloor  = %.6g\n', suggested.alphaValueFloor);
    fprintf('  colorRedAt       = %.6g\n', suggested.colorRedAt);
    fprintf('  cMaxFixed        = %.6g\n', suggested.cMaxFixed);
end

if opt.makePlot
    figure('Color','w');
    tiledlayout(3,1,'Padding','compact','TileSpacing','compact');

    ax1 = nexttile; hold(ax1, 'on');
    plot(summary.K, summary.thresholdPreQ, '-o', 'LineWidth', 1.6, 'Color', [0.1 0.1 0.1], 'DisplayName', 'Threshold (pre quantile)');
    plot(summary.K, summary.alphaFloorSuggest, '-o', 'LineWidth', 1.4, 'Color', [0.45 0.45 0.45], 'DisplayName', 'Alpha floor suggest');
    plot(summary.K, summary.cMaxSuggest, '-o', 'LineWidth', 1.4, 'Color', [0.75 0.25 0.25], 'DisplayName', 'cMax suggest (post)');
    ylabel('|\Delta| scale');
    xlabel('K (pixel-space neighbors)');
    grid on;
    legend('Location','best');
    title(sprintf('Scale suggestions vs K | pre<=%.0f ms, post>=%.0f ms', opt.preEndMs, opt.postStartMs));

    ax2 = nexttile; hold(ax2, 'on');
    plot(summary.K, 100*summary.preExceedFrac, '-o', 'LineWidth', 1.6, 'Color', [0.30 0.30 0.30], 'DisplayName', 'Pre > threshold');
    plot(summary.K, 100*summary.postExceedFrac, '-o', 'LineWidth', 1.8, 'Color', [0.05 0.40 0.85], 'DisplayName', 'Post > threshold');
    ylabel('Exceedance (%)');
    xlabel('K (pixel-space neighbors)');
    yline(100 - opt.preQuantilePct, '--k', '5% target pre exceed');
    grid on;
    legend('Location','best');
    title('|\Delta| threshold hit-rate (noise vs signal)');

    ax3 = nexttile; hold(ax3, 'on');
    plot(summary.K, summary.meanAlongSpanKNN, '-o', 'LineWidth', 1.8, 'Color', [0.20 0.50 0.20], 'DisplayName', 'All bins');
    plot(summary.K, summary.meanAlongSpanKNN_Pre, '-o', 'LineWidth', 1.4, 'Color', [0.55 0.55 0.55], 'DisplayName', 'Pre bins');
    plot(summary.K, summary.meanAlongSpanKNN_Post, '-o', 'LineWidth', 1.4, 'Color', [0.05 0.40 0.85], 'DisplayName', 'Post bins');
    ylabel('Mean along\_GC range');
    xlabel('K (pixel-space neighbors)');
    grid on;
    legend('Location','best');
    title('Average along\_GC neighborhood span vs K');

    linkaxes([ax1 ax2 ax3], 'x');
end

R = struct();
R.meta = struct('created', datestr(now,30), 'nComb', nComb, 'nBins', nBins, ...
    'KList', KList, 'preEndMs', opt.preEndMs, 'postStartMs', opt.postStartMs, ...
    'preQuantilePct', opt.preQuantilePct, 'kRefRequested', opt.kRef, ...
    'useSiteWeights', ~isempty(opt.siteWeightsByIndex));
R.tCenter = tCenter;
R.preMask = preMask;
R.postMask = postMask;
R.summary = summary;
R.suggested = suggested;

saveFile = char(opt.saveFile);
if ~isempty(saveFile)
    save(saveFile, 'R', '-v7.3');
    if opt.verbose
        fprintf('Saved KNN noise/signal prep to: %s\n', saveFile);
    end
end

end

function v = route_by_stream(vSigned, stream)
isT = (stream == 1);
isD = (stream == 2);
v = zeros(size(vSigned));
v(isT) = max(0,  vSigned(isT));
v(isD) = max(0, -vSigned(isD));
isOther = ~(isT | isD);
v(isOther) = abs(vSigned(isOther));
end
