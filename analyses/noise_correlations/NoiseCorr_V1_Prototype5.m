function OUT = NoiseCorr_V1_Prototype5(P)
% NOISECORR_V1_PROTOTYPE5
% Exploratory V1 noise-correlation prototype using the first 5 included V1
% sites. Residuals are computed within stimulus x session, scaled by the
% residual SD within site x session x window, and then pooled across
% sessions to obtain pairwise noise correlations.
%
% This is intentionally lightweight and is meant to validate the
% single-trial logic before building the full all-area cache.

if nargin < 1 || isempty(P)
    P = struct();
end

if ~isfield(P, 'Monkey'), P.Monkey = 1; end                 % 1 = Nilson, 2 = Figaro
if ~isfield(P, 'nSites'), P.nSites = 5; end
if ~isfield(P, 'timeWindows'), P.timeWindows = [-200 0; 40 240; 300 500]; end
if ~isfield(P, 'onlyCorrect'), P.onlyCorrect = true; end
if ~isfield(P, 'correctCol'), P.correctCol = 9; end
if ~isfield(P, 'correctVal'), P.correctVal = 1; end
if ~isfield(P, 'stimCol'), P.stimCol = 1; end
if ~isfield(P, 'dayCol'), P.dayCol = 11; end
if ~isfield(P, 'sessions'), P.sessions = [1 2]; end
if ~isfield(P, 'minTrialsPerStimSess'), P.minTrialsPerStimSess = 2; end
if ~isfield(P, 'useCache'), P.useCache = true; end
if ~isfield(P, 'saveCache'), P.saveCache = true; end
if ~isfield(P, 'plotFigures'), P.plotFigures = true; end
if ~isfield(P, 'nExamplePairs'), P.nExamplePairs = 3; end
if ~isfield(P, 'verbose'), P.verbose = true; end

cfg = config();
[monkeySuffix, monkeyFolder] = local_monkey_info(P.Monkey);
hasSessionExclusions = ~isempty(site_session_exclusions(monkeySuffix));
cacheFile = fullfile(cfg.resultsDir, sprintf('NoiseCorr_V1_Prototype5_%s.mat', char(monkeySuffix)));
cacheParams = local_cache_params(P);

if P.useCache && exist(cacheFile, 'file') == 2
    S = load(cacheFile, 'OUT');
    if isfield(S, 'OUT') && local_cache_matches(S.OUT, cacheParams) && ...
            session_exclusion_cache_matches(S.OUT, monkeySuffix)
        OUT = S.OUT;
        OUT.P = P;
        OUT.cacheFile = cacheFile;
        if P.verbose
            fprintf('Loaded cached V1 noise-correlation prototype from %s\n', cacheFile);
        end
        if P.plotFigures
            siteIds = double(OUT.SiteTable.globalSite);
            pairSiteIdx = [double(OUT.PairTable.site1Idx), double(OUT.PairTable.site2Idx)];
            local_plot_heatmaps(OUT.pairCorr, OUT.pairN, pairSiteIdx, siteIds, OUT.timeWindows);
            local_plot_histograms(OUT.pairCorr, OUT.timeWindows);
            local_plot_example_pairs(OUT.ResZ, OUT.pairCorr, OUT.pairN, pairSiteIdx, siteIds, OUT.timeWindows, P.nExamplePairs);
        end
        return;
    end
end

tallPath = fullfile(cfg.matDir, sprintf('Tall_V1_lines_%s.mat', char(monkeySuffix)));
assert(exist(tallPath, 'file') == 2, 'Missing %s.', tallPath);
Sgeo = load(tallPath);
assert(isfield(Sgeo, 'Tall_V1') && isstruct(Sgeo.Tall_V1), '%s must contain Tall_V1.', tallPath);
Tall_V1 = Sgeo.Tall_V1;

resp3binPath = fullfile(cfg.matDir, sprintf('SNR_capsules_%s_d12.mat', char(monkeySuffix)));
R_snr = load_capsules_struct_exclusion_aware(resp3binPath, monkeySuffix, 'cfg', cfg);

[keepSites, keepInfo] = local_select_v1_sites(cfg, monkeySuffix, Tall_V1, R_snr, hasSessionExclusions);
assert(numel(keepSites) >= P.nSites, 'Only %d V1 sites passed inclusion; need %d.', numel(keepSites), P.nSites);
siteIds = keepSites(1:P.nSites);
nSites = numel(siteIds);
pairSiteIdx = nchoosek(1:nSites, 2);
nPairs = size(pairSiteIdx, 1);
nStim = 384;
nSess = numel(P.sessions);
nWin = size(P.timeWindows, 1);

if P.verbose
    fprintf('NoiseCorr V1 prototype (%s): using sites %s\n', char(monkeySuffix), mat2str(siteIds(:)'));
    fprintf('Included V1 sites available: %d; prototype pairs: %d\n', numel(keepSites), nPairs);
end

trialRespPath = fullfile(cfg.dataRoot, monkeyFolder, 'ObjAtt_lines_normMUA.mat');
trialInfoPath = fullfile(cfg.dataRoot, monkeyFolder, 'ObjAtt_lines_MUA_trials.mat');
assert(exist(trialRespPath, 'file') == 2, 'Missing %s.', trialRespPath);
assert(exist(trialInfoPath, 'file') == 2, 'Missing %s.', trialInfoPath);

m1 = matfile(trialRespPath);
m2 = matfile(trialInfoPath);
ALLMAT = m2.ALLMAT;
tb = double(m2.tb);
tb = tb(:)';
assert(size(ALLMAT,2) >= max([P.stimCol, P.correctCol, P.dayCol]), ...
    'ALLMAT lacks required columns.');

stimPerTrial = double(ALLMAT(:, P.stimCol));
isCorrect = true(size(stimPerTrial));
if P.onlyCorrect
    isCorrect = (double(ALLMAT(:, P.correctCol)) == P.correctVal);
end
isDay = ismember(double(ALLMAT(:, P.dayCol)), P.sessions(:));
stimOK = isfinite(stimPerTrial) & stimPerTrial >= 1 & stimPerTrial <= nStim & ...
    (floor(stimPerTrial) == stimPerTrial);
isInclude = isCorrect & isDay & stimOK;

trialIdx = find(isInclude);
stimIncl = double(stimPerTrial(trialIdx));
dayIncl = double(ALLMAT(trialIdx, P.dayCol));
nTrialsIncl = numel(trialIdx);

if P.verbose
    fprintf('Using %d correct trials from sessions %s.\n', nTrialsIncl, mat2str(P.sessions));
end

siteBlock = siteIds(1):siteIds(end);
Xblock = double(m1.normMUA(siteBlock, :, :));  % contiguous channel read for MatFile compatibility
Xall = Xblock(siteIds - siteIds(1) + 1, :, :);
X = Xall(:, trialIdx, :);                      % [nSites x nTrialsIncl x nTime]
winIdx = local_window_indices(tb, P.timeWindows);

RespWin = nan(nSites, nTrialsIncl, nWin);
for w = 1:nWin
    ii = winIdx(w,1):winIdx(w,2);
    RespWin(:,:,w) = mean(X(:,:,ii), 3, 'omitnan');
end

siteKeepSess = true(nSites, nSess);
Texcl = site_session_exclusions(monkeySuffix);
for i = 1:height(Texcl)
    j = find(siteIds == double(Texcl.siteGlobal(i)), 1, 'first');
    k = find(P.sessions == double(Texcl.day(i)), 1, 'first');
    if ~isempty(j) && ~isempty(k)
        siteKeepSess(j, k) = false;
    end
end

nTrialsStimSess = zeros(nStim, nSess, 'uint16');
for k = 1:nSess
    mSess = dayIncl == P.sessions(k);
    nTrialsStimSess(:, k) = uint16(accumarray(stimIncl(mSess), 1, [nStim 1], @sum, 0));
end

Res = nan(nSites, nTrialsIncl, nWin);
ResZ = nan(nSites, nTrialsIncl, nWin);
nTrialsSiteStimSess = zeros(nSites, nStim, nSess, 'uint16');
meanStimSess = nan(nSites, nStim, nSess, nWin);
sdResidSess = nan(nSites, nSess, nWin);

for s = 1:nSites
    for k = 1:nSess
        if ~siteKeepSess(s, k)
            continue;
        end
        mSess = (dayIncl == P.sessions(k));
        for stim = 1:nStim
            idxStim = find(mSess & (stimIncl == stim));
            if numel(idxStim) < P.minTrialsPerStimSess
                continue;
            end
            for w = 1:nWin
                x = reshape(RespWin(s, idxStim, w), [], 1);
                good = isfinite(x);
                if nnz(good) < P.minTrialsPerStimSess
                    continue;
                end
                mu = mean(x(good));
                meanStimSess(s, stim, k, w) = mu;
                Res(s, idxStim(good), w) = x(good) - mu;
                nTrialsSiteStimSess(s, stim, k) = uint16(nnz(good));
            end
        end
        for w = 1:nWin
            r = reshape(Res(s, mSess, w), [], 1);
            good = isfinite(r);
            if nnz(good) < 2
                continue;
            end
            sig = std(r(good), 0);
            if ~isfinite(sig) || sig <= 0
                continue;
            end
            sdResidSess(s, k, w) = sig;
            idxUse = find(mSess);
            z = reshape(Res(s, idxUse, w), [], 1);
            goodZ = isfinite(z);
            z(goodZ) = z(goodZ) ./ sig;
            ResZ(s, idxUse, w) = z;
        end
    end
end

sumX2 = zeros(nSites, nStim, nSess, nWin, 'single');
sumXY = zeros(nPairs, nStim, nSess, nWin, 'single');
nOverlapPairStimSess = zeros(nPairs, nStim, nSess, nWin, 'uint16');

for k = 1:nSess
    mSess = (dayIncl == P.sessions(k));
    for stim = 1:nStim
        idxStim = find(mSess & (stimIncl == stim));
        if isempty(idxStim)
            continue;
        end
        for w = 1:nWin
            for s = 1:nSites
                x = reshape(ResZ(s, idxStim, w), [], 1);
                good = isfinite(x);
                if any(good)
                    sumX2(s, stim, k, w) = single(sum(x(good).^2));
                end
            end
            for p = 1:nPairs
                s1 = pairSiteIdx(p,1);
                s2 = pairSiteIdx(p,2);
                x = reshape(ResZ(s1, idxStim, w), [], 1);
                y = reshape(ResZ(s2, idxStim, w), [], 1);
                good = isfinite(x) & isfinite(y);
                if any(good)
                    sumXY(p, stim, k, w) = single(sum(x(good) .* y(good)));
                    nOverlapPairStimSess(p, stim, k, w) = uint16(nnz(good));
                end
            end
        end
    end
end

pairCorr = nan(nPairs, nWin);
pairN = zeros(nPairs, nWin);
for p = 1:nPairs
    s1 = pairSiteIdx(p,1);
    s2 = pairSiteIdx(p,2);
    for w = 1:nWin
        x = reshape(ResZ(s1, :, w), [], 1);
        y = reshape(ResZ(s2, :, w), [], 1);
        good = isfinite(x) & isfinite(y);
        pairN(p, w) = nnz(good);
        if pairN(p, w) >= 2
            Rtmp = corrcoef(x(good), y(good));
            pairCorr(p, w) = Rtmp(1,2);
        end
    end
end

SiteTable = table((1:nSites).', repmat("V1", nSites, 1), siteIds(:), siteIds(:), ...
    'VariableNames', {'siteIdx','area','globalSite','areaLocalSite'});
PairTable = table((1:nPairs).', pairSiteIdx(:,1), pairSiteIdx(:,2), repmat("V1-V1", nPairs, 1), ...
    'VariableNames', {'pairIdx','site1Idx','site2Idx','pairClass'});

OUT = struct();
OUT.P = P;
OUT.monkeySuffix = monkeySuffix;
OUT.timeWindows = P.timeWindows;
OUT.sessions = P.sessions(:);
OUT.SiteTable = SiteTable;
OUT.PairTable = PairTable;
OUT.keepInfo = keepInfo;
OUT.nTrials = nTrialsStimSess;
OUT.nTrialsSiteStimSess = nTrialsSiteStimSess;
OUT.siteKeepSess = siteKeepSess;
OUT.sumX2 = sumX2;
OUT.sumXY = sumXY;
OUT.nOverlapPairStimSess = nOverlapPairStimSess;
OUT.pairCorr = pairCorr;
OUT.pairN = pairN;
OUT.meanStimSess = meanStimSess;
OUT.sdResidSess = sdResidSess;
OUT.ResZ = ResZ;
OUT.cacheParams = cacheParams;
OUT.siteSessionExclusions = site_session_exclusions(monkeySuffix);
OUT.cacheFile = cacheFile;

if P.saveCache
    save(cacheFile, 'OUT', '-v7.3');
    if P.verbose
        fprintf('Saved cached V1 noise-correlation prototype to %s\n', cacheFile);
    end
end

if P.plotFigures
    local_plot_heatmaps(pairCorr, pairN, pairSiteIdx, siteIds, P.timeWindows);
    local_plot_histograms(pairCorr, P.timeWindows);
    local_plot_example_pairs(ResZ, pairCorr, pairN, pairSiteIdx, siteIds, P.timeWindows, P.nExamplePairs);
end

end

function S = local_cache_params(P)
S = struct();
S.cacheVersion = 1;
S.Monkey = double(P.Monkey);
S.nSites = double(P.nSites);
S.timeWindows = double(P.timeWindows);
S.onlyCorrect = logical(P.onlyCorrect);
S.correctCol = double(P.correctCol);
S.correctVal = double(P.correctVal);
S.stimCol = double(P.stimCol);
S.dayCol = double(P.dayCol);
S.sessions = double(P.sessions(:)');
S.minTrialsPerStimSess = double(P.minTrialsPerStimSess);
end

function tf = local_cache_matches(OUT, cacheParams)
tf = isstruct(OUT) && isfield(OUT, 'cacheParams') && isequaln(OUT.cacheParams, cacheParams);
end

function [keepSites, info] = local_select_v1_sites(cfg, monkeySuffix, Tall_V1, R_snr, hasSessionExclusions)
snrCachePath = fullfile(cfg.matDir, 'SNR_V1_byColor_byWindow.mat');
colorCachePath = fullfile(cfg.matDir, 'ColorTune_balanced_V1.mat');

if exist(snrCachePath, 'file') == 2 && ~hasSessionExclusions
    S = load(snrCachePath);
    assert(isfield(S, 'SNR'), '%s must contain SNR.', snrCachePath);
    SNR = S.SNR;
else
    SNR = compute_snr_per_color_sites(R_snr, Tall_V1, (1:512).', 'Verbose', false);
end

SNRmat = [SNR.yellowEarly(1:512), SNR.yellowLate(1:512), ...
    SNR.purpleEarly(1:512), SNR.purpleLate(1:512)];
[bestSNR, ~] = max(SNRmat, [], 2, 'omitnan');
bestSNR = bestSNR(:);

SNRthr = 0.7;
pTDthr = 0.05;
NminMatched = 20;
pColorThr = 0.05;

optsTD = struct('timeIdx', 3, 'excludeOverlap', true, 'verbose', false);
OUTtd = attention_modulation_V1_3bin(R_snr, Tall_V1, SNR, optsTD);

if exist(colorCachePath, 'file') == 2 && ~hasSessionExclusions
    S = load(colorCachePath);
    assert(isfield(S, 'ColorTune'), '%s must contain ColorTune.', colorCachePath);
    ColorTune = S.ColorTune;
else
    mainSites = find(bestSNR > SNRthr);
    ColorTune = compute_color_tuning_balanced_sites(R_snr, Tall_V1, (1:512).', mainSites, 'Verbose', false);
    ColorTune.bestSNR = bestSNR;
end

isColorSig = isfinite(ColorTune.early.p(1:512)) & (ColorTune.early.p(1:512) < pColorThr);
isColorSig = isColorSig(:);
matchedN = OUTtd.wY(:) + OUTtd.wP(:);
isMain = isfinite(bestSNR) & (bestSNR > SNRthr);
isRescue = isfinite(OUTtd.pValueTD(:)) & (OUTtd.pValueTD(:) < pTDthr) & (matchedN >= NminMatched);
isKeep = isMain | isRescue | isColorSig;

keepSites = find(isKeep);

info = struct();
info.bestSNR = bestSNR;
info.isMain = isMain;
info.isRescue = isRescue;
info.isColorSig = isColorSig;
info.keepSites = keepSites;
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
        error('P.Monkey must be 1 (Nilson) or 2 (Figaro).');
end
end

function winIdx = local_window_indices(tb, timeWindows)
nWin = size(timeWindows, 1);
winIdx = zeros(nWin, 2);
for w = 1:nWin
    i1 = find(tb >= timeWindows(w,1), 1, 'first');
    i2 = find(tb <= timeWindows(w,2), 1, 'last');
    assert(~isempty(i1) && ~isempty(i2) && i2 >= i1, ...
        'Window [%g %g] ms does not overlap tb.', timeWindows(w,1), timeWindows(w,2));
    winIdx(w,:) = [i1 i2];
end
end

function local_plot_heatmaps(pairCorr, pairN, pairSiteIdx, siteIds, timeWindows)
nSites = numel(siteIds);
nWin = size(timeWindows, 1);
M = nan(nSites, nSites, nWin);
N = nan(nSites, nSites, nWin);
for w = 1:nWin
    M(:,:,w) = nan(nSites);
    N(:,:,w) = nan(nSites);
    for p = 1:size(pairSiteIdx, 1)
        i = pairSiteIdx(p,1);
        j = pairSiteIdx(p,2);
        M(i,j,w) = pairCorr(p,w);
        M(j,i,w) = pairCorr(p,w);
        N(i,j,w) = pairN(p,w);
        N(j,i,w) = pairN(p,w);
    end
end

finiteVals = pairCorr(isfinite(pairCorr));
if isempty(finiteVals)
    clim = [-0.1 0.1];
else
    lo = min(finiteVals);
    hi = max(finiteVals);
    if lo >= 0 || hi <= 0
        pad = max(0.02, 0.05 * max(hi - lo, eps));
        clim = [lo - pad, hi + pad];
    else
        cmax = max(abs([lo hi]));
        clim = [-1.05 * cmax, 1.05 * cmax];
    end
end

figure('Color', 'w', 'Name', 'V1 noise correlations: pairwise heatmaps');
tiledlayout(1, nWin, 'Padding', 'compact', 'TileSpacing', 'compact');
for w = 1:nWin
    ax = nexttile;
    him = imagesc(ax, M(:,:,w), clim);
    axis(ax, 'image');
    set(him, 'AlphaData', isfinite(M(:,:,w)));
    set(ax, 'Color', [0.92 0.92 0.92]);
    colormap(ax, parula(256));
    xticks(ax, 1:nSites);
    yticks(ax, 1:nSites);
    xticklabels(ax, string(siteIds(:)));
    yticklabels(ax, string(siteIds(:)));
    xlabel(ax, 'V1 global site');
    ylabel(ax, 'V1 global site');
    title(ax, sprintf('%d to %d ms', timeWindows(w,1), timeWindows(w,2)));
    end
for w = 1:nWin
    ax = nexttile(w);
    colorbar(ax);
end
end

function local_plot_histograms(pairCorr, timeWindows)
colors = [0.35 0.35 0.35; 0.15 0.45 0.85; 0.85 0.40 0.15];
figure('Color', 'w', 'Name', 'V1 noise correlations: histograms');
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
for w = 1:size(pairCorr, 2)
    ax = nexttile;
    x = pairCorr(:, w);
    histogram(ax, x(isfinite(x)), 'BinWidth', 0.05, 'FaceColor', colors(w,:), 'FaceAlpha', 0.8);
    xlabel(ax, 'Noise correlation');
    ylabel(ax, 'N pairs');
    title(ax, sprintf('%d to %d ms', timeWindows(w,1), timeWindows(w,2)));
    grid(ax, 'on');
end
end

function local_plot_example_pairs(ResZ, pairCorr, pairN, pairSiteIdx, siteIds, timeWindows, nExamplePairs)
nPairs = size(pairSiteIdx, 1);
nExamplePairs = min(nExamplePairs, nPairs);
[~, order] = sort(abs(pairCorr(:, 3)), 'descend', 'MissingPlacement', 'last');
pairSel = order(1:nExamplePairs);

figure('Color', 'w', 'Name', 'V1 noise correlations: example pair scatters');
tiledlayout(nExamplePairs, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
for i = 1:nExamplePairs
    p = pairSel(i);
    s1 = pairSiteIdx(p,1);
    s2 = pairSiteIdx(p,2);
    for w = 1:3
        ax = nexttile;
        x = reshape(ResZ(s1, :, w), [], 1);
        y = reshape(ResZ(s2, :, w), [], 1);
        good = isfinite(x) & isfinite(y);
        scatter(ax, x(good), y(good), 12, 'filled', ...
            'MarkerFaceColor', [0.25 0.25 0.25], 'MarkerFaceAlpha', 0.35);
        hold(ax, 'on');
        if nnz(good) >= 2
            pp = polyfit(x(good), y(good), 1);
            xx = linspace(min(x(good)), max(x(good)), 100);
            plot(ax, xx, polyval(pp, xx), 'r-', 'LineWidth', 1.0);
        end
        xlabel(ax, sprintf('site %d z-resid', siteIds(s1)));
        ylabel(ax, sprintf('site %d z-resid', siteIds(s2)));
        title(ax, sprintf('pair %d | %d-%d ms | r=%.3f | N=%d', ...
            p, timeWindows(w,1), timeWindows(w,2), pairCorr(p,w), pairN(p,w)));
        grid(ax, 'on');
    end
end
end
