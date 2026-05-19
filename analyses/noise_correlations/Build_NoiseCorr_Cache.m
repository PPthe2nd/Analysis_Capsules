function OUT = Build_NoiseCorr_Cache(Monkey, Puser)
% BUILD_NOISECORR_CACHE
% Build a monkey-level noise-correlation cache across the standard kept
% V1, V4, and IT site sets. Trial residuals are computed within
% stimulus x session, normalized by the residual SD within site x session x
% window, and stored as per-site/per-pair moments for later regrouping.

if nargin < 1 || isempty(Monkey)
    Monkey = 1; % 1 = Nilson, 2 = Figaro
end
if nargin < 2 || isempty(Puser)
    Puser = struct();
end

P = struct();
P.timeWindows = [-200 0; 40 240; 300 500];
P.onlyCorrect = true;
P.correctCol = 9;
P.correctVal = 1;
P.stimCol = 1;
P.dayCol = 11;
P.sessions = [1 2];
P.minTrialsPerStimSess = 2;
P.chunkTrials = 200;
P.useCache = true;
P.saveCache = true;
P.plotFigure = false;
P.verbose = true;
P.siteLimitByArea = struct('V1', [], 'V4', [], 'IT', []);
P.cacheTag = "";

if ~isempty(Puser)
    fn = fieldnames(Puser);
    for i = 1:numel(fn)
        P.(fn{i}) = Puser.(fn{i});
    end
end

cfg = config();
[monkeySuffix, monkeyFolder] = local_monkey_info(Monkey);
cacheTag = local_cache_tag(P);
cacheFile = fullfile(cfg.resultsDir, sprintf('NoiseCorr_Cache_%s_%s.mat', char(monkeySuffix), cacheTag));
cacheParams = local_cache_params(Monkey, P, cacheTag);

if P.useCache && exist(cacheFile, 'file') == 2
    S = load(cacheFile, 'OUT');
    if isfield(S, 'OUT') && local_cache_matches(S.OUT, cacheParams) && ...
            session_exclusion_cache_matches(S.OUT, monkeySuffix) && ...
            local_pair_files_exist(S.OUT)
        OUT = S.OUT;
        OUT.P = P;
        OUT.cacheFile = cacheFile;
        if P.verbose
            fprintf('Loaded cached noise-correlation build from %s\n', cacheFile);
        end
        if P.plotFigure
            local_plot_summary(OUT);
        end
        return;
    end
end

resp3binPath = fullfile(cfg.matDir, sprintf('SNR_capsules_%s_d12.mat', char(monkeySuffix)));
assert(exist(resp3binPath, 'file') == 2, 'Missing %s.', resp3binPath);
R3full = load_capsules_struct_exclusion_aware(resp3binPath, monkeySuffix, 'cfg', cfg);

[SiteTable, AreaSummary] = local_build_site_table(cfg, monkeySuffix, R3full, P);
nSites = height(SiteTable);
assert(nSites >= 2, 'Need at least two included sites; found %d.', nSites);

[PairTable, PairClassTable, pairSiteIdx] = local_build_pair_table(SiteTable);
nPairs = height(PairTable);
nStim = 384;
nSess = numel(P.sessions);
nWin = size(P.timeWindows, 1);
tBuild = tic;

if P.verbose
    fprintf('NoiseCorr build (%s): %d sites total | %d pairs\n', char(monkeySuffix), nSites, nPairs);
    disp(AreaSummary);
    disp(PairClassTable);
end

trialRespPath = fullfile(cfg.dataRoot, monkeyFolder, 'ObjAtt_lines_normMUA.mat');
trialInfoPath = fullfile(cfg.dataRoot, monkeyFolder, 'ObjAtt_lines_MUA_trials.mat');
assert(exist(trialRespPath, 'file') == 2, 'Missing %s.', trialRespPath);
assert(exist(trialInfoPath, 'file') == 2, 'Missing %s.', trialInfoPath);

mResp = matfile(trialRespPath);
mInfo = matfile(trialInfoPath);
ALLMAT = mInfo.ALLMAT;
tb = double(mInfo.tb);
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
    fprintf('Using %d correct trials from sessions %s\n', nTrialsIncl, mat2str(P.sessions));
end

globalSite = double(SiteTable.globalSite(:));
assert(numel(unique(globalSite)) == nSites, 'Included globalSite values must be unique.');
siteBlock = min(globalSite):max(globalSite);
siteBlockSel = globalSite - siteBlock(1) + 1;
winIdx = local_window_indices(tb, P.timeWindows);

RespWin = nan(nSites, nTrialsIncl, nWin);
for t0 = 1:P.chunkTrials:nTrialsIncl
    t1 = min(nTrialsIncl, t0 + P.chunkTrials - 1);
    idxLocal = t0:t1;
    trialChunk = trialIdx(idxLocal);
    trialBlock = trialChunk(1):trialChunk(end);
    trialSel = trialChunk - trialBlock(1) + 1;
    Xblk = double(mResp.normMUA(siteBlock, trialBlock, :));
    Xsel = Xblk(siteBlockSel, trialSel, :);
    for w = 1:nWin
        ii = winIdx(w,1):winIdx(w,2);
        RespWin(:, idxLocal, w) = mean(Xsel(:,:,ii), 3, 'omitnan');
    end
    if P.verbose && (t0 == 1 || t1 == nTrialsIncl || mod(t1, 1000) == 0)
        local_print_progress('  Loaded trials', t1, nTrialsIncl, toc(tBuild));
    end
end

siteKeepSess = true(nSites, nSess);
Texcl = site_session_exclusions(monkeySuffix);
for i = 1:height(Texcl)
    j = find(globalSite == double(Texcl.siteGlobal(i)), 1, 'first');
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
            ResZ(s, mSess, w) = Res(s, mSess, w) ./ sig;
        end
    end
end

sumX2 = zeros(nSites, nStim, nSess, nWin, 'single');
sumXY = zeros(nPairs, nStim, nSess, nWin, 'single');
pairLin = sub2ind([nSites nSites], pairSiteIdx(:,1), pairSiteIdx(:,2));
totalMomentCells = nSess * nStim * nWin;
doneMomentCells = 0;

for k = 1:nSess
    mSess = (dayIncl == P.sessions(k));
    for stim = 1:nStim
        idxStim = find(mSess & (stimIncl == stim));
        if isempty(idxStim)
            continue;
        end
        for w = 1:nWin
            Xzw = reshape(ResZ(:, idxStim, w), nSites, []);
            sumX2(:, stim, k, w) = single(sum(Xzw.^2, 2, 'omitnan'));
            Xtmp = Xzw;
            Xtmp(~isfinite(Xtmp)) = 0;
            C = Xtmp * Xtmp.';
            sumXY(:, stim, k, w) = single(C(pairLin));
            doneMomentCells = doneMomentCells + 1;
            if P.verbose && (doneMomentCells == 1 || doneMomentCells == totalMomentCells || mod(doneMomentCells, 200) == 0)
                local_print_progress('  Accumulated moment cells', doneMomentCells, totalMomentCells, toc(tBuild));
            end
        end
    end
end

pairCorr = local_pair_corr_from_moments(sumX2, sumXY, pairSiteIdx);
siteSharedTrialMismatch = local_count_mismatch(nTrialsSiteStimSess, nTrialsStimSess, siteKeepSess, P.minTrialsPerStimSess);
PairClassTable = local_save_pair_class_files(cfg, monkeySuffix, cacheTag, PairClassTable, PairTable, sumXY, pairCorr, P, cacheParams, Texcl);

OUT = struct();
OUT.P = P;
OUT.monkeySuffix = monkeySuffix;
OUT.cacheFile = cacheFile;
OUT.cacheParams = cacheParams;
OUT.siteSessionExclusions = Texcl;
OUT.timeWindows = P.timeWindows;
OUT.sessions = P.sessions(:);
OUT.AreaSummary = AreaSummary;
OUT.SiteTable = SiteTable;
OUT.PairTable = PairTable;
OUT.PairClassTable = PairClassTable;
OUT.siteKeepSess = siteKeepSess;
OUT.nTrials = nTrialsStimSess;
OUT.nTrialsSiteStimSess = nTrialsSiteStimSess;
OUT.sdResidSess = sdResidSess;
OUT.sumX2 = sumX2;
OUT.pairCorr = pairCorr;
OUT.siteSharedTrialMismatch = siteSharedTrialMismatch;
OUT.buildElapsedSec = toc(tBuild);

if P.saveCache
    save(cacheFile, 'OUT', '-v7.3');
    if P.verbose
        fprintf('Saved full noise-correlation cache to %s\n', cacheFile);
    end
end

if P.plotFigure
    local_plot_summary(OUT);
end
end

function PairClassTable = local_save_pair_class_files(cfg, monkeySuffix, cacheTag, PairClassTable, PairTable, sumXY, pairCorr, P, cacheParams, Texcl)
pairFiles = strings(height(PairClassTable), 1);
for i = 1:height(PairClassTable)
    cls = string(PairClassTable.pairClass(i));
    s = PairClassTable.startIdx(i);
    e = PairClassTable.endIdx(i);
    pairFile = fullfile(cfg.resultsDir, sprintf('NoiseCorr_Pairs_%s_%s_%s.mat', char(monkeySuffix), cacheTag, local_pair_class_token(cls)));
    pairFiles(i) = string(pairFile);
    if ~isfinite(s) || ~isfinite(e) || PairClassTable.nPairs(i) == 0
        if P.saveCache
            PairOut = struct();
            PairOut.pairClass = cls;
            PairOut.cacheParams = cacheParams;
            PairOut.siteSessionExclusions = Texcl;
            PairOut.PairTable = PairTable([],:);
            PairOut.sumXY = zeros(0, size(sumXY,2), size(sumXY,3), size(sumXY,4), 'single');
            PairOut.pairCorr = zeros(0, size(pairCorr,2));
            PairOut.timeWindows = P.timeWindows;
            PairOut.sessions = P.sessions(:);
            save(pairFile, 'PairOut', '-v7.3');
        end
        continue;
    end
    PairOut = struct();
    PairOut.pairClass = cls;
    PairOut.cacheParams = cacheParams;
    PairOut.siteSessionExclusions = Texcl;
    PairOut.PairTable = PairTable(s:e, :);
    PairOut.sumXY = sumXY(s:e, :, :, :);
    PairOut.pairCorr = pairCorr(s:e, :);
    PairOut.timeWindows = P.timeWindows;
    PairOut.sessions = P.sessions(:);
    if P.saveCache
        save(pairFile, 'PairOut', '-v7.3');
        if P.verbose
            fprintf('Saved %s pair moments to %s\n', cls, pairFile);
        end
    end
end
PairClassTable.pairFile = pairFiles;
end

function [SiteTable, AreaSummary] = local_build_site_table(cfg, monkeySuffix, R3full, P)
hasSessionExclusions = ~isempty(site_session_exclusions(monkeySuffix));

% V1
tallV1Path = fullfile(cfg.matDir, sprintf('Tall_V1_lines_%s.mat', char(monkeySuffix)));
SV1 = load(tallV1Path, 'Tall_V1');
Tall_V1 = local_sort_tall(SV1.Tall_V1, 'Tall_V1');
[keepV1, infoV1] = local_select_v1_sites(cfg, monkeySuffix, Tall_V1, R3full, hasSessionExclusions); %#ok<ASGLU>
keepV1 = local_apply_limit(keepV1, P.siteLimitByArea, 'V1');

% V4
tallV4Path = fullfile(cfg.matDir, sprintf('Tall_V4_lines_%s.mat', char(monkeySuffix)));
SV4 = load(tallV4Path, 'Tall_V4', 'RFrange');
Tall_V4 = local_sort_tall(SV4.Tall_V4, 'Tall_V4');
RFrangeV4 = SV4.RFrange(:);
[keepV4Local, infoV4] = local_select_v4_sites(Tall_V4, R3full, RFrangeV4); %#ok<ASGLU>
keepV4Local = local_apply_limit(keepV4Local, P.siteLimitByArea, 'V4');
keepV4Global = RFrangeV4(keepV4Local);

% IT
tallITPath = fullfile(cfg.matDir, sprintf('Tall_IT_lines_%s.mat', char(monkeySuffix)));
SIT = load(tallITPath, 'Tall_IT', 'RFrange');
Tall_IT = local_sort_tall(SIT.Tall_IT, 'Tall_IT');
RFrangeIT = SIT.RFrange(:);
[keepITLocal, infoIT] = local_select_it_sites(Tall_IT, R3full, RFrangeIT); %#ok<ASGLU>
keepITLocal = local_apply_limit(keepITLocal, P.siteLimitByArea, 'IT');
keepITGlobal = RFrangeIT(keepITLocal);

siteIdx = (1:(numel(keepV1) + numel(keepV4Local) + numel(keepITLocal))).';
area = [repmat("V1", numel(keepV1), 1); ...
        repmat("V4", numel(keepV4Local), 1); ...
        repmat("IT", numel(keepITLocal), 1)];
globalSite = [double(keepV1(:)); double(keepV4Global(:)); double(keepITGlobal(:))];
areaLocalSite = [double(keepV1(:)); double(keepV4Local(:)); double(keepITLocal(:))];

SiteTable = table(siteIdx, area, globalSite, areaLocalSite, ...
    'VariableNames', {'siteIdx','area','globalSite','areaLocalSite'});

AreaSummary = table( ...
    ["V1"; "V4"; "IT"], ...
    [512; numel(RFrangeV4); numel(RFrangeIT)], ...
    [numel(keepV1); numel(keepV4Local); numel(keepITLocal)], ...
    'VariableNames', {'area','nAvailable','nIncluded'});
end

function [PairTable, PairClassTable, pairSiteIdx] = local_build_pair_table(SiteTable)
idxV1 = find(SiteTable.area == "V1");
idxV4 = find(SiteTable.area == "V4");
idxIT = find(SiteTable.area == "IT");

[pairVV, clsVV] = local_same_area_pairs(idxV1, "V1-V1");
[pairV4V4, clsV4V4] = local_same_area_pairs(idxV4, "V4-V4");
[pairITIT, clsITIT] = local_same_area_pairs(idxIT, "IT-IT");
[pairV1V4, clsV1V4] = local_cross_area_pairs(idxV1, idxV4, "V1-V4");
[pairV1IT, clsV1IT] = local_cross_area_pairs(idxV1, idxIT, "V1-IT");
[pairV4IT, clsV4IT] = local_cross_area_pairs(idxV4, idxIT, "V4-IT");

pairSiteIdx = [pairVV; pairV1V4; pairV1IT; pairV4V4; pairV4IT; pairITIT];
pairClass = [clsVV; clsV1V4; clsV1IT; clsV4V4; clsV4IT; clsITIT];
pairIdx = (1:size(pairSiteIdx,1)).';

PairTable = table(pairIdx, pairSiteIdx(:,1), pairSiteIdx(:,2), pairClass, ...
    'VariableNames', {'pairIdx','site1Idx','site2Idx','pairClass'});

classes = ["V1-V1"; "V1-V4"; "V1-IT"; "V4-V4"; "V4-IT"; "IT-IT"];
startIdx = nan(numel(classes),1);
endIdx = nan(numel(classes),1);
nPairs = zeros(numel(classes),1);
for i = 1:numel(classes)
    idx = find(PairTable.pairClass == classes(i));
    nPairs(i) = numel(idx);
    if ~isempty(idx)
        startIdx(i) = idx(1);
        endIdx(i) = idx(end);
    end
end
PairClassTable = table(classes, startIdx, endIdx, nPairs, ...
    'VariableNames', {'pairClass','startIdx','endIdx','nPairs'});
end

function [pairs, cls] = local_same_area_pairs(idx, label)
if numel(idx) < 2
    pairs = zeros(0,2);
else
    pairs = nchoosek(idx(:).', 2);
end
cls = repmat(string(label), size(pairs,1), 1);
end

function [pairs, cls] = local_cross_area_pairs(idxA, idxB, label)
if isempty(idxA) || isempty(idxB)
    pairs = zeros(0,2);
else
    [A, B] = ndgrid(idxA(:), idxB(:));
    pairs = [A(:), B(:)];
end
cls = repmat(string(label), size(pairs,1), 1);
end

function pairCorr = local_pair_corr_from_moments(sumX2, sumXY, pairSiteIdx)
nPairs = size(sumXY, 1);
nWin = size(sumXY, 4);
nSites = size(sumX2, 1);
pairCorr = nan(nPairs, nWin);

sitePower = zeros(nSites, nWin);
for w = 1:nWin
    sitePower(:, w) = sum(sum(double(sumX2(:,:,:,w)), 3), 2);
end

for w = 1:nWin
    pairCov = sum(sum(double(sumXY(:,:,:,w)), 3), 2);
    for p = 1:nPairs
        i = pairSiteIdx(p,1);
        j = pairSiteIdx(p,2);
        denom = sqrt(sitePower(i,w) * sitePower(j,w));
        if isfinite(denom) && denom > 0
            pairCorr(p,w) = pairCov(p) / denom;
        end
    end
end
end

function mismatch = local_count_mismatch(nTrialsSiteStimSess, nTrialsStimSess, siteKeepSess, minTrialsPerStimSess)
nSites = size(nTrialsSiteStimSess,1);
nSess = size(nTrialsStimSess,2);
mismatch = zeros(nSites, nSess);
for s = 1:nSites
    for k = 1:nSess
        if ~siteKeepSess(s,k)
            continue;
        end
        expCounts = double(nTrialsStimSess(:,k));
        gotCounts = double(squeeze(nTrialsSiteStimSess(s,:,k))).';
        use = expCounts >= minTrialsPerStimSess;
        mismatch(s,k) = nnz(gotCounts(use) ~= expCounts(use));
    end
end
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
end

function [keepLocal, info] = local_select_v4_sites(Tall_V4, R3full, RFrange)
nV4 = numel(RFrange);
siteRows = (1:nV4).';
R3 = localize_response_rows_local(R3full, RFrange);

nObjectStim = zeros(nV4,1);
for stim = 1:numel(Tall_V4)
    T = Tall_V4(stim).T;
    assign = string(T.assignment(siteRows));
    nObjectStim = nObjectStim + (assign == "target") + (assign == "distractor");
end
hasObjectRF = nObjectStim >= 1;

SNR = compute_snr_per_color_sites(R3, Tall_V4, siteRows, 'Verbose', false);
SNRmat = [SNR.yellowEarly(siteRows), SNR.yellowLate(siteRows), ...
          SNR.purpleEarly(siteRows), SNR.purpleLate(siteRows)];
bestSNR = max(SNRmat, [], 2, 'omitnan');
keepMask = hasObjectRF & isfinite(bestSNR) & (bestSNR > 0.7);
keepLocal = find(keepMask);

info = struct();
info.bestSNR = bestSNR;
info.hasObjectRF = hasObjectRF;
end

function [keepLocal, info] = local_select_it_sites(Tall_IT, R3full, RFrange)
nIT = numel(RFrange);
siteRows = (1:nIT).';
Rloc = localize_response_rows_local(R3full, RFrange);

nObjectStim = zeros(nIT,1);
for stim = 1:numel(Tall_IT)
    T = Tall_IT(stim).T;
    assign = string(T.assignment(siteRows));
    nObjectStim = nObjectStim + (assign == "target") + (assign == "distractor");
end
hasObjectRF = nObjectStim >= 1;

SNR = compute_snr_per_color_sites(Rloc, Tall_IT, siteRows, 'Verbose', false);
topKAbsQuartetEarly = local_compute_topk_quartet_response(Rloc, SNR, 5);
keepMask = hasObjectRF & isfinite(topKAbsQuartetEarly) & (topKAbsQuartetEarly > 0.7);
keepLocal = find(keepMask);

info = struct();
info.topKAbsQuartetEarly = topKAbsQuartetEarly;
info.hasObjectRF = hasObjectRF;
end

function topKAbsQuartetEarly = local_compute_topk_quartet_response(Rloc, SNR, topK)
nIT = size(Rloc.meanAct, 1);
nStim = size(Rloc.meanAct, 2);
muSpont = SNR.muSpont(:);
sdSpont = SNR.sdSpont(:);
badNoise = ~isfinite(sdSpont) | (sdSpont <= 0);

nQuartets = floor(nStim / 8) * 2;
quartetMembers = zeros(nQuartets, 4);
q = 0;
for base = 0:8:(nStim - 8)
    q = q + 1;
    quartetMembers(q,:) = base + [1 2 5 6];
    q = q + 1;
    quartetMembers(q,:) = base + [3 4 7 8];
end

if isvector(Rloc.nTrials)
    nTrialsAll = double(Rloc.nTrials(:));
    perSiteTrials = false;
else
    perSiteTrials = true;
end

muQuartetEarly = nan(nIT, nQuartets);
for iSite = 1:nIT
    if perSiteTrials
        nTrSite = double(Rloc.nTrials(iSite, :)).';
    else
        nTrSite = nTrialsAll;
    end
    rEarly = squeeze(Rloc.meanAct(iSite,:,2)).';
    for qIdx = 1:nQuartets
        stimQ = quartetMembers(qIdx,:);
        nTrQ = nTrSite(stimQ);
        rQ = rEarly(stimQ);
        idx = isfinite(rQ) & isfinite(nTrQ) & (nTrQ > 0);
        if any(idx)
            muQuartetEarly(iSite, qIdx) = sum(nTrQ(idx) .* rQ(idx)) / sum(nTrQ(idx));
        end
    end
end

signedQuartetEarly = bsxfun(@rdivide, bsxfun(@minus, muQuartetEarly, muSpont), sdSpont);
signedQuartetEarly(badNoise, :) = NaN;

topKAbsQuartetEarly = nan(nIT,1);
for iSite = 1:nIT
    vals = abs(signedQuartetEarly(iSite,:));
    vals = vals(isfinite(vals));
    if isempty(vals)
        continue;
    end
    vals = sort(vals, 'descend');
    k = min(topK, numel(vals));
    topKAbsQuartetEarly(iSite) = mean(vals(1:k));
end
end

function R_loc = localize_response_rows_local(R_full, siteGlobal)
R_loc = R_full;
R_loc.meanAct = R_full.meanAct(siteGlobal, :, :);
R_loc.meanSqAct = R_full.meanSqAct(siteGlobal, :, :);
if ismatrix(R_full.nTrials) && size(R_full.nTrials,1) >= max(siteGlobal)
    R_loc.nTrials = R_full.nTrials(siteGlobal, :);
else
    R_loc.nTrials = R_full.nTrials;
end
end

function TallSorted = local_sort_tall(Tall, tallFieldName)
stimNums = arrayfun(@(x) x.stimNum, Tall(:));
[stimNumsSorted, ordStim] = sort(stimNums(:));
assert(all(stimNumsSorted(:).' == 1:numel(Tall)), ...
    '%s.stimNum must cover 1..%d exactly.', tallFieldName, numel(Tall));
TallSorted = Tall(ordStim);
end

function keep = local_apply_limit(keep, siteLimitByArea, areaName)
keep = keep(:);
if ~isstruct(siteLimitByArea) || ~isfield(siteLimitByArea, areaName)
    return;
end
lim = siteLimitByArea.(areaName);
if isempty(lim) || ~isfinite(lim)
    return;
end
lim = max(0, floor(double(lim)));
if numel(keep) > lim
    keep = keep(1:lim);
end
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

function cacheTag = local_cache_tag(P)
if isfield(P, 'cacheTag') && strlength(string(P.cacheTag)) > 0
    cacheTag = char(string(P.cacheTag));
    return;
end

lims = P.siteLimitByArea;
if isempty(lims) || (~local_has_limit(lims, 'V1') && ~local_has_limit(lims, 'V4') && ~local_has_limit(lims, 'IT'))
    cacheTag = 'all';
else
    cacheTag = sprintf('v1%s_v4%s_it%s', local_limit_str(lims, 'V1'), local_limit_str(lims, 'V4'), local_limit_str(lims, 'IT'));
end
end

function tf = local_has_limit(S, fieldName)
tf = isstruct(S) && isfield(S, fieldName) && ~isempty(S.(fieldName)) && isfinite(S.(fieldName));
end

function s = local_limit_str(S, fieldName)
if local_has_limit(S, fieldName)
    s = sprintf('%d', floor(double(S.(fieldName))));
else
    s = 'all';
end
end

function S = local_cache_params(Monkey, P, cacheTag)
S = struct();
S.cacheVersion = 1;
S.Monkey = double(Monkey);
S.timeWindows = double(P.timeWindows);
S.onlyCorrect = logical(P.onlyCorrect);
S.correctCol = double(P.correctCol);
S.correctVal = double(P.correctVal);
S.stimCol = double(P.stimCol);
S.dayCol = double(P.dayCol);
S.sessions = double(P.sessions(:)');
S.minTrialsPerStimSess = double(P.minTrialsPerStimSess);
S.chunkTrials = double(P.chunkTrials);
S.cacheTag = string(cacheTag);
if isstruct(P.siteLimitByArea)
    S.siteLimitByArea = P.siteLimitByArea;
else
    S.siteLimitByArea = struct('V1', [], 'V4', [], 'IT', []);
end
end

function tf = local_cache_matches(OUT, cacheParams)
tf = isstruct(OUT) && isfield(OUT, 'cacheParams') && isequaln(OUT.cacheParams, cacheParams);
end

function tf = local_pair_files_exist(OUT)
tf = isstruct(OUT) && isfield(OUT, 'PairClassTable') && istable(OUT.PairClassTable);
if ~tf
    return;
end
if ~ismember('pairFile', OUT.PairClassTable.Properties.VariableNames)
    tf = false;
    return;
end
pairFiles = string(OUT.PairClassTable.pairFile);
tf = all(arrayfun(@(x) exist(char(x), 'file') == 2, pairFiles));
end

function tok = local_pair_class_token(cls)
tok = char(strrep(string(cls), '-', '_'));
end

function local_print_progress(label, done, total, elapsedSec)
frac = done / max(total, 1);
if done < max(10, ceil(0.02 * max(total, 1)))
    fprintf('%s %d / %d | elapsed %.1fs | ETA n/a\n', label, done, total, elapsedSec);
else
    etaSec = elapsedSec * (max(total, 1) - done) / max(done, 1);
    fprintf('%s %d / %d | elapsed %.1fs | ETA %.1fs\n', label, done, total, elapsedSec, etaSec);
end
end

function local_plot_summary(OUT)
pairCorr = OUT.pairCorr;
PairClassTable = OUT.PairClassTable;
tw = OUT.timeWindows;
colors = [0.35 0.35 0.35; 0.15 0.45 0.85; 0.85 0.40 0.15];

figure('Color', 'w', 'Name', sprintf('Noise correlations %s', char(OUT.monkeySuffix)));
tiledlayout(height(PairClassTable), 1, 'Padding', 'compact', 'TileSpacing', 'compact');
for i = 1:height(PairClassTable)
    ax = nexttile;
    hold(ax, 'on');
    s = PairClassTable.startIdx(i);
    e = PairClassTable.endIdx(i);
    if ~isfinite(s) || ~isfinite(e)
        title(ax, sprintf('%s (0 pairs)', PairClassTable.pairClass(i)));
        continue;
    end
    vals = pairCorr(s:e, :);
    for w = 1:size(vals, 2)
        histogram(ax, vals(:,w), 'DisplayStyle', 'stairs', 'EdgeColor', colors(w,:), 'LineWidth', 1.5);
    end
    title(ax, sprintf('%s (N=%d pairs)', PairClassTable.pairClass(i), PairClassTable.nPairs(i)));
    xlabel(ax, 'Noise correlation');
    ylabel(ax, 'N pairs');
    grid(ax, 'on');
end
legend(compose('%d to %d ms', tw(:,1), tw(:,2)), 'Location', 'best');
end
