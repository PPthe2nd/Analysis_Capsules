function OUT = Attention_TargetSide_Timecourse_IT(Puser)
% ATTENTION_TARGETSIDE_TIMECOURSE_IT
% Time-resolved IT responses for late target-direction and target-distance
% site classes, using the high-resolution 10 ms response file
% Resp_capsules_*_d12.mat.
%
% Four site groups are shown:
%   1) distance-growth-significant sites   (late distance effect > prestim)
%   2) direction-significant sites         (dir|dist)
%   3) distance-growth-only sites          (late distance effect > prestim & ~dir|dist)
%   4) direction-only sites                (dir|dist & ~dist|dir)
%
% Within each site, direction traces are aligned by the fitted late
% direction preference from Attention_TargetSide_DistanceControl_IT.
%
% Distance traces are NOT aligned by activity. They are defined purely from
% the independently measured RF geometry:
%   - target-closer configuration     : target arm closer to RF than the
%                                       distractor arm for that quartet
%   - distractor-closer configuration : opposite signed geometry

%% Settings
P = struct();
P.Monkey = 1; % 1 = Nilson, 2 = Figaro
P.respCapsulesFileN = 'Resp_capsules_N_d12.mat'; % 10 ms bins
P.respCapsulesFileF = 'Resp_capsules_F_d12.mat';
P.resp3binFileN = 'SNR_capsules_N_d12.mat';
P.resp3binFileF = 'SNR_capsules_F_d12.mat';
P.normalizeToSpontSd = true;
P.minQuartetsPerSite = 20;
P.plotSem = true;
P.plotFigure = true;
P.plotDifferenceFigure = true;
P.saveResult = true;
P.forceRefit = false;
P.fitTiming = true;
P.plotTimingFitFigure = false;
P.forceTimingRefit = false;
P.timingFitStartMs = 0;
P.timingFitEndMs = 500;
P.timingFitSmoothW = 3;
P.timingFitMinPoints = 6;
P.timingFitMinTau = 5;
P.timingFitMaxTau = 400;
P.timingFitInitTau = 40;
P.timingFitT50PadMs = 100;
P.timingFitMaxAmpFactor = 3;
P.timingFitVerbose = true;
P.diffAlpha = 0.05;
P.diffTest = 'signrank'; % 'signrank' or 'ttest'
P.diffFdr = true;
P.diffWindowMs = [100 400];
P.distSelectPreWindowMs = [-200 0];
P.distSelectLateWindowMs = [300 500];
P.distSelectAlpha = 0.05;
P.distSelectTest = 'signrank'; % 'signrank' or 'ttest'

if nargin >= 1 && ~isempty(Puser)
    fn = fieldnames(Puser);
    for i = 1:numel(fn)
        P.(fn{i}) = Puser.(fn{i});
    end
end

cfg = config();

%% Monkey-specific files
if P.Monkey == 1
    monkeySuffix = "N";
    tallFile = 'Tall_IT_lines_N.mat';
    respFile = P.respCapsulesFileN;
    resp3binFile = P.resp3binFileN;
elseif P.Monkey == 2
    monkeySuffix = "F";
    tallFile = 'Tall_IT_lines_F.mat';
    respFile = P.respCapsulesFileF;
    resp3binFile = P.resp3binFileF;
else
    error('Attention_TargetSide_Timecourse_IT:InvalidMonkey', ...
        'P.Monkey must be 1 (Nilson) or 2 (Figaro).');
end

tallPath = fullfile(cfg.matDir, tallFile);
respPath = fullfile(cfg.matDir, respFile);
resp3binPath = fullfile(cfg.matDir, resp3binFile);
basePath = fullfile(cfg.matDir, sprintf('Attention_TargetSide_Tuning_IT_directiondelta_%s.mat', char(monkeySuffix)));
distPath = fullfile(cfg.matDir, sprintf('Attention_TargetSide_DistanceControl_IT_%s.mat', char(monkeySuffix)));
outPath = fullfile(cfg.resultsDir, sprintf('Attention_TargetSide_Timecourse_IT_%s_d12.mat', char(monkeySuffix)));
currentExclusions = site_session_exclusions(monkeySuffix);
hasSessionExclusions = ~isempty(currentExclusions);

assert(exist(tallPath, 'file') == 2, 'Missing %s.', tallPath);
assert(exist(respPath, 'file') == 2, 'Missing %s.', respPath);
assert(exist(resp3binPath, 'file') == 2, 'Missing %s.', resp3binPath);

useCached = exist(outPath, 'file') == 2 && ~P.forceRefit;
if useCached
    S = load(outPath, 'OUT');
    if session_exclusion_cache_matches(S, monkeySuffix)
        assert(isfield(S, 'OUT') && isstruct(S.OUT), '%s must contain OUT.', outPath);
        OUT = S.OUT;
        cacheUpdated = false;
        if P.fitTiming
            if ~isfield(OUT, 'TimingFit') || ~isstruct(OUT.TimingFit) || P.forceTimingRefit
                OUT.TimingFit = fit_timing_sigmoids_local(OUT.Groups, OUT.timeWindows, P);
                cacheUpdated = true;
            end
            if P.plotTimingFitFigure
                plot_timing_fit_figure_local(OUT.TimingFit, char(OUT.monkeySuffix));
            end
        end
        if cacheUpdated && P.saveResult
            save(outPath, 'OUT', '-v7.3');
        end
        if P.plotFigure
            make_timecourse_figure(OUT);
        end
        if P.plotDifferenceFigure
            make_difference_figure(OUT);
        end
        return;
    end
    if hasSessionExclusions
        fprintf(['Cached IT target-side timecourse does not match the active session exclusions ' ...
                 'for monkey %s; recomputing.\n'], char(monkeySuffix));
    else
        fprintf('Cached IT target-side timecourse is from an exclusion-aware run; recomputing canonical cache.\n');
    end
end
%% Load geometry and cached analyses
Sgeo = load(tallPath, 'Tall_IT', 'RFrange');
assert(isfield(Sgeo, 'Tall_IT') && isstruct(Sgeo.Tall_IT), '%s must contain Tall_IT.', tallPath);
assert(isfield(Sgeo, 'RFrange') && ~isempty(Sgeo.RFrange), '%s must contain RFrange.', tallPath);
Tall_IT = Sgeo.Tall_IT;
RFrange = Sgeo.RFrange(:);
nIT = numel(RFrange);
itSites = (1:nIT).';

if hasSessionExclusions
    fprintf(['Session exclusions are active for monkey %s; refreshing IT target-side ' ...
             'fits before timecourse extraction.\n'], char(monkeySuffix));
    BASE = Attention_TargetSide_Tuning_IT(struct( ...
        'makeSummaryFigures', false, ...
        'makeExampleFigures', false, ...
        'makeAngleReferenceFigure', false));
    DIST = Attention_TargetSide_DistanceControl_IT(struct('makeSummaryFigure', false));
else
    assert(exist(basePath, 'file') == 2, 'Missing %s. Run Attention_TargetSide_Tuning_IT first.', basePath);
    assert(exist(distPath, 'file') == 2, 'Missing %s. Run Attention_TargetSide_DistanceControl_IT first.', distPath);
    Sbase = load(basePath, 'OUT');
    Sdist = load(distPath, 'OUT');
    BASE = Sbase.OUT;
    DIST = Sdist.OUT;
end

reqBase = {'QuartetTable','deltaQuartetLate','varQuartetLate','RFrange'};
for k = 1:numel(reqBase)
    assert(isfield(BASE, reqBase{k}), 'Base OUT missing field %s.', reqBase{k});
end
reqDist = {'RegressionLate','ctrlSig','distGivenDirSig','targetAdvPx'};
for k = 1:numel(reqDist)
    assert(isfield(DIST, reqDist{k}), 'Distance OUT missing field %s.', reqDist{k});
end
assert(numel(BASE.RFrange) == nIT && all(BASE.RFrange(:) == RFrange(:)), ...
    'Base RFrange does not match Tall_IT RFrange.');

QuartetTable = BASE.QuartetTable;
thetaDeg = double(QuartetTable.targetDirDeg(:));
stimRef = double(QuartetTable.stimRef(:));
nQuartets = height(QuartetTable);
pairA = [stimRef, stimRef + 4];
pairB = [stimRef + 1, stimRef + 5];

%% Load 3-bin responses for spontaneous normalization
R3full = load_capsules_struct_exclusion_aware(resp3binPath, monkeySuffix, 'cfg', cfg);
R3 = R3full;
R3.meanAct = R3full.meanAct(RFrange, :, :);
R3.meanSqAct = R3full.meanSqAct(RFrange, :, :);
if ismatrix(R3full.nTrials) && size(R3full.nTrials,1) >= max(RFrange)
    R3.nTrials = R3full.nTrials(RFrange, :);
else
    R3.nTrials = R3full.nTrials;
end
SNR = compute_snr_per_color_sites(R3, Tall_IT, itSites, 'Verbose', false);
muSpont = double(SNR.muSpont(itSites));
sdSpont = double(SNR.sdSpont(itSites));

%% Load high-resolution responses
Rfull = load_capsules_struct_exclusion_aware(respPath, monkeySuffix, 'cfg', cfg);
Rresp = Rfull;
Rresp.meanAct = Rfull.meanAct(RFrange, :, :);
if isfield(Rfull, 'meanSqAct')
    Rresp.meanSqAct = Rfull.meanSqAct(RFrange, :, :);
end
if ismatrix(Rfull.nTrials) && size(Rfull.nTrials,1) >= max(RFrange)
    Rresp.nTrials = Rfull.nTrials(RFrange, :);
else
    Rresp.nTrials = Rfull.nTrials;
end

[nCh, nStim, nBins] = size(Rresp.meanAct);
assert(nCh == nIT, 'Expected %d localized IT rows, got %d.', nIT, nCh);
assert(nStim == 384, 'Expected 384 stimuli in %s.', respFile);
assert(size(Rresp.timeWindows,1) == nBins && size(Rresp.timeWindows,2) == 2, ...
    'Rresp.timeWindows must be [nBins x 2].');
tCenters = mean(double(Rresp.timeWindows), 2);
distPreMask = (tCenters >= P.distSelectPreWindowMs(1)) & (tCenters < P.distSelectPreWindowMs(2));
distLateMask = (tCenters >= P.distSelectLateWindowMs(1)) & (tCenters < P.distSelectLateWindowMs(2));
assert(any(distPreMask), 'No bins found in distance-selection prestim window [%d %d).', ...
    P.distSelectPreWindowMs(1), P.distSelectPreWindowMs(2));
assert(any(distLateMask), 'No bins found in distance-selection late window [%d %d).', ...
    P.distSelectLateWindowMs(1), P.distSelectLateWindowMs(2));

if isvector(Rresp.nTrials)
    nTrialsShared = double(Rresp.nTrials(:)');
    perSiteTrials = false;
elseif ismatrix(Rresp.nTrials) && size(Rresp.nTrials,2) == nStim
    perSiteTrials = true;
    nTrialsShared = [];
else
    error('Rresp.nTrials must be vector(384) or matrix(nChannels x 384).');
end

%% Site-wise direction- and geometry-aligned time courses
tcDirPref = nan(nIT, nBins);
tcDirNon = nan(nIT, nBins);
tcDistTar = nan(nIT, nBins);
tcDistDis = nan(nIT, nBins);
nQuartetDir = zeros(nIT,1);
nQuartetDist = zeros(nIT,1);
distPreMean = nan(nIT,1);
distLateMean = nan(nIT,1);
distLateMinusPreMean = nan(nIT,1);
distLateMinusPreP = nan(nIT,1);
distLateMinusPreN = zeros(nIT,1);

targetAdvPx = double(DIST.targetAdvPx);
deltaLate = double(BASE.deltaQuartetLate);
varLate = double(BASE.varQuartetLate);
Reg = DIST.RegressionLate;

prefDegDir = nan(nIT,1);
betaTargetAdv = nan(nIT,1);
for s = 1:nIT
    if isfield(Reg(s), 'betaCos') && isfield(Reg(s), 'betaSin') && ...
            isfinite(Reg(s).betaCos) && isfinite(Reg(s).betaSin)
        prefDegDir(s) = mod(atan2d(Reg(s).betaSin, Reg(s).betaCos), 360);
    end
    if isfield(Reg(s), 'betaTargetAdv')
        betaTargetAdv(s) = double(Reg(s).betaTargetAdv);
    end
end

for s = 1:nIT
    validQ = isfinite(deltaLate(s,:)) & isfinite(varLate(s,:)) & isfinite(targetAdvPx(s,:));
    if ~any(validQ)
        continue;
    end

    resp = squeeze(double(Rresp.meanAct(s,:,:))); % [nStim x nBins]
    if size(resp,1) ~= nStim && size(resp,2) == nStim
        resp = resp.';
    end
    assert(size(resp,1) == nStim && size(resp,2) == nBins, ...
        'Unexpected response shape for IT site %d.', s);

    if perSiteTrials
        nTrSite = double(Rresp.nTrials(s,:));
    else
        nTrSite = nTrialsShared;
    end
    nTrSite(~isfinite(nTrSite) | nTrSite < 0) = 0;

    dirPrefByQ = nan(nQuartets, nBins);
    dirNonByQ = nan(nQuartets, nBins);
    distTarByQ = nan(nQuartets, nBins);
    distDisByQ = nan(nQuartets, nBins);

    for q = find(validQ)
        tcA = pair_mean_timecourse(resp, nTrSite, pairA(q,:));
        tcB = pair_mean_timecourse(resp, nTrSite, pairB(q,:));
        if ~any(isfinite(tcA)) || ~any(isfinite(tcB))
            continue;
        end

        if P.normalizeToSpontSd
            tcA = normalize_spont_sd(tcA, muSpont(s), sdSpont(s));
            tcB = normalize_spont_sd(tcB, muSpont(s), sdSpont(s));
        end

        if isfinite(prefDegDir(s))
            dA = abs(local_angdiff_deg(thetaDeg(q), prefDegDir(s)));
            dB = abs(local_angdiff_deg(mod(thetaDeg(q) + 180, 360), prefDegDir(s)));
            if dA <= dB
                dirPrefByQ(q,:) = tcA;
                dirNonByQ(q,:) = tcB;
            else
                dirPrefByQ(q,:) = tcB;
                dirNonByQ(q,:) = tcA;
            end
        end

        if isfinite(targetAdvPx(s,q)) && (targetAdvPx(s,q) ~= 0)
            if targetAdvPx(s,q) >= 0
                distTarByQ(q,:) = tcA;
                distDisByQ(q,:) = tcB;
            else
                distTarByQ(q,:) = tcB;
                distDisByQ(q,:) = tcA;
            end
        end
    end

    goodDir = all(isfinite(dirPrefByQ), 2) & all(isfinite(dirNonByQ), 2);
    goodDist = all(isfinite(distTarByQ), 2) & all(isfinite(distDisByQ), 2);

    nQuartetDir(s) = nnz(goodDir);
    nQuartetDist(s) = nnz(goodDist);

    if nQuartetDir(s) >= P.minQuartetsPerSite
        tcDirPref(s,:) = mean(dirPrefByQ(goodDir,:), 1, 'omitnan');
        tcDirNon(s,:) = mean(dirNonByQ(goodDir,:), 1, 'omitnan');
    end
    if nQuartetDist(s) >= P.minQuartetsPerSite
        tcDistTar(s,:) = mean(distTarByQ(goodDist,:), 1, 'omitnan');
        tcDistDis(s,:) = mean(distDisByQ(goodDist,:), 1, 'omitnan');

        dPre = mean(distTarByQ(goodDist, distPreMask), 2, 'omitnan') - ...
            mean(distDisByQ(goodDist, distPreMask), 2, 'omitnan');
        dLate = mean(distTarByQ(goodDist, distLateMask), 2, 'omitnan') - ...
            mean(distDisByQ(goodDist, distLateMask), 2, 'omitnan');
        dDelta = dLate - dPre;
        okDelta = isfinite(dPre) & isfinite(dLate) & isfinite(dDelta);
        if any(okDelta)
            dPre = dPre(okDelta);
            dLate = dLate(okDelta);
            dDelta = dDelta(okDelta);
            distPreMean(s) = mean(dPre, 'omitnan');
            distLateMean(s) = mean(dLate, 'omitnan');
            distLateMinusPreMean(s) = mean(dDelta, 'omitnan');
            distLateMinusPreP(s) = one_sample_p_local(dDelta, P.distSelectAlpha, P.distSelectTest);
            distLateMinusPreN(s) = numel(dDelta);
        end
    end
end

%% Site groups
isDirSig = logical(DIST.ctrlSig(:)) & (nQuartetDir >= P.minQuartetsPerSite);
isDistSig = isfinite(distLateMinusPreP) & (distLateMinusPreP < P.distSelectAlpha) & ...
    (distLateMinusPreMean > 0) & (nQuartetDist >= P.minQuartetsPerSite);
isDistOnly = isDistSig & ~isDirSig;
isDirOnly = isDirSig & ~isDistSig;

fprintf('IT target-side time courses using %s (%s)\n', respFile, char(monkeySuffix));
fprintf('  Distance-growth sites (late-pre): %d\n', nnz(isDistSig));
fprintf('  Direction-significant sites (dir|dist): %d\n', nnz(isDirSig));
fprintf('  Distance-growth-only sites: %d\n', nnz(isDistOnly));
fprintf('  Direction-only sites: %d\n', nnz(isDirOnly));

G1 = build_group_struct('Distance growth sig (late-pre)', isDistSig, tcDistTar, tcDistDis, ...
    'Target-closer config', 'Distractor-closer config');
G2 = build_group_struct('Direction sig (dir|dist)', isDirSig, tcDirPref, tcDirNon, ...
    'Preferred direction config', 'Opposite direction config');
G3 = build_group_struct('Distance growth only', isDistOnly, tcDistTar, tcDistDis, ...
    'Target-closer config', 'Distractor-closer config');
G4 = build_group_struct('Direction only', isDirOnly, tcDirPref, tcDirNon, ...
    'Preferred direction config', 'Opposite direction config');
G = [G1 G2 G3 G4];

%% Output
OUT = struct();
OUT.P = P;
OUT.monkeySuffix = monkeySuffix;
OUT.RFrange = RFrange;
OUT.tCenters = tCenters;
OUT.timeWindows = double(Rresp.timeWindows);
OUT.prefDegDir = prefDegDir;
OUT.betaTargetAdv = betaTargetAdv;
OUT.nQuartetDir = nQuartetDir;
OUT.nQuartetDist = nQuartetDist;
OUT.distSelectPreWindowMs = P.distSelectPreWindowMs;
OUT.distSelectLateWindowMs = P.distSelectLateWindowMs;
OUT.distPreMean = distPreMean;
OUT.distLateMean = distLateMean;
OUT.distLateMinusPreMean = distLateMinusPreMean;
OUT.distLateMinusPreP = distLateMinusPreP;
OUT.distLateMinusPreN = distLateMinusPreN;
OUT.tcDirPref = tcDirPref;
OUT.tcDirNon = tcDirNon;
OUT.tcDistTargetCloser = tcDistTar;
OUT.tcDistDistractorCloser = tcDistDis;
OUT.tcDistPref = tcDistTar;
OUT.tcDistNon = tcDistDis;
OUT.isDirSig = isDirSig;
OUT.isDistSig = isDistSig;
OUT.isDistOnly = isDistOnly;
OUT.isDirOnly = isDirOnly;
OUT.Groups = G;
OUT.siteSessionExclusions = currentExclusions;
if P.fitTiming
    OUT.TimingFit = fit_timing_sigmoids_local(G, OUT.timeWindows, P);
end

if P.saveResult
    save(outPath, 'OUT', '-v7.3');
    fprintf('Saved IT target-side time courses to %s\n', outPath);
end

if P.plotFigure
    make_timecourse_figure(OUT);
end
if P.plotDifferenceFigure
    make_difference_figure(OUT);
end
if P.fitTiming && P.plotTimingFitFigure
    plot_timing_fit_figure_local(OUT.TimingFit, char(OUT.monkeySuffix));
end
end

function G = build_group_struct(label, mask, tcCond1All, tcCond2All, cond1Label, cond2Label)
mask = logical(mask(:));
G = struct();
G.label = label;
G.mask = mask;
G.nSites = nnz(mask);
G.cond1Label = cond1Label;
G.cond2Label = cond2Label;
G.tcPref = tcCond1All(mask, :);
G.tcNon = tcCond2All(mask, :);
G.meanPref = mean(G.tcPref, 1, 'omitnan');
G.meanNon = mean(G.tcNon, 1, 'omitnan');
G.semPref = std(G.tcPref, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(G.tcPref),1));
G.semNon = std(G.tcNon, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(G.tcNon),1));
G.tcDiff = G.tcPref - G.tcNon;
G.meanDiff = mean(G.tcDiff, 1, 'omitnan');
G.semDiff = std(G.tcDiff, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(G.tcDiff),1));
end

function make_timecourse_figure(OUT)
t = OUT.tCenters(:);
G = OUT.Groups;

cPref = [0.85 0.20 0.15];
cNon = [0.25 0.25 0.25];

fig = figure('Color', 'w', 'Name', sprintf('IT target-side time courses (%s)', char(OUT.monkeySuffix)));
useTiled = exist('tiledlayout', 'file') == 2;
if useTiled
    tiledlayout(2,2, 'TileSpacing', 'compact', 'Padding', 'compact');
end

for i = 1:numel(G)
    if useTiled
        ax = nexttile;
    else
        ax = subplot(2,2,i); %#ok<LAXES>
    end
    hold(ax, 'on');

    if OUT.P.plotSem
        add_sem_patch(ax, t, G(i).meanPref, G(i).semPref, cPref, 0.18);
        add_sem_patch(ax, t, G(i).meanNon, G(i).semNon, cNon, 0.14);
    end

    plot(ax, t, G(i).meanPref, '-', 'Color', cPref, 'LineWidth', 2.2);
    plot(ax, t, G(i).meanNon, '--', 'Color', cNon, 'LineWidth', 2.0);
    xline(ax, 0, 'k-');
    xlabel(ax, 'Time from stimulus onset (ms)');
    ylabel(ax, 'Response (spont. SD units)');
    title(ax, sprintf('%s (N=%d)', G(i).label, G(i).nSites));
    grid(ax, 'on');
    if i == 1
        legend(ax, {G(i).cond1Label, G(i).cond2Label}, ...
            'Location', 'best', 'Box', 'off');
    end
end

annotation(fig, 'textbox', [0.08 0.94 0.84 0.05], ...
    'String', sprintf(['IT target-side time courses (%s, 10 ms bins). ' ...
    'Direction panels use the fitted direction preference; distance panels use target-closer vs distractor-closer configurations defined only by RF geometry, and the distance site set is selected by a larger 300-500 ms than prestim effect.'], ...
    char(OUT.monkeySuffix)), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'FontSize', 12, 'FontWeight', 'bold');
end

function make_difference_figure(OUT)
t = OUT.tCenters(:);
G = OUT.Groups;
cDiff = [0.10 0.10 0.10];
cPatch = [0.45 0.45 0.45];
fitCols = [ ...
    0.0000    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560];

fig = figure('Color', 'w', 'Name', sprintf('IT target-side paired differences (%s)', char(OUT.monkeySuffix)));
useTiled = exist('tiledlayout', 'file') == 2;
if useTiled
    tiledlayout(2,2, 'TileSpacing', 'compact', 'Padding', 'compact');
end

for i = 1:numel(G)
    if useTiled
        ax = nexttile;
    else
        ax = subplot(2,2,i); %#ok<LAXES>
    end
    hold(ax, 'on');

    [~, ~, sigBin, pWin] = diff_stats_local(G(i).tcDiff, t, OUT.P);
    if OUT.P.plotSem
        add_sem_patch(ax, t, G(i).meanDiff, G(i).semDiff, cPatch, 0.20);
    end
    hDiff = plot(ax, t, G(i).meanDiff, '-', 'Color', cDiff, 'LineWidth', 2.2);
    hFit = [];
    if isfield(OUT, 'TimingFit') && isstruct(OUT.TimingFit) && ...
            isfield(OUT.TimingFit, 'yFit') && size(OUT.TimingFit.yFit,1) >= i
        yFit = double(OUT.TimingFit.yFit(i,:));
        if any(isfinite(yFit))
            hFit = plot(ax, t, yFit, '-', 'Color', fitCols(i,:), 'LineWidth', 2.0);
        end
        if isfield(OUT.TimingFit, 'summary') && height(OUT.TimingFit.summary) >= i && ...
                isfinite(OUT.TimingFit.summary.t50(i))
            plot(ax, OUT.TimingFit.summary.t50(i), 0, 'o', 'MarkerSize', 5.5, ...
                'MarkerFaceColor', fitCols(i,:), 'MarkerEdgeColor', [0.1 0.1 0.1], 'LineWidth', 0.8);
        end
    end
    xline(ax, 0, 'k-');
    yline(ax, 0, 'k:');
    xlabel(ax, 'Time from stimulus onset (ms)');
    ylabel(ax, '\Delta response');
    title(ax, sprintf('%s (N=%d, win p=%.3g)', G(i).label, G(i).nSites, pWin));
    grid(ax, 'on');

    yL = ylim(ax);
    ySig = yL(1) + 0.08 * range(yL);
    if any(sigBin)
        plot(ax, t(sigBin), ySig * ones(nnz(sigBin),1), 'o', ...
            'MarkerSize', 3.8, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none');
    end
    text(ax, t(2), ySig, sig_label_local(OUT.P), ...
        'FontSize', 8, 'VerticalAlignment', 'bottom', 'Color', [0.2 0.2 0.2]);
    if isfield(OUT, 'TimingFit') && isstruct(OUT.TimingFit) && isfield(OUT.TimingFit, 'summary') && ...
            height(OUT.TimingFit.summary) >= i && isfinite(OUT.TimingFit.summary.t50(i))
        xText = t(end) - 0.02 * range(t);
        yText = yL(1) + 0.16 * range(yL);
        text(ax, xText, yText, sprintf('t_{50} = %.0f ms', OUT.TimingFit.summary.t50(i)), ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
            'FontSize', 10, 'FontWeight', 'bold', 'Color', fitCols(i,:), ...
            'BackgroundColor', 'w', 'Margin', 1.5);
    end
    if i == 1 && ~isempty(hFit) && isgraphics(hFit)
        legend(ax, [hDiff hFit], {'mean diff','sigmoid fit'}, 'Location', 'northwest', 'Box', 'off');
    end
end

annotation(fig, 'textbox', [0.05 0.94 0.90 0.05], ...
    'String', sprintf(['IT paired difference traces (%s, 10 ms bins). ' ...
    'Direction panels show preferred-direction minus opposite-direction; distance panels show target-closer minus distractor-closer, defined only by RF geometry, with distance sites selected by a larger 300-500 ms than prestim effect. Colored line = independent sigmoid fit; t_{50} is shown in each panel.'], ...
    char(OUT.monkeySuffix)), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'FontSize', 12, 'FontWeight', 'bold');
end

function tc = pair_mean_timecourse(respStimByTime, nTrialsStim, stimSet)
stimIdx = stimSet(:)';
nTrials = double(nTrialsStim(:))';
mask = ismember(1:size(respStimByTime,1), stimIdx);
mask = mask & isfinite(nTrials) & (nTrials > 0);

nBins = size(respStimByTime, 2);
tc = nan(1, nBins);
if ~any(mask)
    return;
end

w = nTrials(mask);
X = respStimByTime(mask, :);
W = repmat(w(:), 1, nBins);
ok = isfinite(X);
num = sum((X .* W) .* ok, 1, 'omitnan');
den = sum(W .* ok, 1, 'omitnan');
use = den > 0;
tc(use) = num(use) ./ den(use);
end

function tc = normalize_spont_sd(tc, muSpont, sdSpont)
if ~(isfinite(muSpont) && isfinite(sdSpont) && sdSpont > 0)
    tc(:) = NaN;
else
    tc = (tc - muSpont) ./ sdSpont;
end
end

function add_sem_patch(ax, t, m, s, col, faceAlpha)
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
    'FaceAlpha', faceAlpha, 'EdgeColor', 'none');
end

function F = fit_timing_sigmoids_local(Groups, timeWindows, P)
nGroups = numel(Groups);
nBins = numel(mean(double(timeWindows), 2));
t = mean(double(timeWindows), 2);
t = t(:);
fitMaskBase = isfinite(t) & (t >= P.timingFitStartMs) & (t <= P.timingFitEndMs);

groupIdx = (1:nGroups).';
groupLabel = strings(nGroups,1);
Ahat = nan(nGroups,1);
t50hat = nan(nGroups,1);
tauhat = nan(nGroups,1);
rmse = nan(nGroups,1);
nFitPts = zeros(nGroups,1);
residSD = nan(nGroups,1);
slopeAtT50 = nan(nGroups,1);
t50SE = nan(nGroups,1);
status = strings(nGroups,1);
status(:) = "not_fitted";
yRaw = nan(nGroups, nBins);
ySmooth = nan(nGroups, nBins);
yFit = nan(nGroups, nBins);

for g = 1:nGroups
    groupLabel(g) = string(Groups(g).label);
    ygRaw = double(Groups(g).meanDiff(:));
    if numel(ygRaw) ~= nBins
        error('Timing fit expects %d bins for group %d, got %d.', nBins, g, numel(ygRaw));
    end
    yRaw(g,:) = ygRaw(:).';
    yg = ygRaw;
    if P.timingFitSmoothW > 1
        yg = smooth_vec_movmean_omitnan_local(yg, round(P.timingFitSmoothW));
    end
    ySmooth(g,:) = yg(:).';

    fitMask = fitMaskBase & isfinite(yg);
    nFitPts(g) = nnz(fitMask);
    if nFitPts(g) < P.timingFitMinPoints
        status(g) = "too_few_points";
        continue;
    end

    tFit = t(fitMask);
    yFitData = yg(fitMask);
    A0 = max(yFitData, [], 'omitnan');
    if ~isfinite(A0) || A0 <= 0
        status(g) = "nonpositive";
        continue;
    end

    yHalf = 0.5 * A0;
    iHalf = find(yFitData >= yHalf, 1, 'first');
    if isempty(iHalf)
        t500 = median(tFit, 'omitnan');
    else
        t500 = tFit(iHalf);
    end

    AGlobal = max(yFitData, [], 'omitnan');
    if ~isfinite(AGlobal) || AGlobal <= 0
        AGlobal = A0;
    end
    AUpper = max(1e-3, P.timingFitMaxAmpFactor * AGlobal);
    p0 = [A0; t500; min(max(P.timingFitInitTau, P.timingFitMinTau), P.timingFitMaxTau)];
    lb = [0; min(tFit) - P.timingFitT50PadMs; P.timingFitMinTau];
    ub = [AUpper; max(tFit) + P.timingFitT50PadMs; P.timingFitMaxTau];

    optsLSQ = optimoptions('lsqcurvefit', 'Display', 'off');
    modelFun = @(pp,xx) pp(1) ./ (1 + exp(-(xx - pp(2)) ./ pp(3)));
    try
        pHat = lsqcurvefit(modelFun, p0, tFit, yFitData, lb, ub, optsLSQ);
    catch
        status(g) = "fit_failed";
        continue;
    end

    Ahat(g) = pHat(1);
    t50hat(g) = pHat(2);
    tauhat(g) = pHat(3);
    yFit(g,:) = modelFun(pHat, t(:)).';
    e = yFit(g,fitMask).' - yFitData;
    rmse(g) = sqrt(mean(e.^2, 'omitnan'));
    dof = nFitPts(g) - 2;
    if dof > 0
        rss = sum(e.^2, 'omitnan');
        if isfinite(rss) && rss >= 0
            residSD(g) = sqrt(rss / dof);
        end
    end
    if isfinite(Ahat(g)) && isfinite(tauhat(g)) && tauhat(g) > 0
        slopeAtT50(g) = abs(Ahat(g)) / (4 * tauhat(g));
    end
    if isfinite(residSD(g)) && isfinite(slopeAtT50(g)) && slopeAtT50(g) > 0
        t50SE(g) = residSD(g) / slopeAtT50(g);
    end
    status(g) = "ok";
end

summary = table(groupIdx, groupLabel, Ahat, t50hat, tauhat, rmse, nFitPts, ...
    residSD, slopeAtT50, t50SE, status, ...
    'VariableNames', {'groupIdx','groupLabel','A','t50','tau','rmse','nFitPoints', ...
    'residSD','slopeAtT50','t50SE','status'});

F = struct();
F.summary = summary;
F.t = t;
F.yRaw = yRaw;
F.ySmooth = ySmooth;
F.yFit = yFit;
F.fitMask = fitMaskBase;
F.options = P;

if P.timingFitVerbose
    fprintf('Independent sigmoid timing fits for IT target-side effects:\n');
    disp(summary(:, {'groupLabel','A','t50','t50SE','tau','rmse','nFitPoints','status'}));
end
end

function plot_timing_fit_figure_local(F, monkeySuffix)
if isempty(F) || ~isstruct(F) || ~isfield(F, 'summary')
    return;
end

t = F.t(:);
S = F.summary;
nGroups = height(S);
fig = figure('Color', 'w', 'Name', sprintf('IT target-side timing fits (%s)', monkeySuffix), ...
    'NumberTitle', 'off');
useTiled = exist('tiledlayout', 'file') == 2;
if useTiled
    tiledlayout(2,2, 'TileSpacing', 'compact', 'Padding', 'compact');
end

cmap = lines(max(nGroups, 4));
for g = 1:nGroups
    if useTiled
        ax = nexttile;
    else
        ax = subplot(2,2,g); %#ok<LAXES>
    end
    hold(ax, 'on');
    rawCol = 0.65 * [1 1 1] + 0.35 * cmap(g,:);
    plot(ax, t, F.yRaw(g,:), '-', 'Color', rawCol, 'LineWidth', 1.1);
    plot(ax, t, F.ySmooth(g,:), '-', 'Color', cmap(g,:), 'LineWidth', 1.8);
    if all(isfinite(F.yFit(g,:)))
        plot(ax, t, F.yFit(g,:), '-', 'Color', [0.1 0.1 0.1], 'LineWidth', 2.1);
    end
    if isfinite(S.t50(g))
        plot(ax, S.t50(g), 0, 'o', 'MarkerSize', 6.5, ...
            'MarkerFaceColor', cmap(g,:), 'MarkerEdgeColor', [0.1 0.1 0.1], 'LineWidth', 0.8);
    end
    xline(ax, 0, 'k--');
    yline(ax, 0, 'k:');
    grid(ax, 'on');
    box(ax, 'off');
    xlabel(ax, 'Time from stimulus onset (ms)');
    ylabel(ax, '\Delta response');
    title(ax, sprintf('%s | t50=%s ms | tau=%s ms', ...
        char(S.groupLabel(g)), num2str_or_na_local(S.t50(g), '%.1f'), num2str_or_na_local(S.tau(g), '%.1f')));
end

annotation(fig, 'textbox', [0.06 0.94 0.88 0.05], ...
    'String', sprintf('Independent sigmoid timing fits for IT target-side effect traces (%s). Gray=raw, color=smoothed, black=sigmoid, circle=t50.', monkeySuffix), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontSize', 12, 'FontWeight', 'bold');
end

function y = smooth_vec_movmean_omitnan_local(x, w)
x = double(x(:));
n = numel(x);
w = max(1, round(w));
if w <= 1 || n == 0
    y = x;
    return;
end
halfLo = floor((w-1)/2);
halfHi = ceil((w-1)/2);
y = nan(size(x));
for i = 1:n
    i1 = max(1, i-halfLo);
    i2 = min(n, i+halfHi);
    xi = x(i1:i2);
    if any(isfinite(xi))
        y(i) = mean(xi, 'omitnan');
    end
end
end

function s = num2str_or_na_local(x, fmt)
if nargin < 2 || isempty(fmt)
    fmt = '%.3g';
end
if isfinite(x)
    s = sprintf(fmt, x);
else
    s = 'NA';
end
end

function [pBin, qBin, sigBin, pWin] = diff_stats_local(tcDiff, tCenters, P)
nBins = size(tcDiff, 2);
pBin = nan(1, nBins);
for b = 1:nBins
    xb = tcDiff(:, b);
    xb = xb(isfinite(xb));
    if numel(xb) < 3
        continue;
    end
    pBin(b) = one_sample_p_local(xb, P.diffAlpha, P.diffTest);
end

if P.diffFdr
    [qBin, sigBin] = fdr_bh_local(pBin, P.diffAlpha);
else
    qBin = nan(size(pBin));
    sigBin = isfinite(pBin) & (pBin < P.diffAlpha);
end

wMask = (tCenters >= P.diffWindowMs(1)) & (tCenters <= P.diffWindowMs(2));
xw = mean(tcDiff(:, wMask), 2, 'omitnan');
xw = xw(isfinite(xw));
if numel(xw) >= 3
    pWin = one_sample_p_local(xw, P.diffAlpha, P.diffTest);
else
    pWin = NaN;
end
end

function label = sig_label_local(P)
if P.diffFdr
    label = sprintf('dots: FDR q<%.2f', P.diffAlpha);
else
    label = sprintf('dots: p<%.2f', P.diffAlpha);
end
end

function p = one_sample_p_local(x, alpha, prefTest)
x = x(:);
x = x(isfinite(x));
if numel(x) < 3
    p = NaN;
    return;
end

pref = lower(string(prefTest));
if pref == "ttest"
    [~, p] = ttest(x, 0, 'Alpha', alpha);
    return;
end

if exist('signrank', 'file') == 2
    p = signrank(x, 0, 'alpha', alpha);
else
    [~, p] = ttest(x, 0, 'Alpha', alpha);
end
end

function [qVals, sig] = fdr_bh_local(pVals, alpha)
pVals = double(pVals(:));
qVals = nan(size(pVals));
sig = false(size(pVals));
good = isfinite(pVals);
if ~any(good)
    return;
end
p = pVals(good);
[ps, ord] = sort(p);
m = numel(ps);
thr = alpha * (1:m)' / m;
pass = ps <= thr;
if any(pass)
    k = find(pass, 1, 'last');
    sigSorted = false(m,1);
    sigSorted(1:k) = true;
else
    sigSorted = false(m,1);
end

qSorted = nan(m,1);
qTmp = (m ./ (1:m)') .* ps;
qTmp = flipud(cummin(flipud(qTmp)));
qSorted(:) = min(qTmp, 1);

idx = find(good);
sigIdx = idx(ord(sigSorted));
sig(sigIdx) = true;
qVals(idx(ord)) = qSorted;
qVals = reshape(qVals, size(pVals));
sig = reshape(sig, size(pVals));
end

function d = local_angdiff_deg(a, b)
d = mod((double(a) - double(b)) + 180, 360) - 180;
end
