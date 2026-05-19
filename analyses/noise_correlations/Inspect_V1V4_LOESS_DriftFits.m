function OUT = Inspect_V1V4_LOESS_DriftFits(Monkey, Puser)
% INSPECT_V1V4_LOESS_DRIFTFITS
% Visualize raw prestim baseline traces and LOESS fits for example V1/V4
% sites from the V1-V4 drift-corrected noise-correlation analysis.

if nargin < 1 || isempty(Monkey)
    Monkey = 1;
end
if nargin < 2 || isempty(Puser)
    Puser = struct();
end

P = struct();
P.cacheTag = "all";
P.siteGlobalV1 = [];
P.siteGlobalV4 = [];
P.spans = [0.2 0.3 0.4];
P.preWindowMs = [-200 0];
P.sessions = [1 2];
P.onlyCorrect = true;
P.correctCol = 9;
P.correctVal = 1;
P.stimCol = 1;
P.dayCol = 11;
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
loessPath = fullfile(cfg.resultsDir, sprintf('NoiseCorr_V1V4_RFDistance_TargetDistractor_LOESS_%s_%s.mat', char(monkeySuffix), char(P.cacheTag)));
assert(exist(loessPath, 'file') == 2, 'Missing %s. Run the LOESS analysis first.', loessPath);
S = load(loessPath, 'OUT');
assert(isfield(S, 'OUT'), '%s must contain OUT.', loessPath);
LO = S.OUT;

if isempty(P.siteGlobalV1)
    P.siteGlobalV1 = double(LO.loess.exampleV1.globalSite);
end
if isempty(P.siteGlobalV4)
    P.siteGlobalV4 = double(LO.loess.exampleV4.globalSite);
end

trialRespPath = fullfile(cfg.dataRoot, monkeyFolder, 'ObjAtt_lines_normMUA.mat');
trialInfoPath = fullfile(cfg.dataRoot, monkeyFolder, 'ObjAtt_lines_MUA_trials.mat');
assert(exist(trialRespPath, 'file') == 2, 'Missing %s.', trialRespPath);
assert(exist(trialInfoPath, 'file') == 2, 'Missing %s.', trialInfoPath);
m1 = matfile(trialRespPath);
m2 = matfile(trialInfoPath);
ALLMAT = double(m2.ALLMAT);
tb = double(m2.tb);
tb = tb(:);

stimAll = double(ALLMAT(:, P.stimCol));
dayAll = double(ALLMAT(:, P.dayCol));
trialInclude = isfinite(stimAll) & stimAll >= 1 & stimAll <= 384 & ...
    (floor(stimAll) == stimAll);
if P.onlyCorrect
    trialInclude = trialInclude & (double(ALLMAT(:, P.correctCol)) == P.correctVal);
end
trialInclude = trialInclude & ismember(dayAll, P.sessions(:));
trialIdx = find(trialInclude);
dayIncl = dayAll(trialIdx);

preIdx = find(tb >= P.preWindowMs(1) & tb < P.preWindowMs(2));
assert(~isempty(preIdx), 'No prestim samples in [%d %d).', P.preWindowMs(1), P.preWindowMs(2));

siteList = [P.siteGlobalV1, P.siteGlobalV4];
siteLabels = {sprintf('V1 global %d', P.siteGlobalV1), sprintf('V4 global %d', P.siteGlobalV4)};
siteRaw = nan(2, numel(trialIdx));
for i = 1:2
    siteRaw(i, :) = local_load_prestim_site(m1, siteList(i), preIdx, trialIdx);
end

fits = nan(2, numel(trialIdx), numel(P.spans));
for i = 1:2
    for s = 1:numel(P.spans)
        fits(i, :, s) = local_fit_loess(siteRaw(i, :), dayIncl, P.sessions, P.spans(s));
    end
end

OUT = struct();
OUT.P = P;
OUT.monkeySuffix = monkeySuffix;
OUT.dayIncl = dayIncl;
OUT.siteList = siteList;
OUT.siteLabels = siteLabels;
OUT.siteRaw = siteRaw;
OUT.fits = fits;

if P.verbose
    fprintf('Inspecting LOESS drift fits for %s and %s\n', siteLabels{1}, siteLabels{2});
end

if P.plotFigure
    local_plot_summary(OUT);
end
end

function x = local_load_prestim_site(m1, siteGlobal, preIdx, trialIdx)
X = double(m1.normMUA(siteGlobal, :, preIdx));
X = squeeze(mean(X, 3, 'omitnan'));
x = reshape(X(trialIdx), 1, []);
end

function fitVals = local_fit_loess(x, dayIncl, sessions, frac)
fitVals = nan(size(x));
for k = 1:numel(sessions)
    idx = find(dayIncl == sessions(k));
    if numel(idx) < 5
        continue;
    end
    y = reshape(x(idx), [], 1);
    good = isfinite(y);
    if nnz(good) < 5
        continue;
    end
    spanN = min(nnz(good), max(5, ceil(frac * nnz(good))));
    if spanN >= nnz(good)
        yfit = repmat(mean(y(good), 'omitnan'), nnz(good), 1);
    else
        yfit = smoothdata(y(good), 'loess', spanN);
    end
    tmp = nan(size(y));
    tmp(good) = yfit;
    fitVals(idx) = tmp;
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

function local_plot_summary(OUT)
colors = [0.20 0.60 0.80; 0.85 0.45 0.20; 0.30 0.65 0.30];

figure('Color', 'w', 'Name', sprintf('LOESS drift fits (%s)', char(OUT.monkeySuffix)));
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

for i = 1:2
    for k = 1:2
        ax = nexttile;
        hold(ax, 'on');
        idx = find(OUT.dayIncl == OUT.P.sessions(k));
        x = 1:numel(idx);
        y = reshape(OUT.siteRaw(i, idx), [], 1);
        plot(ax, x, y, '.', 'Color', [0.75 0.75 0.75], 'MarkerSize', 6, 'DisplayName', 'raw');
        for s = 1:numel(OUT.P.spans)
            yfit = reshape(OUT.fits(i, idx, s), [], 1);
            plot(ax, x, yfit, '-', 'Color', colors(s, :), 'LineWidth', 1.8, ...
                'DisplayName', sprintf('loess %.1f', OUT.P.spans(s)));
        end
        xlabel(ax, sprintf('Session %d trial order', OUT.P.sessions(k)));
        ylabel(ax, 'Prestim baseline');
        title(ax, sprintf('%s | session %d', OUT.siteLabels{i}, OUT.P.sessions(k)));
        grid(ax, 'on');
        if i == 1 && k == 1
            legend(ax, 'Location', 'best', 'Box', 'off');
        end
    end
end
end
