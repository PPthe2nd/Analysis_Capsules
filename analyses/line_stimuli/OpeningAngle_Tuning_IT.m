function OUT = OpeningAngle_Tuning_IT()
% OPENINGANGLE_TUNING_IT
% Fit signed IT tuning to the opening angle of each quartet.
%
% The response for each site and quartet is the signed trial-weighted mean
% response across the 4 stimuli in that quartet, normalized as:
%   (muQuartet - muSpont) / sdSpont
%
% A weighted quadratic model is fit separately in the early and late windows:
%   y = b + a1*(ang-ang0) + a2*(ang-ang0)^2
%
% Output metrics include:
% - weighted in-sample variance explained (%)
% - approximate weighted F-test p-value versus a constant model
% - fitted effect size (observed fitted modulation range) in spont. SD units
% - preferred observed opening angle

%% Settings
P = struct();
P.Monkey = 1;           % 1 = Nilson, 2 = Figaro
P.minQuartets = 20;     % minimum number of finite quartet responses required per site
P.makeSummaryFigures = true;
P.makeExampleFigures = true;
P.nExampleSites = 4;    % per window
P.saveResult = true;
P.forceRefit = false;
P.progressEvery = 10;
P.vePlotFloorPct = -100;
P.sigAlpha = 0.05;      % per-site significance threshold on the approximate p-value
P.varFloorFrac = 1e-3;  % floor on quartet-response variance relative to the median

cfg = config();

%% Monkey-specific files
if P.Monkey == 1
    monkeySuffix = "N";
    tallFile = 'Tall_IT_lines_N.mat';
    resp3binFile = 'SNR_capsules_N_d12.mat';
elseif P.Monkey == 2
    monkeySuffix = "F";
    tallFile = 'Tall_IT_lines_F.mat';
    resp3binFile = 'SNR_capsules_F_d12.mat';
else
    error('OpeningAngle_Tuning_IT:InvalidMonkey', ...
        'P.Monkey must be 1 (Nilson) or 2 (Figaro).');
end

tallPath = fullfile(cfg.matDir, tallFile);
resp3binPath = fullfile(cfg.matDir, resp3binFile);
outPath = fullfile(cfg.matDir, sprintf('OpeningAngle_Tuning_IT_%s.mat', char(monkeySuffix)));

assert(exist(tallPath, 'file') == 2, ...
    'Missing %s. Run Line_Stimuli_IT.m first.', tallPath);
assert(exist(resp3binPath, 'file') == 2, ...
    'Missing %s. Create the 3-bin response summary first.', resp3binPath);

useCached = exist(outPath, 'file') == 2 && ~P.forceRefit;
needsSave = ~useCached;

if useCached
    fprintf('Loading cached IT opening-angle tuning from %s\n', outPath);
    S = load(outPath);
    assert(isfield(S, 'OUT') && isstruct(S.OUT), ...
        '%s must contain struct OUT.', outPath);
    OUT = S.OUT;
    assert(isfield(OUT, 'QuartetTable') && isfield(OUT, 'FitEarly') && ...
        isfield(OUT, 'FitLate') && isfield(OUT, 'signedQuartetEarly') && ...
        isfield(OUT, 'signedQuartetLate') && isfield(OUT, 'varQuartetEarly') && ...
        isfield(OUT, 'varQuartetLate') && isfield(OUT, 'RFrange'), ...
        'Cached OUT is missing required fields.');

    QuartetTable = OUT.QuartetTable;
    openingAngleDeg = double(QuartetTable.openingAngleDeg);
    FitEarly = OUT.FitEarly;
    FitLate = OUT.FitLate;
    signedQuartetEarly = OUT.signedQuartetEarly;
    signedQuartetLate = OUT.signedQuartetLate;
    varQuartetEarly = OUT.varQuartetEarly;
    varQuartetLate = OUT.varQuartetLate;
    RFrange = OUT.RFrange(:);
    nIT = numel(RFrange);
else
    %% Load geometry
    Sgeo = load(tallPath);
    assert(isfield(Sgeo, 'Tall_IT') && isstruct(Sgeo.Tall_IT), ...
        '%s must contain struct Tall_IT.', tallFile);
    assert(isfield(Sgeo, 'ALLCOORDS') && isfield(Sgeo, 'RFrange'), ...
        '%s must contain ALLCOORDS and RFrange.', tallFile);

    Tall_IT = Sgeo.Tall_IT;
    ALLCOORDS = Sgeo.ALLCOORDS;
    RFrange = Sgeo.RFrange(:);
    nIT = numel(RFrange);
    siteRows = (1:nIT).';

    stimNums = arrayfun(@(x) x.stimNum, Tall_IT(:));
    [stimNumsSorted, ordStim] = sort(stimNums(:));
    assert(all(stimNumsSorted(:).' == 1:numel(Tall_IT)), ...
        'Tall_IT.stimNum must cover 1..%d exactly.', numel(Tall_IT));
    Tall_IT = Tall_IT(ordStim);
    nStim = numel(Tall_IT);

    [quartetMembers, openingAngleDeg, stimRef, cueXY] = build_quartet_members_and_angle(ALLCOORDS, nStim);
    nQuartets = size(quartetMembers, 1);

    QuartetTable = table((1:nQuartets).', stimRef(:), cueXY(:,1), cueXY(:,2), openingAngleDeg(:), ...
        'VariableNames', {'quartetIdx','stimRef','cueX','cueY','openingAngleDeg'});

    %% Load 3-bin responses and localize to IT rows
    Sresp = load(resp3binPath);
    assert(isfield(Sresp, 'R') && isstruct(Sresp.R), ...
        '%s must contain struct R.', resp3binFile);

    R3_full = Sresp.R;
    R3 = R3_full;
    R3.meanAct = R3_full.meanAct(RFrange, :, :);
    R3.meanSqAct = R3_full.meanSqAct(RFrange, :, :);
    if ismatrix(R3_full.nTrials) && size(R3_full.nTrials,1) >= max(RFrange)
        R3.nTrials = R3_full.nTrials(RFrange, :);
    else
        R3.nTrials = R3_full.nTrials;
    end

    assert(size(R3.meanAct,1) == nIT, 'Localized IT response rows do not match Tall_IT.');
    assert(size(R3.meanAct,2) == nStim, 'Localized IT response stimuli do not match Tall_IT.');

    SNR = compute_snr_per_color_sites(R3, Tall_IT, siteRows, 'Verbose', false);
    muSpont = SNR.muSpont(siteRows);
    sdSpont = SNR.sdSpont(siteRows);

    [signedQuartetEarly, signedQuartetLate, varQuartetEarly, varQuartetLate] = ...
        compute_signed_quartet_responses_with_var(R3, quartetMembers, muSpont, sdSpont);

    %% Fit opening-angle model
    fields = {'status','nObs','baseline','coefLin','coefQuad','angleCenter', ...
        'preferredObservedAngleDeg','rmse','r2TrainPct','effectRangeObs', ...
        'effectAbsPeakObs','fStatApprox','pValueApprox','dfModel','dfError', ...
        'rssNullWeighted','rssFullWeighted'};
    FitEarly = init_fit_struct(nIT, fields);
    FitLate = init_fit_struct(nIT, fields);

    fprintf('Fitting IT opening-angle tuning for monkey %s (%d sites, %d quartets)\n', ...
        char(monkeySuffix), nIT, nQuartets);
    tFit = tic;

    for iSite = 1:nIT
        yEarly = signedQuartetEarly(iSite, :).';
        yLate = signedQuartetLate(iSite, :).';
        vEarly = varQuartetEarly(iSite, :).';
        vLate = varQuartetLate(iSite, :).';

        FitEarly(iSite) = fit_weighted_quadratic(openingAngleDeg, yEarly, vEarly, P);
        FitLate(iSite) = fit_weighted_quadratic(openingAngleDeg, yLate, vLate, P);

        if P.progressEvery > 0 && (iSite == 1 || mod(iSite, P.progressEvery) == 0 || iSite == nIT)
            elapsedSec = toc(tFit);
            rate = iSite / max(elapsedSec, eps);
            etaSec = (nIT - iSite) / max(rate, eps);
            fprintf('  Site %d / %d | elapsed %.1fs | ETA %.1fs\n', ...
                iSite, nIT, elapsedSec, etaSec);
        end
    end
end

veEarly = [FitEarly.r2TrainPct].';
veLate = [FitLate.r2TrainPct].';
effEarly = [FitEarly.effectRangeObs].';
effLate = [FitLate.effectRangeObs].';
prefEarly = [FitEarly.preferredObservedAngleDeg].';
prefLate = [FitLate.preferredObservedAngleDeg].';
pEarly = [FitEarly.pValueApprox].';
pLate = [FitLate.pValueApprox].';

isAngleTunedEarly = isfinite(pEarly) & (pEarly < P.sigAlpha);
isAngleTunedLate = isfinite(pLate) & (pLate < P.sigAlpha);

fprintf('Usable IT sites for opening-angle fitting (early): %d / %d\n', ...
    nnz(isfinite(veEarly)), nIT);
fprintf('Usable IT sites for opening-angle fitting (late):  %d / %d\n', ...
    nnz(isfinite(veLate)), nIT);
fprintf('Opening-angle tuned sites at p < %.3f | early=%d late=%d\n', ...
    P.sigAlpha, nnz(isAngleTunedEarly), nnz(isAngleTunedLate));

%% Summary figures
if P.makeSummaryFigures
    figure('Color', 'w');
    useTiled = exist('tiledlayout', 'file') == 2;
    if useTiled
        tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');
    end

    if useTiled, nexttile; else, subplot(1,2,1); end
    plot_ve_histogram(veEarly, veLate, P.vePlotFloorPct, ...
        sprintf('Opening-angle variance explained (%s)', char(monkeySuffix)));

    if useTiled, nexttile; else, subplot(1,2,2); end
    plot_effect_histogram(effEarly, effLate, ...
        sprintf('Opening-angle effect size (%s)', char(monkeySuffix)));

    figure('Color', 'w');
    histogram(prefEarly(isfinite(prefEarly) & isAngleTunedEarly), 18, ...
        'FaceColor', [0.25 0.45 0.85], 'EdgeColor', 'none');
    hold on;
    histogram(prefLate(isfinite(prefLate) & isAngleTunedLate), 18, ...
        'FaceColor', [0.85 0.35 0.25], 'EdgeColor', 'none');
    xlabel('Preferred observed opening angle (deg)');
    ylabel('N sites');
    title(sprintf('Opening-angle preference among significant sites (%s)', char(monkeySuffix)));
    legend('Early','Late');
    grid on;
end

%% Example figures
if P.makeExampleFigures
    siteSelEarly = select_example_sites(isAngleTunedEarly, veEarly, P.nExampleSites);
    siteSelLate = select_example_sites(isAngleTunedLate, veLate, P.nExampleSites);

    make_example_angle_fit_figure(openingAngleDeg, signedQuartetEarly, FitEarly, ...
        sprintf('IT opening-angle fits early (%s)', char(monkeySuffix)), siteSelEarly);
    make_example_angle_fit_figure(openingAngleDeg, signedQuartetLate, FitLate, ...
        sprintf('IT opening-angle fits late (%s)', char(monkeySuffix)), siteSelLate);
end

%% Tables
Tearly = build_fit_table(FitEarly, RFrange, 'r2TrainPct');
Tlate = build_fit_table(FitLate, RFrange, 'r2TrainPct');

disp('Top IT opening-angle sites by variance explained (early):');
disp(Tearly(1:min(10,height(Tearly)), :));
disp('Top IT opening-angle sites by variance explained (late):');
disp(Tlate(1:min(10,height(Tlate)), :));

%% Pack output
OUT.P = P;
OUT.monkeySuffix = monkeySuffix;
OUT.RFrange = RFrange;
OUT.QuartetTable = QuartetTable;
OUT.signedQuartetEarly = signedQuartetEarly;
OUT.signedQuartetLate = signedQuartetLate;
OUT.varQuartetEarly = varQuartetEarly;
OUT.varQuartetLate = varQuartetLate;
OUT.FitEarly = FitEarly;
OUT.FitLate = FitLate;
OUT.isAngleTunedEarly = isAngleTunedEarly;
OUT.isAngleTunedLate = isAngleTunedLate;
OUT.TableEarly = Tearly;
OUT.TableLate = Tlate;

if needsSave && P.saveResult
    save(outPath, 'OUT', '-v7.3');
    fprintf('Saved IT opening-angle tuning to %s\n', outPath);
end
end

function [quartetMembers, openingAngleDeg, stimRef, cueXY] = build_quartet_members_and_angle(ALLCOORDS, nStim)
nQuartets = floor(nStim / 8) * 2;
quartetMembers = zeros(nQuartets, 4);
openingAngleDeg = nan(nQuartets, 1);
stimRef = nan(nQuartets, 1);
cueXY = nan(nQuartets, 2);
q = 0;

for base = 0:8:(nStim - 8)
    q = q + 1;
    quartetMembers(q,:) = base + [1 2 5 6];
    [openingAngleDeg(q), stimRef(q), cueXY(q,:)] = angle_from_quartet(ALLCOORDS, quartetMembers(q,:));

    q = q + 1;
    quartetMembers(q,:) = base + [3 4 7 8];
    [openingAngleDeg(q), stimRef(q), cueXY(q,:)] = angle_from_quartet(ALLCOORDS, quartetMembers(q,:));
end
end

function [angleDeg, stimRef, cue] = angle_from_quartet(ALLCOORDS, stimQ)
stimRef = stimQ(1);
[angle0, cue0] = opening_angle_from_stim(ALLCOORDS, stimRef);

for k = 2:numel(stimQ)
    [angleK, cueK] = opening_angle_from_stim(ALLCOORDS, stimQ(k));
    assert(norm(cueK - cue0) < 1e-6, ...
        'Cue position differs within quartet [%s].', num2str(stimQ));
    assert(abs(angleK - angle0) < 1e-6, ...
        'Opening angle differs within quartet [%s].', num2str(stimQ));
end

angleDeg = angle0;
cue = cue0;
end

function [angleDeg, cue] = opening_angle_from_stim(ALLCOORDS, stimNum)
f = sprintf('stim_%d', stimNum);
s = double(ALLCOORDS.(f).s(:))';
tFig = double(ALLCOORDS.(f).t_fig(:))';
tBack = double(ALLCOORDS.(f).t_back(:))';

v1 = tFig - s;
v2 = tBack - s;
c = dot(v1, v2) / (norm(v1) * norm(v2));
c = min(max(c, -1), 1);
angleDeg = acosd(c);
cue = s;
end

function [signedEarly, signedLate, varEarly, varLate] = compute_signed_quartet_responses_with_var(R3, quartetMembers, muSpont, sdSpont)
[nSites, ~, ~] = size(R3.meanAct);
nQuartets = size(quartetMembers, 1);

if isvector(R3.nTrials)
    nTrialsAll = double(R3.nTrials(:));
    perSiteTrials = false;
else
    perSiteTrials = true;
end

muQuartetEarly = nan(nSites, nQuartets);
muQuartetLate = nan(nSites, nQuartets);
varQuartetEarly = nan(nSites, nQuartets);
varQuartetLate = nan(nSites, nQuartets);

for iSite = 1:nSites
    if perSiteTrials
        nTrSite = double(R3.nTrials(iSite, :)).';
    else
        nTrSite = nTrialsAll;
    end

    rEarly = squeeze(R3.meanAct(iSite,:,2)).';
    rLate = squeeze(R3.meanAct(iSite,:,3)).';
    rSqEarly = squeeze(R3.meanSqAct(iSite,:,2)).';
    rSqLate = squeeze(R3.meanSqAct(iSite,:,3)).';

    for qIdx = 1:nQuartets
        stimQ = quartetMembers(qIdx,:);
        nTrQ = nTrSite(stimQ);

        rEarlyQ = rEarly(stimQ);
        rSqEarlyQ = rSqEarly(stimQ);
        idxEarly = isfinite(rEarlyQ) & isfinite(rSqEarlyQ) & isfinite(nTrQ) & (nTrQ > 1);
        if any(idxEarly)
            nUse = nTrQ(idxEarly);
            muUse = rEarlyQ(idxEarly);
            msqUse = rSqEarlyQ(idxEarly);
            varStim = max(0, msqUse - muUse.^2);
            varStim = varStim .* (nUse ./ max(nUse - 1, 1));
            muQuartetEarly(iSite, qIdx) = sum(nUse .* muUse) / sum(nUse);
            varQuartetEarly(iSite, qIdx) = sum(nUse .* varStim) / (sum(nUse)^2);
        end

        rLateQ = rLate(stimQ);
        rSqLateQ = rSqLate(stimQ);
        idxLate = isfinite(rLateQ) & isfinite(rSqLateQ) & isfinite(nTrQ) & (nTrQ > 1);
        if any(idxLate)
            nUse = nTrQ(idxLate);
            muUse = rLateQ(idxLate);
            msqUse = rSqLateQ(idxLate);
            varStim = max(0, msqUse - muUse.^2);
            varStim = varStim .* (nUse ./ max(nUse - 1, 1));
            muQuartetLate(iSite, qIdx) = sum(nUse .* muUse) / sum(nUse);
            varQuartetLate(iSite, qIdx) = sum(nUse .* varStim) / (sum(nUse)^2);
        end
    end
end

signedEarly = bsxfun(@rdivide, bsxfun(@minus, muQuartetEarly, muSpont), sdSpont);
signedLate = bsxfun(@rdivide, bsxfun(@minus, muQuartetLate, muSpont), sdSpont);
varEarly = bsxfun(@rdivide, varQuartetEarly, sdSpont.^2);
varLate = bsxfun(@rdivide, varQuartetLate, sdSpont.^2);

badNoise = ~isfinite(sdSpont) | (sdSpont <= 0);
signedEarly(badNoise, :) = NaN;
signedLate(badNoise, :) = NaN;
varEarly(badNoise, :) = NaN;
varLate(badNoise, :) = NaN;
end

function F = init_fit_struct(nSites, fields)
tmp = struct();
for i = 1:numel(fields)
    tmp.(fields{i}) = NaN;
end
tmp.status = "not_fit";
F = repmat(tmp, nSites, 1);
end

function F = fit_weighted_quadratic(angleDeg, y, varY, P)
fields = {'status','nObs','baseline','coefLin','coefQuad','angleCenter', ...
    'preferredObservedAngleDeg','rmse','r2TrainPct','effectRangeObs', ...
    'effectAbsPeakObs','fStatApprox','pValueApprox','dfModel','dfError', ...
    'rssNullWeighted','rssFullWeighted'};
F = struct();
for i = 1:numel(fields)
    F.(fields{i}) = NaN;
end
F.status = "not_fit";

valid = isfinite(angleDeg) & isfinite(y) & isfinite(varY) & (varY > 0);
ang = angleDeg(valid);
yv = y(valid);
vv = varY(valid);
F.nObs = numel(yv);

if numel(yv) < P.minQuartets
    F.status = "too_few_quartets";
    return;
end

varFloor = max(median(vv(isfinite(vv) & vv > 0)) * P.varFloorFrac, 1e-6);
vv = max(vv, varFloor);
sqrtW = sqrt(1 ./ vv);

ang0 = mean(ang);
x = ang - ang0;
X = [ones(numel(x),1), x, x.^2];
Xw = bsxfun(@times, X, sqrtW);
yw = sqrtW .* yv;

beta = Xw \ yw;
yHat = X * beta;
if any(~isfinite(beta))
    F.status = "fit_failed";
    return;
end

wMean = sum((1 ./ vv) .* yv) / sum(1 ./ vv);
rssNull = sum((1 ./ vv) .* (yv - wMean).^2);
rssFull = sum((1 ./ vv) .* (yv - yHat).^2);
dfModel = 2;
dfError = numel(yv) - size(X,2);

F.status = "ok";
F.baseline = beta(1);
F.coefLin = beta(2);
F.coefQuad = beta(3);
F.angleCenter = ang0;
F.rmse = sqrt(mean((yv - yHat).^2));
if rssNull > 0
    F.r2TrainPct = 100 * (1 - rssFull / rssNull);
end
F.effectRangeObs = max(yHat) - min(yHat);
F.effectAbsPeakObs = max(abs(yHat));
F.rssNullWeighted = rssNull;
F.rssFullWeighted = rssFull;
F.dfModel = dfModel;
F.dfError = dfError;

[~, iMax] = max(yHat);
F.preferredObservedAngleDeg = ang(iMax);

if dfError > 0 && rssNull > rssFull
    F.fStatApprox = ((rssNull - rssFull) / dfModel) / (rssFull / dfError);
    F.pValueApprox = 1 - fcdf_local(F.fStatApprox, dfModel, dfError);
else
    F.fStatApprox = 0;
    F.pValueApprox = 1;
end
end

function idx = select_example_sites(isSig, ve, nSelect)
idxAll = find(isSig & isfinite(ve));
if isempty(idxAll)
    idxAll = find(isfinite(ve));
end
if isempty(idxAll)
    idx = [];
    return;
end
[~, ord] = sort(ve(idxAll), 'descend');
idx = idxAll(ord(1:min(nSelect, numel(ord))));
end

function plot_ve_histogram(veEarly, veLate, floorPct, ttl)
vePlotEarly = veEarly;
vePlotLate = veLate;
vePlotEarly(isfinite(vePlotEarly) & (vePlotEarly < floorPct)) = floorPct;
vePlotLate(isfinite(vePlotLate) & (vePlotLate < floorPct)) = floorPct;

histogram(vePlotEarly(isfinite(vePlotEarly)), 30, 'FaceColor', [0.25 0.45 0.85], 'EdgeColor', 'none');
hold on;
histogram(vePlotLate(isfinite(vePlotLate)), 30, 'FaceColor', [0.85 0.35 0.25], 'EdgeColor', 'none');
xlabel('Variance explained (%)');
ylabel('N sites');
title(ttl);
legend('Early','Late');
grid on;
end

function plot_effect_histogram(effEarly, effLate, ttl)
histogram(effEarly(isfinite(effEarly)), 30, 'FaceColor', [0.25 0.45 0.85], 'EdgeColor', 'none');
hold on;
histogram(effLate(isfinite(effLate)), 30, 'FaceColor', [0.85 0.35 0.25], 'EdgeColor', 'none');
xlabel('Fitted modulation range (spont. SD units)');
ylabel('N sites');
title(ttl);
legend('Early','Late');
grid on;
end

function make_example_angle_fit_figure(angleDeg, Ymat, FitStruct, figName, siteSel)
if isempty(siteSel)
    return;
end

nShow = numel(siteSel);
fig = figure('Name', figName, 'NumberTitle', 'off', 'Color', 'w');
useTiled = exist('tiledlayout', 'file') == 2;
if useTiled
    tl = tiledlayout(fig, 1, nShow, 'TileSpacing', 'compact', 'Padding', 'compact');
end
angleGrid = linspace(min(angleDeg), max(angleDeg), 200).';

for i = 1:nShow
    s = siteSel(i);
    if useTiled
        ax = nexttile(tl);
    else
        ax = subplot(1, nShow, i);
    end

    y = Ymat(s, :).';
    valid = isfinite(y) & isfinite(angleDeg);
    ang = angleDeg(valid);
    yv = y(valid);
    [angSort, ord] = sort(ang);
    ySort = yv(ord);

    scatter(ax, angSort, ySort, 45, 'filled');
    hold(ax, 'on');
    fitY = quadratic_curve(angleGrid, FitStruct(s));
    plot(ax, angleGrid, fitY, 'k-', 'LineWidth', 1.5);
    xlabel(ax, 'Opening angle (deg)');
    ylabel(ax, 'Signed response');
    title(ax, sprintf('site %d\nVE %.1f%% | eff %.2f | p %.3g', s, ...
        FitStruct(s).r2TrainPct, FitStruct(s).effectRangeObs, FitStruct(s).pValueApprox), ...
        'FontSize', 10);
    grid(ax, 'on');
end
if exist('sgtitle', 'file') == 2
    sgtitle(figName);
end
end

function yHat = quadratic_curve(angleDeg, F)
x = angleDeg - F.angleCenter;
yHat = F.baseline + F.coefLin .* x + F.coefQuad .* (x.^2);
end

function T = build_fit_table(FitStruct, RFrange, sortField)
n = numel(FitStruct);
statusText = string({FitStruct.status}).';
T = table((1:n).', RFrange(:), [FitStruct.nObs].', ...
    [FitStruct.preferredObservedAngleDeg].', [FitStruct.effectRangeObs].', ...
    [FitStruct.r2TrainPct].', [FitStruct.pValueApprox].', statusText, ...
    'VariableNames', {'localSite','globalSiteInR','nObs','preferredObservedAngleDeg', ...
    'effectRangeObs','r2TrainPct','pValueApprox','status'});

ok = isfinite(T.(sortField));
T = sortrows(T(ok, :), sortField, 'descend');
end

function p = fcdf_local(x, df1, df2)
z = (df1 .* x) ./ (df1 .* x + df2);
z = min(max(z, 0), 1);
p = betainc(z, df1/2, df2/2);
end
