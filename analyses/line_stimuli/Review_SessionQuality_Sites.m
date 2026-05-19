function OUT = Review_SessionQuality_Sites(Puser)
% REVIEW_SESSIONQUALITY_SITES
% Interactive reviewer for site-by-session quality issues flagged by large
% day-1/day-2 changes in late-window trial variability.
%
% For each candidate site, the script shows:
%   1) all-stimulus session mean timecourses (session 1 vs session 2)
%   2) single-trial prestim means across trial order, colored by session
%
% Review controls:
%   Keep      -> keep both sessions
%   Excl sugg -> exclude the suggested bad session
%   Excl s1   -> exclude session 1
%   Excl s2   -> exclude session 2
%   Quit      -> quit review

P = struct();
P.Monkey = 1;                    % 1 = Nilson, 2 = Figaro
P.Areas = {'IT','V1','V4'};
P.UseResponsiveOnly = false;
P.LateSdRatioLow = 0.5;
P.LateSdRatioHigh = 2.0;
P.onlyCorrect = true;
P.correctCol = 9;
P.correctVal = 1;
P.dayCol = 11;
P.sessions = [1 2];
P.expectedMaxStim = 384;
P.chunkTrials = 200;
P.preWindowMs = [-200 0];
P.StartIdx = 1;
P.MaxCandidates = inf;
P.Interactive = true;
P.PlotFigures = true;
P.UseDialog = true;

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
elseif P.Monkey == 2
    monkeySuffix = "F";
    monkeyFolder = 'Figaro';
else
    error('Review_SessionQuality_Sites:InvalidMonkey', ...
        'Monkey must be 1 (Nilson) or 2 (Figaro).');
end

capsPath = fullfile(cfg.matDir, sprintf('SNR_capsules_%s_d12.mat', char(monkeySuffix)));
trialRespPath = fullfile(cfg.dataRoot, monkeyFolder, 'ObjAtt_lines_normMUA.mat');
trialInfoPath = fullfile(cfg.dataRoot, monkeyFolder, 'ObjAtt_lines_MUA_trials.mat');
assert(exist(capsPath, 'file') == 2, 'Missing %s.', capsPath);
assert(exist(trialRespPath, 'file') == 2, 'Missing %s.', trialRespPath);
assert(exist(trialInfoPath, 'file') == 2, 'Missing %s.', trialInfoPath);

[Split, ~, ~] = load_capsules_day_split(capsPath, monkeySuffix, ...
    'cfg', cfg, ...
    'sessions', P.sessions, ...
    'onlyCorrect', P.onlyCorrect, ...
    'correctCol', P.correctCol, ...
    'correctVal', P.correctVal, ...
    'dayCol', P.dayCol, ...
    'chunkTrials', P.chunkTrials, ...
    'expectedMaxStim', P.expectedMaxStim, ...
    'useCache', true);

areaList = string(P.Areas(:));
Cand = table();
for iA = 1:numel(areaList)
    Ci = build_area_candidates_local(char(areaList(iA)), monkeySuffix, Split, P, cfg);
    Cand = [Cand; Ci]; %#ok<AGROW>
end

if isempty(Cand)
    fprintf('No session-quality candidates found for monkey %s.\n', char(monkeySuffix));
    OUT = struct('P', P, 'monkeySuffix', monkeySuffix, 'Candidates', Cand, 'Decisions', table());
    return;
end

Cand = sortrows(Cand, {'severity','areaName','siteGlobal'}, {'descend','ascend','ascend'});
if isfinite(P.MaxCandidates)
    lastIdx = min(height(Cand), P.StartIdx + P.MaxCandidates - 1);
else
    lastIdx = height(Cand);
end
Cand = Cand(P.StartIdx:lastIdx, :);

fprintf('Reviewing %d session-quality candidates for monkey %s\n', height(Cand), char(monkeySuffix));
disp(Cand(:, {'areaName','siteGlobal','lateSdRatio','preSdRatio','lateEvSess1','lateEvSess2','suggestedExcludeDay'}));

m1 = matfile(trialRespPath);
m2 = matfile(trialInfoPath);
ALLMAT = m2.ALLMAT;
tb = double(m2.tb);
tb = tb(:);
preIdx = find(tb >= P.preWindowMs(1) & tb < P.preWindowMs(2));
assert(~isempty(preIdx), 'No prestim raw samples found in [%d %d).', P.preWindowMs(1), P.preWindowMs(2));

stimPerTrial = double(ALLMAT(:,1)); %#ok<NASGU>
sessionPerTrial = double(ALLMAT(:,P.dayCol));
trialInclude = true(size(sessionPerTrial));
if P.onlyCorrect
    trialInclude = trialInclude & (double(ALLMAT(:,P.correctCol)) == P.correctVal);
end
trialInclude = trialInclude & ismember(sessionPerTrial, P.sessions(:));

Decision = table( ...
    strings(0,1), zeros(0,1), zeros(0,1), strings(0,1), strings(0,1), ...
    'VariableNames', {'areaName','siteGlobal','excludedDay','action','note'});

for iC = 1:height(Cand)
    row = Cand(iC,:);
    diag = build_site_diagnostic_local(row, m1, tb, preIdx, trialInclude, sessionPerTrial, P);
    if P.PlotFigures
        fig = make_review_figure_local(row, diag, P);
        drawnow;
    else
        fig = [];
    end

    if P.Interactive
        if P.UseDialog
            resp = ask_review_decision_local(row, iC, height(Cand), fig);
        else
            prompt = sprintf(['[%d/%d] %s site %d | suggested exclude session %d ' ...
                '[y keep / n exclude suggested / 1 / 2 / q]: '], ...
                iC, height(Cand), char(row.areaName), row.siteGlobal, row.suggestedExcludeDay);
            resp = lower(strtrim(input(prompt, 's')));
        end
    else
        resp = 'y';
    end

    action = "keep";
    exclDay = 0;
    note = "";
    if strcmp(resp, 'q')
        if ~isempty(fig) && isvalid(fig)
            close(fig);
        end
        break;
    elseif strcmp(resp, 'n')
        action = "exclude";
        exclDay = row.suggestedExcludeDay;
    elseif strcmp(resp, '1') || strcmp(resp, '2')
        action = "exclude";
        exclDay = str2double(resp);
    else
        action = "keep";
    end

    if action == "exclude"
        note = sprintf('%s review: late SD ratio s2/s1=%.3f, late evoked [%.3f %.3f], pre means [%.3f %.3f]', ...
            char(row.areaName), row.lateSdRatio, row.lateEvSess1, row.lateEvSess2, row.preMeanSess1, row.preMeanSess2);
        upsert_site_session_exclusions_user(monkeySuffix, row.siteGlobal, exclDay, note);
        fprintf('Excluded %s site %d session %d\n', char(row.areaName), row.siteGlobal, exclDay);
    else
        fprintf('Kept %s site %d unchanged\n', char(row.areaName), row.siteGlobal);
    end

    Decision = [Decision; {row.areaName, row.siteGlobal, exclDay, action, note}]; %#ok<AGROW>

    if ~isempty(fig) && isvalid(fig)
        close(fig);
    end
end

OUT = struct();
OUT.P = P;
OUT.monkeySuffix = monkeySuffix;
OUT.Candidates = Cand;
OUT.Decisions = Decision;
end

function Cand = build_area_candidates_local(areaName, monkeySuffix, Split, P, cfg)
switch upper(areaName)
    case 'IT'
        resultPath = fullfile(cfg.resultsDir, sprintf('Compare_IT_SessionResponses_3bin_%s.mat', char(monkeySuffix)));
    case 'V1'
        resultPath = fullfile(cfg.resultsDir, sprintf('Compare_V1_SessionResponses_3bin_%s.mat', char(monkeySuffix)));
    case 'V4'
        resultPath = fullfile(cfg.resultsDir, sprintf('Compare_V4_SessionResponses_3bin_%s.mat', char(monkeySuffix)));
    otherwise
        error('Unknown area %s.', areaName);
end
assert(exist(resultPath, 'file') == 2, ...
    'Missing %s. Run the session-quality check for %s first.', resultPath, areaName);

S = load(resultPath, 'OUT');
O = S.OUT;
rows = double(O.siteGlobal(:));
keep = true(size(rows));
if P.UseResponsiveOnly && isfield(O, 'keepResponsive')
    keep = logical(O.keepResponsive(:));
end

nSite = numel(rows);
mu = nan(nSite, 3, 2);
sd = nan(nSite, 3, 2);
for iSess = 1:2
    n = double(Split(iSess).nTrials(:));
    for s = 1:nSite
        for w = 1:3
            x = squeeze(double(Split(iSess).meanAct(rows(s), :, w))).';
            y = squeeze(double(Split(iSess).meanSqAct(rows(s), :, w))).';
            idx = isfinite(x) & isfinite(y) & isfinite(n) & (n > 0);
            if any(idx)
                ww = n(idx);
                mu0 = sum(ww .* x(idx)) / sum(ww);
                ex2 = sum(ww .* y(idx)) / sum(ww);
                mu(s,w,iSess) = mu0;
                sd(s,w,iSess) = sqrt(max(0, ex2 - mu0.^2));
            end
        end
    end
end

preMean = squeeze(mu(:,1,:));
lateEv = squeeze(mu(:,3,:) - mu(:,1,:));
preSd = squeeze(sd(:,1,:));
lateSd = squeeze(sd(:,3,:));
preSdRatio = preSd(:,2) ./ preSd(:,1);
lateSdRatio = lateSd(:,2) ./ lateSd(:,1);
good = keep & isfinite(lateSdRatio) & ((lateSdRatio < P.LateSdRatioLow) | (lateSdRatio > P.LateSdRatioHigh));

if ~any(good)
    Cand = table();
    return;
end

suggested = repmat(P.sessions(1), nSite, 1);
isSess2Worse = lateEv(:,2) < lateEv(:,1);
suggested(isSess2Worse) = P.sessions(2);
tieMask = abs(lateEv(:,2) - lateEv(:,1)) < 1e-9;
sess2MoreExtremePre = abs(preMean(:,2)) > abs(preMean(:,1));
suggested(tieMask & sess2MoreExtremePre) = P.sessions(2);

severity = abs(log2(lateSdRatio));

ng = nnz(good);
areaCol = repmat(string(upper(areaName)), ng, 1);
siteCol = rows(good); siteCol = siteCol(:);
keepCol = keep(good); keepCol = keepCol(:);
pre1Col = preMean(good,1); pre1Col = pre1Col(:);
pre2Col = preMean(good,2); pre2Col = pre2Col(:);
preSdCol = preSdRatio(good); preSdCol = preSdCol(:);
lateSdCol = lateSdRatio(good); lateSdCol = lateSdCol(:);
late1Col = lateEv(good,1); late1Col = late1Col(:);
late2Col = lateEv(good,2); late2Col = late2Col(:);
suggestCol = suggested(good); suggestCol = suggestCol(:);
severityCol = severity(good); severityCol = severityCol(:);
Cand = table( ...
    areaCol, ...
    siteCol, ...
    keepCol, ...
    pre1Col, ...
    pre2Col, ...
    preSdCol, ...
    lateSdCol, ...
    late1Col, ...
    late2Col, ...
    suggestCol, ...
    severityCol, ...
    'VariableNames', { ...
        'areaName','siteGlobal','keepResponsive', ...
        'preMeanSess1','preMeanSess2','preSdRatio','lateSdRatio', ...
        'lateEvSess1','lateEvSess2','suggestedExcludeDay','severity'});
end

function z = robust_abs_z_local(x, ref)
med = median(ref, 'omitnan');
madv = median(abs(ref - med), 'omitnan');
scale = max(1.4826 * madv, 1e-6);
z = abs(x - med) ./ scale;
z(~isfinite(z)) = 0;
end

function z = robust_low_z_local(x, ref)
med = median(ref, 'omitnan');
madv = median(abs(ref - med), 'omitnan');
scale = max(1.4826 * madv, 1e-6);
z = max(0, (med - x) ./ scale);
z(~isfinite(z)) = 0;
end

function D = build_site_diagnostic_local(row, m1, tb, preIdx, trialInclude, sessionPerTrial, P)
siteGlobal = double(row.siteGlobal);
xPre = squeeze(mean(m1.normMUA(siteGlobal, :, preIdx), 3, 'omitnan'));
xPre = double(xPre(:));
tcAll = squeeze(double(m1.normMUA(siteGlobal, :, :)));
if isvector(tcAll)
    tcAll = reshape(tcAll, 1, []);
end

allUsed = find(trialInclude);
usedVals = xPre(allUsed);
usedSess = sessionPerTrial(allUsed);

sessionIds = P.sessions(:).';
nSess = numel(sessionIds);
sessionMeanTc = nan(nSess, numel(tb));
sessionSemTc = nan(nSess, numel(tb));
sessionN = zeros(nSess,1);
for iSess = 1:nSess
    idx = trialInclude & (sessionPerTrial == sessionIds(iSess));
    X = tcAll(idx, :);
    sessionN(iSess) = size(X, 1);
    if ~isempty(X)
        sessionMeanTc(iSess,:) = mean(X, 1, 'omitnan');
        sessionSemTc(iSess,:) = std(X, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(X), 1));
    end
end

D = struct();
D.tb = tb(:);
D.allUsed = allUsed(:);
D.usedVals = usedVals(:);
D.usedSess = usedSess(:);
D.sessionIds = sessionIds(:);
D.sessionMeanTc = sessionMeanTc;
D.sessionSemTc = sessionSemTc;
D.sessionN = sessionN;
end

function fig = make_review_figure_local(row, D, P)
fig = figure('Color', 'w', 'Name', sprintf('%s site %d session review', char(row.areaName), row.siteGlobal));
tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile;
hold(ax1, 'on');
cols = [0.15 0.35 0.75; 0.80 0.25 0.15];
for iSess = 1:numel(D.sessionIds)
    c = cols(iSess,:);
    add_sem_patch_local(ax1, D.tb, D.sessionMeanTc(iSess,:).', D.sessionSemTc(iSess,:).', c, 0.16);
    plot(ax1, D.tb, D.sessionMeanTc(iSess,:), '-', 'Color', c, 'LineWidth', 2.0);
end
yl = ylim(ax1);
line(ax1, [0 0], yl, 'Color', [0 0 0], 'LineStyle', '-', 'LineWidth', 1.0);
xlabel(ax1, 'Time from stimulus onset (ms)');
ylabel(ax1, 'Mean response (normMUA)');
title(ax1, sprintf('%s site %d | suggested exclude session %d', ...
    char(row.areaName), row.siteGlobal, row.suggestedExcludeDay), 'Interpreter', 'none');
legend(ax1, ...
    {sprintf('session %d (N=%d)', D.sessionIds(1), D.sessionN(1)), ...
     sprintf('session %d (N=%d)', D.sessionIds(2), D.sessionN(2))}, ...
    'Location', 'best', 'Box', 'off');
grid(ax1, 'on');

ax2 = nexttile;
hold(ax2, 'on');
trialAxis = 1:numel(D.allUsed);
isSess1 = D.usedSess == D.sessionIds(1);
isSess2 = D.usedSess == D.sessionIds(2);
scatter(ax2, trialAxis(isSess1), D.usedVals(isSess1), 12, ...
    'MarkerFaceColor', cols(1,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.38);
scatter(ax2, trialAxis(isSess2), D.usedVals(isSess2), 12, ...
    'MarkerFaceColor', cols(2,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.38);
sessChange = find(diff(D.usedSess(:)) ~= 0);
for i = 1:numel(sessChange)
    xline(ax2, sessChange(i) + 0.5, ':', 'Color', [0.70 0.70 0.70]);
end
xlabel(ax2, 'Included trial order');
ylabel(ax2, 'Prestim trial mean (normMUA)');
title(ax2, 'Prestim trial order (dotted lines = session change)');
grid(ax2, 'on');

annotation(fig, 'textbox', [0.05 0.95 0.90 0.045], ...
    'String', sprintf(['late SD ratio s2/s1 = %.3f | pre SD ratio s2/s1 = %.3f | ' ...
    'late evoked [%.3f %.3f] | pre mean [%.3f %.3f]'], ...
    row.lateSdRatio, row.preSdRatio, row.lateEvSess1, row.lateEvSess2, row.preMeanSess1, row.preMeanSess2), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'FontWeight', 'bold');
end

function add_sem_patch_local(ax, x, mu, sem, color, alphaVal)
x = x(:);
mu = mu(:);
sem = sem(:);
good = isfinite(x) & isfinite(mu) & isfinite(sem);
if ~any(good)
    return;
end
x = x(good);
mu = mu(good);
sem = sem(good);
patch(ax, [x; flipud(x)], [mu-sem; flipud(mu+sem)], color, ...
    'FaceAlpha', alphaVal, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

function resp = ask_review_decision_local(row, idx, nTotal, fig)
if ~isempty(fig) && isvalid(fig)
    figure(fig);
    drawnow;
end

msg = sprintf(['[%d/%d] %s site %d\n' ...
    'Suggested exclude session: %d\n\n' ...
    'Choose how to handle this site.'], ...
    idx, nTotal, char(row.areaName), row.siteGlobal, row.suggestedExcludeDay);

choice = questdlg(msg, 'Review Session Quality Site', ...
    'Keep', 'Excl sugg', 'Excl s1', 'Keep');

if isempty(choice)
    choice = 'Quit';
end

switch choice
    case 'Keep'
        resp = 'y';
    case 'Excl sugg'
        resp = 'n';
    case 'Excl s1'
        resp = '1';
    otherwise
        choice2 = questdlg(msg, 'Review Session Quality Site', ...
            'Excl s2', 'Quit', 'Back', 'Excl s2');
        if isempty(choice2)
            choice2 = 'Quit';
        end
        switch choice2
            case 'Excl s2'
                resp = '2';
            case 'Back'
                resp = ask_review_decision_local(row, idx, nTotal, fig);
            otherwise
                resp = 'q';
        end
end
end
