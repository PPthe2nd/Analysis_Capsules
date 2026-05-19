function OUT = Inspect_Prestim_TargetSide_Distribution_IT()
% INSPECT_PRESTIM_TARGETSIDE_DISTRIBUTION_IT
% Site-level diagnostic for pre-stimulus target-side activity in IT.
%
% This uses the saved output from Attention_TargetSide_Timecourse_IT and
% summarizes the pre-stimulus window (t < 0 ms) for each site:
%   - paired site means for the two target-side configurations
%   - histogram of the paired difference
%
% Note: this diagnostic is site-level, not trial-level, because the
% Resp_capsules_* files contain stimulus-wise means rather than single-trial
% responses.

%% Settings
P = struct();
P.Monkey = 1; % 1 = Nilson, 2 = Figaro
P.makeFigure = true;
P.useSavedTimecourse = true;

cfg = config();

if P.Monkey == 1
    monkeySuffix = "N";
else
    monkeySuffix = "F";
end

tcPath = fullfile(cfg.resultsDir, sprintf('Attention_TargetSide_Timecourse_IT_%s_d12.mat', char(monkeySuffix)));
hasSessionExclusions = ~isempty(site_session_exclusions(monkeySuffix));
if exist(tcPath, 'file') == 2 && P.useSavedTimecourse && ~hasSessionExclusions
    S = load(tcPath, 'OUT');
    assert(isfield(S, 'OUT') && isstruct(S.OUT), '%s must contain OUT.', tcPath);
    T = S.OUT;
else
    T = Attention_TargetSide_Timecourse_IT(struct('plotFigure', false, 'plotDifferenceFigure', false));
end

assert(isfield(T, 'Groups') && isfield(T, 'tCenters'), 'Timecourse OUT is missing Groups/tCenters.');
preMask = T.tCenters < 0;

G = T.Groups;
for i = 1:numel(G)
    cond1 = mean(G(i).tcPref(:, preMask), 2, 'omitnan');
    cond2 = mean(G(i).tcNon(:, preMask), 2, 'omitnan');
    diff = cond1 - cond2;
    keep = isfinite(cond1) & isfinite(cond2) & isfinite(diff);
    cond1 = cond1(keep);
    cond2 = cond2(keep);
    diff = diff(keep);

    Sg = struct();
    Sg.label = G(i).label;
    if isfield(G(i), 'cond1Label')
        Sg.cond1Label = string(G(i).cond1Label);
    else
        Sg.cond1Label = "Cond 1";
    end
    if isfield(G(i), 'cond2Label')
        Sg.cond2Label = string(G(i).cond2Label);
    else
        Sg.cond2Label = "Cond 2";
    end
    Sg.cond1 = cond1;
    Sg.cond2 = cond2;
    Sg.diff = diff;
    Sg.nSites = numel(diff);
    Sg.meanDiff = mean(diff, 'omitnan');
    Sg.medianDiff = median(diff, 'omitnan');
    Sg.trim20Diff = trimmean(diff, 20);
    Sg.fracPositive = mean(diff > 0, 'omitnan');
    if numel(diff) >= 3 && exist('signrank', 'file') == 2
        Sg.pSignrank = signrank(diff, 0);
    elseif numel(diff) >= 3
        [~, Sg.pSignrank] = ttest(diff, 0);
    else
        Sg.pSignrank = NaN;
    end
    Stats(i) = Sg; %#ok<AGROW>
end

OUT = struct();
OUT.P = P;
OUT.monkeySuffix = monkeySuffix;
OUT.preMask = preMask;
OUT.Stats = Stats;

fprintf('IT pre-stim site distributions (%s)\n', char(monkeySuffix));
for i = 1:numel(Stats)
    fprintf('  %s | N=%d | mean=%.4f | median=%.4f | trim20=%.4f | fracPos=%.3f | p=%.3g\n', ...
        Stats(i).label, Stats(i).nSites, Stats(i).meanDiff, Stats(i).medianDiff, ...
        Stats(i).trim20Diff, Stats(i).fracPositive, Stats(i).pSignrank);
end

if P.makeFigure
    make_distribution_figure(Stats, char(monkeySuffix));
end
end

function make_distribution_figure(Stats, monkeySuffix)
fig = figure('Color', 'w', ...
    'Name', sprintf('IT pre-stim target-side conditions (%s)', monkeySuffix), ...
    'Position', [60 100 1400 700]);
useTiled = exist('tiledlayout', 'file') == 2;
if useTiled
    tiledlayout(2, numel(Stats), 'TileSpacing', 'compact', 'Padding', 'compact');
end

cPref = [0.85 0.20 0.15];
cNon = [0.25 0.25 0.25];

for i = 1:numel(Stats)
    S = Stats(i);

    if useTiled
        ax1 = nexttile;
    else
        ax1 = subplot(2, numel(Stats), i); %#ok<LAXES>
    end
    hold(ax1, 'on');
    xCond1 = 1 + 0.06 * randn(S.nSites, 1);
    xCond2 = 2 + 0.06 * randn(S.nSites, 1);
    for k = 1:S.nSites
        plot(ax1, [xCond1(k) xCond2(k)], [S.cond1(k) S.cond2(k)], '-', ...
            'Color', [0.80 0.80 0.80], 'LineWidth', 0.8);
    end
    plot(ax1, xCond1, S.cond1, 'o', 'MarkerFaceColor', cPref, 'MarkerEdgeColor', 'none', 'MarkerSize', 4.5);
    plot(ax1, xCond2, S.cond2, 'o', 'MarkerFaceColor', cNon, 'MarkerEdgeColor', 'none', 'MarkerSize', 4.5);
    plot(ax1, [0.85 1.15], [mean(S.cond1,'omitnan') mean(S.cond1,'omitnan')], '-', 'Color', cPref, 'LineWidth', 2.5);
    plot(ax1, [1.85 2.15], [mean(S.cond2,'omitnan') mean(S.cond2,'omitnan')], '-', 'Color', cNon, 'LineWidth', 2.5);
    set(ax1, 'XLim', [0.6 2.4], 'XTick', [1 2], 'XTickLabel', cellstr([S.cond1Label; S.cond2Label]));
    ylabel(ax1, 'Pre-stim mean activity');
    title(ax1, sprintf('%s (N=%d)', S.label, S.nSites));
    grid(ax1, 'on');

    txt = sprintf('mean %.3f | med %.3f\ntrim20 %.3f | p %.2g', ...
        S.meanDiff, S.medianDiff, S.trim20Diff, S.pSignrank);
    text(ax1, 0.02, 0.98, txt, 'Units', 'normalized', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
        'FontSize', 8, 'Color', [0.25 0.25 0.25]);

    if useTiled
        ax2 = nexttile;
    else
        ax2 = subplot(2, numel(Stats), numel(Stats) + i); %#ok<LAXES>
    end
    hold(ax2, 'on');
    histogram(ax2, S.diff, 14, 'FaceColor', [0.55 0.55 0.55], 'EdgeColor', 'none');
    xline(ax2, 0, 'k:');
    xline(ax2, S.meanDiff, '-', 'Color', [0.85 0.20 0.15], 'LineWidth', 2);
    xline(ax2, S.medianDiff, '--', 'Color', [0.15 0.35 0.75], 'LineWidth', 1.8);
    xlabel(ax2, sprintf('Pre-stim %s - %s', char(S.cond1Label), char(S.cond2Label)));
    ylabel(ax2, 'N sites');
    grid(ax2, 'on');
    legend(ax2, {'', '', 'mean', 'median'}, 'Location', 'best', 'Box', 'off');
end

annotation(fig, 'textbox', [0.05 0.95 0.90 0.04], ...
    'String', sprintf(['IT pre-stimulus site distributions (%s). Top row: paired site means for the two target-side conditions. ' ...
    'Bottom row: histogram of the paired pre-stim difference. This figure is site-level, not trial-level.'], monkeySuffix), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'FontSize', 12, 'FontWeight', 'bold');
end
