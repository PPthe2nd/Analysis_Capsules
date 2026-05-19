function OUT = Show_TargetSide_Timecourse_Examples_IT(Puser)
% SHOW_TARGETSIDE_TIMECOURSE_EXAMPLES_IT
% Show example IT sites spanning strong, weak, and typical absolute
% response levels in the target-side timecourse groups from
% Attention_TargetSide_Timecourse_IT.
%
% For each target-side timecourse group, the script selects:
%   - 2 strong-response sites
%   - 1 typical-response site
%   - 2 weak-response sites
%
% The response level is the mean of the two absolute-response traces
% (target-closer and distractor-closer) over the requested response window.

P = struct();
P.Monkey = 1; % 1 = Nilson, 2 = Figaro
P.responseWindowMs = [50 500];
P.lateWindowMs = [300 500];
P.nStrong = 2;
P.nWeak = 2;
P.plotFigure = true;
P.saveResult = false;

if nargin >= 1 && ~isempty(Puser)
    fn = fieldnames(Puser);
    for i = 1:numel(fn)
        P.(fn{i}) = Puser.(fn{i});
    end
end

cfg = config();
if P.Monkey == 1
    monkeySuffix = "N";
else
    monkeySuffix = "F";
end

TIME = Attention_TargetSide_Timecourse_IT(struct( ...
    'Monkey', P.Monkey, ...
    'plotFigure', false, ...
    'plotDifferenceFigure', false, ...
    'plotTimingFitFigure', false));

t = TIME.tCenters(:);
respMask = (t >= P.responseWindowMs(1)) & (t <= P.responseWindowMs(2));
lateMask = (t >= P.lateWindowMs(1)) & (t <= P.lateWindowMs(2));
assert(any(respMask), 'No bins in response window [%d %d].', P.responseWindowMs(1), P.responseWindowMs(2));
assert(any(lateMask), 'No bins in late window [%d %d].', P.lateWindowMs(1), P.lateWindowMs(2));

nGroups = numel(TIME.Groups);
Examples = repmat(struct( ...
    'groupLabel', "", ...
    'responseWindowMs', P.responseWindowMs, ...
    'lateWindowMs', P.lateWindowMs, ...
    'siteLocal', zeros(0,1), ...
    'siteGlobal', zeros(0,1), ...
    'responseLevel', zeros(0,1), ...
    'kind', strings(0,1), ...
    'latePref', zeros(0,1), ...
    'lateNon', zeros(0,1), ...
    'lateDiff', zeros(0,1)), nGroups, 1);

for g = 1:nGroups
    G = TIME.Groups(g);
    siteLocalAll = find(G.mask);
    siteGlobalAll = TIME.RFrange(G.mask);
    nSites = numel(siteLocalAll);
    if nSites == 0
        Examples(g).groupLabel = string(G.label);
        continue;
    end

    pref = double(G.tcPref);
    non = double(G.tcNon);
    responseLevel = 0.5 * ( ...
        mean(pref(:, respMask), 2, 'omitnan') + ...
        mean(non(:, respMask), 2, 'omitnan'));
    responseLevel(~isfinite(responseLevel)) = NaN;

    [~, orderStrong] = sort(responseLevel, 'descend', 'MissingPlacement', 'last');
    [~, orderWeak] = sort(responseLevel, 'ascend', 'MissingPlacement', 'last');
    strongIdx = orderStrong(1:min(P.nStrong, nSites));

    weakIdx = zeros(0,1);
    for k = 1:numel(orderWeak)
        idx = orderWeak(k);
        if ~ismember(idx, strongIdx)
            weakIdx(end+1,1) = idx; %#ok<AGROW>
        end
        if numel(weakIdx) >= min(P.nWeak, max(0, nSites - numel(strongIdx) - 1))
            break;
        end
    end

    medLevel = median(responseLevel, 'omitnan');
    levelTmp = responseLevel;
    levelTmp(unique([strongIdx(:); weakIdx(:)])) = NaN;
    [~, iTypical] = min(abs(levelTmp - medLevel));
    if isempty(iTypical) || ~isfinite(levelTmp(iTypical))
        remaining = setdiff(find(isfinite(responseLevel)), unique([strongIdx(:); weakIdx(:)]), 'stable');
        if isempty(remaining)
            remaining = find(isfinite(responseLevel), 1, 'first');
        end
        iTypical = remaining(1);
    end

    pickIdx = [strongIdx(:); iTypical; weakIdx(:)];
    kinds = [ ...
        "strong 1"; ...
        "strong 2"; ...
        "typical"; ...
        "weak 1"; ...
        "weak 2"];
    kinds = kinds(1:numel(pickIdx));
    latePref = mean(pref(:, lateMask), 2, 'omitnan');
    lateNon = mean(non(:, lateMask), 2, 'omitnan');

    Examples(g).groupLabel = string(G.label);
    Examples(g).siteLocal = siteLocalAll(pickIdx);
    Examples(g).siteGlobal = siteGlobalAll(pickIdx);
    Examples(g).responseLevel = responseLevel(pickIdx);
    Examples(g).kind = kinds;
    Examples(g).latePref = latePref(pickIdx);
    Examples(g).lateNon = lateNon(pickIdx);
    Examples(g).lateDiff = latePref(pickIdx) - lateNon(pickIdx);
end

OUT = struct();
OUT.P = P;
OUT.monkeySuffix = monkeySuffix;
OUT.TIME = TIME;
OUT.Examples = Examples;

if P.saveResult
    outPath = fullfile(cfg.resultsDir, sprintf('Show_TargetSide_Timecourse_Examples_IT_%s.mat', char(monkeySuffix)));
    save(outPath, 'OUT', '-v7.3');
    fprintf('Saved IT target-side SEM example summary to %s\n', outPath);
end

if P.plotFigure
    make_examples_figure_local(OUT);
end
end

function make_examples_figure_local(OUT)
TIME = OUT.TIME;
t = TIME.tCenters(:);
Examples = OUT.Examples;
nGroups = numel(Examples);
nCols = max(cellfun(@numel, {Examples.siteLocal}));

cPref = [0.90 0.25 0.18];
cNon = [0.25 0.25 0.25];
cPrefMean = [1.00 0.75 0.72];
cNonMean = [0.65 0.65 0.65];

fig = figure('Color', 'w', 'Name', sprintf('IT target-side SEM examples (%s)', char(OUT.monkeySuffix)), ...
    'NumberTitle', 'off');
useTiled = exist('tiledlayout', 'file') == 2;
if useTiled
    tiledlayout(nGroups, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
end

for g = 1:nGroups
    G = TIME.Groups(g);
    ex = Examples(g);
    nThis = numel(ex.siteLocal);
    for j = 1:nCols
        if useTiled
            ax = nexttile;
        else
            ax = subplot(nGroups, nCols, (g-1)*nCols + j); %#ok<LAXES>
        end
        hold(ax, 'on');
        box(ax, 'off');
        grid(ax, 'on');
        xline(ax, 0, 'k-');
        yline(ax, 0, 'k:');
        xlabel(ax, 'Time from stimulus onset (ms)');
        ylabel(ax, 'Response (spont. SD units)');

        if j <= nThis
            siteLocal = ex.siteLocal(j);
            siteGlobal = ex.siteGlobal(j);
            rowInGroup = find(find(G.mask) == siteLocal, 1, 'first');
            pref = double(G.tcPref(rowInGroup,:));
            non = double(G.tcNon(rowInGroup,:));

            plot(ax, t, G.meanPref, '--', 'Color', cPrefMean, 'LineWidth', 1.5);
            plot(ax, t, G.meanNon, '--', 'Color', cNonMean, 'LineWidth', 1.5);
            plot(ax, t, pref, '-', 'Color', cPref, 'LineWidth', 1.7);
            plot(ax, t, non, '-', 'Color', cNon, 'LineWidth', 1.7);

            kindLabel = char(ex.kind(j));
            title(ax, sprintf('%s\nsite %d (%s)', short_group_label_local(G.label), siteGlobal, kindLabel), ...
                'FontSize', 10, 'Interpreter', 'none');

            yL = ylim(ax);
            xText = t(end) - 0.02 * range(t);
            yText = yL(1) + 0.08 * range(yL);
            text(ax, xText, yText, sprintf('resp=%.3f | late \\Delta=%.3f', ...
                ex.responseLevel(j), ex.lateDiff(j)), ...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
                'FontSize', 8.5, 'BackgroundColor', 'w', 'Margin', 1.5);

            if g == 1 && j == 1
                legend(ax, {'group mean pref', 'group mean non', 'site pref', 'site non'}, ...
                    'Location', 'northwest', 'Box', 'off');
            end
        else
            axis(ax, 'off');
        end
    end
end

annotation(fig, 'textbox', [0.05 0.955 0.90 0.04], ...
    'String', sprintf(['IT target-side example sites contributing to broad SEM bands (%s). ' ...
    'Each group shows strong-response sites, a typical site, and weak-response sites, using the mean of the two absolute-response traces in [%d, %d] ms.'], ...
    char(OUT.monkeySuffix), OUT.P.responseWindowMs(1), OUT.P.responseWindowMs(2)), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'FontSize', 11, 'FontWeight', 'bold');
end

function s = short_group_label_local(label)
lab = string(label);
lab = replace(lab, "Distance growth sig (late-pre)", "Distance growth sig");
lab = replace(lab, "Direction sig (dir|dist)", "Direction sig");
lab = replace(lab, "Distance growth only", "Distance only");
lab = replace(lab, "Direction only", "Direction only");
s = char(lab);
end
