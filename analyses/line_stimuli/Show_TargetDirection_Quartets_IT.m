function OUT = Show_TargetDirection_Quartets_IT()
% SHOW_TARGETDIRECTION_QUARTETS_IT
% Illustrate a few IT sites with strong target-direction tuning by showing
% the 5 quartets with the strongest late attention-direction effect.
%
% For each selected site:
% - choose quartets by |deltaQuartetLate|
% - show one representative stimulus per quartet
% - display the member of each opposite-direction pair that drives the site
%   more strongly, while showing the full signed 0..360 deg direction view
% - render target arm light gray and distractor arm dark gray so color does
%   not visually drive the impression
% - annotate each panel with the two opposite directions and the signed
%   quartet-level attention-direction effect

%% Settings
P = struct();
P.Monkey = 1;            % 1 = Nilson, 2 = Figaro
P.nSites = 4;
P.nDirNotDistExtras = 2;
P.nQuartetsPerSite = 5;
P.timeBin = 3;           % 2 = early, 3 = late
P.requireSig = true;
P.makeFigures = true;
P.siteSelectionMode = 'cardinal_plus_dirnotdist';   % 'cardinal', 'topve', or 'cardinal_plus_dirnotdist'
P.cardinalTargetsDeg = [90 180 270 0];  % up, left, down, right

cfg = config();

if P.Monkey == 1
    monkeySuffix = "N";
    tallFile = 'Tall_IT_lines_N.mat';
    resp3binFile = 'SNR_capsules_N_d12.mat';
elseif P.Monkey == 2
    monkeySuffix = "F";
    tallFile = 'Tall_IT_lines_F.mat';
    resp3binFile = 'SNR_capsules_F_d12.mat';
else
    error('Show_TargetDirection_Quartets_IT:InvalidMonkey', ...
        'P.Monkey must be 1 (Nilson) or 2 (Figaro).');
end

basePath = fullfile(cfg.matDir, sprintf('Attention_TargetSide_Tuning_IT_directiondelta_%s.mat', char(monkeySuffix)));
distCtrlPath = fullfile(cfg.matDir, sprintf('Attention_TargetSide_DistanceControl_IT_%s.mat', char(monkeySuffix)));
tallPath = fullfile(cfg.matDir, tallFile);
respPath = fullfile(cfg.matDir, resp3binFile);

assert(exist(basePath, 'file') == 2, ...
    'Missing %s. Run Attention_TargetSide_Tuning_IT first.', basePath);
if strcmpi(P.siteSelectionMode, 'cardinal_plus_dirnotdist')
    assert(exist(distCtrlPath, 'file') == 2, ...
        'Missing %s. Run Attention_TargetSide_DistanceControl_IT first.', distCtrlPath);
end
assert(exist(tallPath, 'file') == 2, ...
    'Missing %s. Run Line_Stimuli_IT first.', tallPath);
assert(exist(respPath, 'file') == 2, ...
    'Missing %s. Create the 3-bin response summary first.', respPath);

%% Load required data
Sbase = load(basePath, 'OUT');
BASE = Sbase.OUT;
assert(isfield(BASE, 'FitLate') && isfield(BASE, 'QuartetTable') && ...
    isfield(BASE, 'deltaQuartetLate') && isfield(BASE, 'RFrange'), ...
    'Base OUT is missing required fields.');

Stall = load(tallPath, 'ALLCOORDS', 'RTAB384', 'RFrange', 'Tall_IT');
assert(isfield(Stall, 'ALLCOORDS') && isfield(Stall, 'RTAB384') && isfield(Stall, 'RFrange') && ...
    isfield(Stall, 'Tall_IT'), '%s must contain ALLCOORDS, RTAB384, RFrange, and Tall_IT.', tallPath);
ALLCOORDS = Stall.ALLCOORDS;
RTAB384 = Stall.RTAB384;
RFrange = Stall.RFrange(:);
Tall_IT = Stall.Tall_IT;

Sresp = load(respPath, 'R');
assert(isfield(Sresp, 'R') && isstruct(Sresp.R), '%s must contain struct R.', respPath);
R3_full = Sresp.R;
R3.meanAct = R3_full.meanAct(RFrange, :, :);
R3.meanSqAct = R3_full.meanSqAct(RFrange, :, :);
if isvector(R3_full.nTrials)
    R3.nTrials = R3_full.nTrials(:);
elseif ismatrix(R3_full.nTrials) && size(R3_full.nTrials,1) >= max(RFrange) && ...
        size(R3_full.nTrials,2) == size(R3_full.meanAct,2)
    R3.nTrials = R3_full.nTrials(RFrange, :);
else
    R3.nTrials = R3_full.nTrials;
end

siteRows = (1:numel(RFrange)).';
SNR = compute_snr_per_color_sites(R3, Tall_IT, siteRows, 'Verbose', false);
muSpont = SNR.muSpont(siteRows);
sdSpont = SNR.sdSpont(siteRows);

QuartetTable = BASE.QuartetTable;
thetaA = double(QuartetTable.targetDirDeg(:));
stimRepA = double(QuartetTable.stimRef(:));
stimRepB = stimRepA + 1;
deltaLate = BASE.deltaQuartetLate;
FitLate = BASE.FitLate;

veLate = [FitLate.r2TrainPct].';
pLate = [FitLate.pValueApprox].';
keep = isfinite(veLate);
if P.requireSig
    keep = keep & isfinite(pLate) & (pLate < 0.05);
end
siteIdx = find(keep);
assert(~isempty(siteIdx), 'No IT sites pass the selected late target-direction criterion.');
prefLate = mod([FitLate.prefDeg].', 360);
selectionLabel = strings(0,1);
if strcmpi(P.siteSelectionMode, 'cardinal')
    siteSel = select_sites_by_cardinal_direction(siteIdx, prefLate, veLate, P.cardinalTargetsDeg);
    siteSel = siteSel(1:min(P.nSites, numel(siteSel)));
    selectionLabel = repmat("cardinal", numel(siteSel), 1);
elseif strcmpi(P.siteSelectionMode, 'cardinal_plus_dirnotdist')
    siteSelCard = select_sites_by_cardinal_direction(siteIdx, prefLate, veLate, P.cardinalTargetsDeg);
    siteSelCard = siteSelCard(1:min(P.nSites, numel(siteSelCard)));
    selectionLabelCard = repmat("cardinal", numel(siteSelCard), 1);

    Sdist = load(distCtrlPath, 'OUT');
    D = Sdist.OUT;
    assert(isfield(D, 'ctrlSig') && isfield(D, 'distGivenDirSig') && isfield(D, 'uniqueR2Dir'), ...
        '%s must contain ctrlSig, distGivenDirSig, and uniqueR2Dir.', distCtrlPath);
    extraPool = find(D.ctrlSig & ~D.distGivenDirSig & isfinite(D.uniqueR2Dir) & isfinite(veLate));
    extraPool = setdiff(extraPool(:), siteSelCard(:), 'stable');
    extraSel = select_sites_by_unique_direction(extraPool, D.uniqueR2Dir, veLate, P.nDirNotDistExtras);
    selectionLabelExtra = repmat("dir-not-dist", numel(extraSel), 1);

    siteSel = [siteSelCard(:); extraSel(:)];
    selectionLabel = [selectionLabelCard; selectionLabelExtra];
else
    [~, ordSites] = sort(veLate(siteIdx), 'descend');
    siteSel = siteIdx(ordSites(1:min(P.nSites, numel(ordSites))));
    selectionLabel = repmat("topve", numel(siteSel), 1);
end

OUT = struct();
OUT.P = P;
OUT.monkeySuffix = monkeySuffix;
OUT.siteSel = siteSel;
OUT.globalSiteSel = RFrange(siteSel);
OUT.prefDegSel = prefLate(siteSel);
OUT.selectionLabel = selectionLabel;

for ii = 1:numel(siteSel)
    s = siteSel(ii);
    deltaS = deltaLate(s, :).';
    okQ = isfinite(deltaS) & isfinite(thetaA);
    assert(any(okQ), 'Selected site %d has no finite quartets.', s);
    [~, ordQ] = sort(abs(deltaS(okQ)), 'descend');
    idxQall = find(okQ);
    qSel = idxQall(ordQ(1:min(P.nQuartetsPerSite, numel(ordQ))));

    Site(ii).localSite = s; %#ok<AGROW>
    Site(ii).globalSiteInR = RFrange(s); %#ok<AGROW>
    Site(ii).selectionLabel = selectionLabel(ii); %#ok<AGROW>
    Site(ii).fitLate = FitLate(s); %#ok<AGROW>
    Site(ii).quartetIdx = qSel(:); %#ok<AGROW>
    Site(ii).stimShown = nan(numel(qSel), 1); %#ok<AGROW>
    Site(ii).shownDirDeg = nan(numel(qSel), 1); %#ok<AGROW>
    Site(ii).deltaQuartetLate = nan(numel(qSel), 1); %#ok<AGROW>
    Site(ii).prefPairMeanLate = nan(numel(qSel), 1); %#ok<AGROW>
    Site(ii).oppPairMeanLate = nan(numel(qSel), 1); %#ok<AGROW>
    Site(ii).rfXpx = nan(numel(qSel), 1); %#ok<AGROW>
    Site(ii).rfYpx = nan(numel(qSel), 1); %#ok<AGROW>
    Site(ii).thetaDegAll = thetaA; %#ok<AGROW>
    Site(ii).deltaQuartetAll = deltaS; %#ok<AGROW>

    for jj = 1:numel(qSel)
        q = qSel(jj);
        thetaCanon = mod(thetaA(q), 360);
        pairA = [stimRepA(q), stimRepA(q) + 4];
        pairB = [stimRepA(q) + 1, stimRepA(q) + 5];
        if deltaS(q) >= 0
            stimShow = stimRepA(q);
            dirShow = thetaCanon;
            deltaShow = deltaS(q);
            pairPref = pairA;
            pairOpp = pairB;
        else
            stimShow = stimRepB(q);
            dirShow = mod(thetaCanon + 180, 360);
            deltaShow = -deltaS(q);
            pairPref = pairB;
            pairOpp = pairA;
        end

        prefPairMean = pair_mean_norm(R3, s, pairPref, P.timeBin, muSpont(s), sdSpont(s));
        oppPairMean = pair_mean_norm(R3, s, pairOpp, P.timeBin, muSpont(s), sdSpont(s));

        Site(ii).stimShown(jj) = stimShow;
        Site(ii).shownDirDeg(jj) = dirShow;
        Site(ii).deltaQuartetLate(jj) = deltaShow;
        Site(ii).prefPairMeanLate(jj) = prefPairMean;
        Site(ii).oppPairMeanLate(jj) = oppPairMean;
        Site(ii).rfXpx(jj) = Tall_IT(stimShow).T.x_px(s);
        Site(ii).rfYpx(jj) = Tall_IT(stimShow).T.y_px(s);
    end

    [~, ordShown] = sort(Site(ii).shownDirDeg, 'ascend');
    Site(ii).quartetIdx = Site(ii).quartetIdx(ordShown);
    Site(ii).stimShown = Site(ii).stimShown(ordShown);
    Site(ii).shownDirDeg = Site(ii).shownDirDeg(ordShown);
    Site(ii).deltaQuartetLate = Site(ii).deltaQuartetLate(ordShown);
    Site(ii).prefPairMeanLate = Site(ii).prefPairMeanLate(ordShown);
    Site(ii).oppPairMeanLate = Site(ii).oppPairMeanLate(ordShown);
    Site(ii).rfXpx = Site(ii).rfXpx(ordShown);
    Site(ii).rfYpx = Site(ii).rfYpx(ordShown);

    if P.makeFigures
        make_site_figure(ALLCOORDS, RTAB384, Site(ii), char(monkeySuffix));
    end
end

OUT.Site = Site;
end

function siteSel = select_sites_by_cardinal_direction(siteIdx, prefDeg, veLate, targetDeg)
targetDeg = targetDeg(:);
siteIdx = siteIdx(:);
prefUse = prefDeg(siteIdx);
veUse = veLate(siteIdx);
chosen = false(size(siteIdx));
siteSel = nan(numel(targetDeg), 1);

for k = 1:numel(targetDeg)
    d = abs(local_angdiff_deg(prefUse, targetDeg(k)));
    inSector = d <= 45 & ~chosen & isfinite(veUse);
    if any(inSector)
        cand = find(inSector);
        [~, ord] = sort(veUse(cand), 'descend');
        pick = cand(ord(1));
    else
        cand = find(~chosen & isfinite(veUse) & isfinite(d));
        if isempty(cand)
            continue;
        end
        score = [d(cand), -veUse(cand)];
        [~, ord] = sortrows(score, [1 2]);
        pick = cand(ord(1));
    end
    chosen(pick) = true;
    siteSel(k) = siteIdx(pick);
end

siteSel = siteSel(isfinite(siteSel));
end

function siteSel = select_sites_by_unique_direction(sitePool, uniqueR2Dir, veLate, nPick)
sitePool = sitePool(:);
if isempty(sitePool) || nPick <= 0
    siteSel = zeros(0,1);
    return;
end
score = [-uniqueR2Dir(sitePool), -veLate(sitePool)];
[~, ord] = sortrows(score, [1 2]);
siteSel = sitePool(ord(1:min(nPick, numel(ord))));
end

function d = local_angdiff_deg(a, b)
d = mod((a - b) + 180, 360) - 180;
end

function pairMeanNorm = pair_mean_norm(R3, iSite, stimSet, timeBin, muSpont, sdSpont)
pairMeanNorm = NaN;
if ~(isfinite(muSpont) && isfinite(sdSpont) && sdSpont > 0)
    return;
end

if isvector(R3.nTrials)
    nUse = double(R3.nTrials(stimSet));
elseif ismatrix(R3.nTrials) && size(R3.nTrials,1) >= iSite && size(R3.nTrials,2) >= max(stimSet)
    nUse = double(R3.nTrials(iSite, stimSet));
else
    nUse = NaN(size(stimSet));
end
nUse = nUse(:);
muUse = squeeze(double(R3.meanAct(iSite, stimSet, timeBin)));
muUse = muUse(:);
good = isfinite(nUse) & (nUse > 0) & isfinite(muUse);
if ~any(good)
    return;
end

nUse = nUse(good);
muUse = muUse(good);
muPair = sum(nUse(:) .* muUse(:)) / sum(nUse(:));
pairMeanNorm = (muPair - muSpont) / sdSpont;
end

function make_site_figure(ALLCOORDS, RTAB384, S, monkeySuffix)
nQ = numel(S.quartetIdx);
figW = max(1500, 245 * nQ + 360);
figH = 500;
prefDeg = mod(S.fitLate.prefDeg, 360);
fig = figure('Color', 'w', 'Name', sprintf('IT target-direction quartets site %d', S.globalSiteInR), ...
    'NumberTitle', 'off', 'Position', [80 120 figW figH]);

left = 0.035;
right = 0.025;
qGap = 0.014;
fitGap = 0.032;
fitW = 0.17;
panelY = 0.20;
panelH = 0.62;
quartetW = (1 - left - right - fitW - fitGap - qGap * max(nQ - 1, 0)) / max(nQ, 1);
pairScale = compute_site_bar_scale(S.prefPairMeanLate, S.oppPairMeanLate);

for j = 1:nQ
    x0 = left + (j - 1) * (quartetW + qGap);
    ax = axes('Parent', fig, 'Position', [x0 panelY quartetW panelH]);
    [img, H, W, footerH, xRF, yRF, xFig, yFig, xBack, yBack] = render_target_gray_stim( ...
        ALLCOORDS, RTAB384, S.stimShown(j), S.rfXpx(j), S.rfYpx(j));
    imshow(img, 'Parent', ax, 'InitialMagnification', 'fit');
    axis(ax, 'image');
    axis(ax, 'off');
    hold(ax, 'on');
    rectangle(ax, 'Position', [1 1 W-1 H-1], 'EdgeColor', [0.78 0.78 0.78], 'LineWidth', 1);
    draw_pref_arrow(ax, W, H, prefDeg);
    draw_between_capsules_arrow(ax, xBack, yBack, xFig, yFig);
    plot(ax, xRF, yRF, 'o', ...
        'MarkerSize', 6, 'MarkerFaceColor', [0.95 0.15 0.15], ...
        'MarkerEdgeColor', 'w', 'LineWidth', 0.8);
    plot(ax, [1 W], [H+1 H+1], '-', 'Color', [0.88 0.88 0.88], 'LineWidth', 1);

    title(ax, sprintf('stim %d', S.stimShown(j)), ...
        'FontSize', 10, 'FontWeight', 'bold');
    footerFrac = footerH / (H + footerH);
    oppDirDeg = mod(S.shownDirDeg(j) + 180, 360);
    text(ax, 0.5, 0.11 + footerFrac/2, sprintf('%.0f%c vs %.0f%c | \\Delta_q = +%.2f', ...
        S.shownDirDeg(j), char(176), oppDirDeg, char(176), S.deltaQuartetLate(j)), ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 10, 'FontWeight', 'bold', ...
        'Interpreter', 'tex', 'Color', [0.15 0.15 0.15]);
    draw_pair_bar_footer(ax, W, H, footerH, S.prefPairMeanLate(j), S.oppPairMeanLate(j), ...
        pairScale, j == nQ);
end

xFit = 1 - right - fitW;
paxFit = polaraxes('Parent', fig, 'Position', [xFit panelY fitW panelH]);
plot_fit_panel(paxFit, S);

annotation(fig, 'textbox', [0.18 0.87 0.64 0.07], ...
    'String', sprintf('IT site %d (%s) | %s | late VE %.1f%% | pref %.0f%c', ...
        S.globalSiteInR, monkeySuffix, char(S.selectionLabel), S.fitLate.r2TrainPct, prefDeg, char(176)), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'FontSize', 17, 'FontWeight', 'bold', ...
    'Color', [0.08 0.08 0.08]);

annotation(fig, 'textbox', [0.035 0.045 0.72 0.04], ...
    'String', 'light = stronger direction of the shown quartet, dark = opposite direction, red dot = RF center, yellow arrow = fitted preferred direction, black arrow = shown pair direction, yellow points = same quartets shown as symmetric theta/+Delta and theta+180/-Delta pairs; black curve = circular moving average; dashed ring = zero', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'middle', 'FontSize', 9, 'Color', [0.35 0.35 0.35]);
end

function draw_pref_arrow(ax, W, H, prefDeg)
L = 0.16 * min(W, H);
x0 = 0.16 * W;
y0 = 0.18 * H;
dx = L * cosd(prefDeg);
dy = -L * sind(prefDeg); % image y-axis points downward

quiver(ax, x0, y0, dx, dy, 0, ...
    'Color', 'w', 'LineWidth', 3.0, 'MaxHeadSize', 0.9);
quiver(ax, x0, y0, dx, dy, 0, ...
    'Color', [0.98 0.88 0.12], 'LineWidth', 1.8, 'MaxHeadSize', 0.9);
end

function draw_between_capsules_arrow(ax, xBack, yBack, xFig, yFig)
v = [xFig - xBack, yFig - yBack];
if ~all(isfinite(v)) || norm(v) < 1
    return;
end
mid = 0.5 * ([xBack, yBack] + [xFig, yFig]);
seg = 0.22 * v;
x0 = mid(1) - 0.5 * seg(1);
y0 = mid(2) - 0.5 * seg(2);

quiver(ax, x0, y0, seg(1), seg(2), 0, ...
    'Color', 'w', 'LineWidth', 2.6, 'MaxHeadSize', 1.2);
quiver(ax, x0, y0, seg(1), seg(2), 0, ...
    'Color', [0.12 0.12 0.12], 'LineWidth', 1.4, 'MaxHeadSize', 1.2);
end

function draw_pair_bar_footer(ax, W, H, footerH, prefMean, oppMean, scale, showScaleBar)
barTop = H + 28;
barBottom = H + footerH - 8;
yZero = round((barTop + barBottom) / 2);
xPref = 0.43 * W;
xOpp = 0.57 * W;
barW = 0.08 * W;

plot(ax, [0.30*W 0.70*W], [yZero yZero], '-', 'Color', [0.78 0.78 0.78], 'LineWidth', 0.8);
draw_footer_bar(ax, xPref, yZero, barW, prefMean, scale, barTop, barBottom, [0.82 0.82 0.82]);
draw_footer_bar(ax, xOpp, yZero, barW, oppMean, scale, barTop, barBottom, [0.32 0.32 0.32]);

text(ax, xPref, barBottom + 2, 'pref', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', 'FontSize', 8.5, 'Color', [0.25 0.25 0.25]);
text(ax, xOpp, barBottom + 2, 'opp', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', 'FontSize', 8.5, 'Color', [0.25 0.25 0.25]);

if showScaleBar
    draw_footer_scale_bar(ax, W, yZero, barTop, barBottom, scale);
end
end

function draw_footer_bar(ax, xCenter, yZero, barW, value, scale, barTop, barBottom, faceColor)
if ~(isfinite(value) && isfinite(scale) && scale > 0)
    return;
end
barHalfRange = 0.44 * (barBottom - barTop);
barH = barHalfRange * (value / scale);
if barH >= 0
    yRect = yZero - barH;
    hRect = barH;
else
    yRect = yZero;
    hRect = -barH;
end
rectangle(ax, 'Position', [xCenter - barW/2, yRect, barW, max(hRect, 1)], ...
    'FaceColor', faceColor, 'EdgeColor', [0.25 0.25 0.25], 'LineWidth', 0.8);
end

function draw_footer_scale_bar(ax, W, yZero, barTop, barBottom, scale)
if ~(isfinite(scale) && scale > 0)
    return;
end

xS = 0.83 * W;
yTop = barTop + 2;
yBot = barBottom - 2;

plot(ax, [xS xS], [yTop yBot], '-', 'Color', [0.35 0.35 0.35], 'LineWidth', 0.9);
plot(ax, [xS-6 xS+6], [yTop yTop], '-', 'Color', [0.35 0.35 0.35], 'LineWidth', 0.9);
plot(ax, [xS-6 xS+6], [yZero yZero], '-', 'Color', [0.35 0.35 0.35], 'LineWidth', 0.9);
plot(ax, [xS-6 xS+6], [yBot yBot], '-', 'Color', [0.35 0.35 0.35], 'LineWidth', 0.9);

text(ax, xS + 10, yTop, sprintf('+%.1f', scale), ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'FontSize', 8.5, 'Color', [0.25 0.25 0.25]);
text(ax, xS + 10, yZero, '0', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'FontSize', 8.5, 'Color', [0.25 0.25 0.25]);
text(ax, xS + 10, yBot, sprintf('-%.1f', scale), ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'FontSize', 8.5, 'Color', [0.25 0.25 0.25]);
text(ax, xS - 4, barBottom + 2, 'SD', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
    'FontSize', 8, 'Color', [0.25 0.25 0.25]);
end

function scale = compute_site_bar_scale(prefVals, oppVals)
vals = abs([prefVals(:); oppVals(:)]);
vals = vals(isfinite(vals));
if isempty(vals)
    scale = 0.75;
else
    scale = max([0.75; vals]);
end
end

function plot_fit_panel(ax, S)
theta = double(S.thetaDegAll(:));
delta = double(S.deltaQuartetAll(:));
smoothWin = 10;

ok = isfinite(theta) & isfinite(delta);
thetaOk = mod(theta(ok), 360);
deltaOk = delta(ok);

thetaPlotAll = [thetaOk; mod(thetaOk + 180, 360)];
deltaPlotAll = [deltaOk; -deltaOk];
thetaShow = mod(double(S.shownDirDeg(:)), 360);
deltaShow = double(S.deltaQuartetLate(:));
thetaShowAll = [thetaShow; mod(thetaShow + 180, 360)];
deltaShowAll = [deltaShow; -deltaShow];

thetaGrid = linspace(0, 360, 361).';
fitY = harmonic_curve(thetaGrid, S.fitLate);
maxAbs = max(abs([deltaPlotAll; deltaShowAll; fitY]));
if ~(isfinite(maxAbs) && maxAbs > 0)
    maxAbs = 1;
end

r0 = 1.25 * maxAbs + 0.05;
rMin = max(0, r0 - 1.10 * maxAbs);
rMax = r0 + 1.10 * maxAbs;
[thetaSmooth, deltaSmooth] = circular_moving_average_curve(thetaPlotAll, deltaPlotAll, smoothWin);

hold(ax, 'on');
polarplot(ax, deg2rad(linspace(0, 360, 361)), r0 * ones(361,1), '--', ...
    'Color', [0.65 0.65 0.65], 'LineWidth', 1.0);
polarplot(ax, deg2rad([S.fitLate.prefDeg S.fitLate.prefDeg]), [rMin rMax], '--', ...
    'Color', [0.20 0.20 0.20], 'LineWidth', 1.0);
polarplot(ax, deg2rad([mod(S.fitLate.prefDeg + 180, 360) mod(S.fitLate.prefDeg + 180, 360)]), [rMin rMax], '--', ...
    'Color', [0.20 0.20 0.20], 'LineWidth', 1.0);
polarplot(ax, deg2rad(thetaPlotAll), r0 + deltaPlotAll, 'o', ...
    'MarkerSize', 4.2, 'MarkerFaceColor', [0.72 0.72 0.72], ...
    'MarkerEdgeColor', 'none', 'Color', [0.72 0.72 0.72]);
polarplot(ax, deg2rad(thetaShowAll), r0 + deltaShowAll, 'o', ...
    'MarkerSize', 5.6, 'MarkerFaceColor', [0.98 0.88 0.12], ...
    'MarkerEdgeColor', [0.20 0.20 0.20], 'LineWidth', 0.8, ...
    'Color', [0.98 0.88 0.12]);
polarplot(ax, deg2rad(thetaSmooth), r0 + deltaSmooth, 'k-', 'LineWidth', 1.8);

ax.ThetaLim = [0 360];
ax.ThetaTick = [0 90 180 270];
ax.ThetaZeroLocation = 'right';
ax.ThetaDir = 'counterclockwise';
ax.RLim = [rMin rMax];
ax.RAxisLocation = 135;
ax.RTick = [r0 - maxAbs, r0, r0 + maxAbs];
ax.RTickLabel = {sprintf('%.1f', -maxAbs), '0', sprintf('%.1f', maxAbs)};
ax.GridAlpha = 0.22;
ax.MinorGridAlpha = 0.12;
ax.ThetaColor = [0.25 0.25 0.25];
ax.RColor = [0.25 0.25 0.25];
title(ax, 'Late Polar Summary', 'FontSize', 10, 'FontWeight', 'bold');
end

function yHat = harmonic_curve(thetaDeg, F)
yHat = F.baseline + F.coefCos * cosd(thetaDeg) + F.coefSin * sind(thetaDeg);
end

function [thetaCurve, deltaCurve] = circular_moving_average_curve(thetaDeg, deltaVal, winN)
thetaDeg = mod(thetaDeg(:), 360);
deltaVal = deltaVal(:);
good = isfinite(thetaDeg) & isfinite(deltaVal);
thetaDeg = thetaDeg(good);
deltaVal = deltaVal(good);

if isempty(thetaDeg)
    thetaCurve = linspace(0, 360, 361).';
    deltaCurve = zeros(size(thetaCurve));
    return;
end

[thetaSort, ord] = sort(thetaDeg);
deltaSort = deltaVal(ord);
n = numel(thetaSort);
thetaExt = [thetaSort - 360; thetaSort; thetaSort + 360];
deltaExt = [deltaSort; deltaSort; deltaSort];
deltaSmExt = movmean(deltaExt, winN, 'Endpoints', 'shrink');
deltaSm = deltaSmExt(n+1:2*n);

thetaCurve = [thetaSort; thetaSort(1) + 360];
deltaCurve = [deltaSm; deltaSm(1)];
end

function [imgGray, H, W, footerH, xRFout, yRFout, xFigOut, yFigOut, xBackOut, yBackOut] = render_target_gray_stim(ALLCOORDS, RTAB384, stimNum, xRFpx, yRFpx)
[~, masks, meta] = render_stim_with_masks2(ALLCOORDS, RTAB384, stimNum, 'DrawDots', false);
H0 = meta.imageSize(2);
W0 = meta.imageSize(1);
dispW = 320;
dispH = 240;
footerH = 96;

bg = 0.55;
targetGray = 0.82;
distrGray = 0.32;

img0 = repmat(bg, [H0 W0 3]);
for c = 1:3
    chan = img0(:,:,c);
    chan(masks.backArm) = distrGray;
    chan(masks.figArm) = targetGray;
    img0(:,:,c) = chan;
end

if exist('imresize', 'file') == 2
    imgDisp = imresize(img0, [dispH dispW], 'bilinear');
else
    imgDisp = img0;
    dispH = H0;
    dispW = W0;
end

imgGray = ones(dispH + footerH, dispW, 3);
imgGray(1:dispH,:,:) = imgDisp;

xRFout = xRFpx * (dispW / W0);
yRFout = yRFpx * (dispH / H0);
tFig = double(meta.t_fig(:))';
tBack = double(meta.t_back(:))';
toImg = @(p) [p(1) + W0/2, H0/2 - p(2)];
pFig = toImg(tFig);
pBack = toImg(tBack);
xFigOut = pFig(1) * (dispW / W0);
yFigOut = pFig(2) * (dispH / H0);
xBackOut = pBack(1) * (dispW / W0);
yBackOut = pBack(2) * (dispH / H0);
H = dispH;
W = dispW;
end
