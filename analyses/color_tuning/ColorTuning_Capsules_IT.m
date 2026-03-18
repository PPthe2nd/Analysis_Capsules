% Color tuning with complementary balancing for IT sites in the line task.
%
% This mirrors the V4 color-tuning workflow, but uses Tall_IT so the
% yellow/purple label is read at the pre-existing RF center loaded via RFs.m
% when Tall_IT_lines_*.mat was built.

Monkey = 1; % 1 = Nilson, 2 = Figaro
SNRthr = 0.7; % minimum onset SNR for inclusion
alpha = 0.05; % significance threshold for the histogram overlay
SaveResults = true;

cfg = config();

if Monkey == 1
    monkeySuffix = "N";
    tallFile = 'Tall_IT_lines_N.mat';
    resp3binFile = 'SNR_capsules_N_d12.mat';
elseif Monkey == 2
    monkeySuffix = "F";
    tallFile = 'Tall_IT_lines_F.mat';
    resp3binFile = 'SNR_capsules_F_d12.mat';
else
    error('ColorTuning_Capsules_IT:InvalidMonkey', ...
        'Monkey must be 1 (Nilson) or 2 (Figaro).');
end

tallPath = fullfile(cfg.matDir, tallFile);
resp3binPath = fullfile(cfg.matDir, resp3binFile);
outFile = fullfile(cfg.matDir, sprintf('ColorTune_balanced_IT_%s.mat', char(monkeySuffix)));

assert(exist(tallPath, 'file') == 2, ...
    'Missing %s. Run Line_Stimuli_IT.m first.', tallPath);
assert(exist(resp3binPath, 'file') == 2, ...
    'Missing %s. Create the 3-bin response summary first.', resp3binPath);

Sgeo = load(tallPath);
assert(isfield(Sgeo, 'Tall_IT') && isstruct(Sgeo.Tall_IT), ...
    '%s must contain struct Tall_IT.', tallFile);
assert(isfield(Sgeo, 'RFrange') && ~isempty(Sgeo.RFrange), ...
    '%s must contain RFrange.', tallFile);

Tall_IT = Sgeo.Tall_IT;
RFrange = Sgeo.RFrange(:);
nIT = numel(RFrange);
siteRows = (1:nIT).';

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

SNR = compute_snr_per_color_sites(R3, Tall_IT, siteRows, 'Verbose', true);
SNRmat = [SNR.yellowEarly(siteRows), SNR.yellowLate(siteRows), ...
          SNR.purpleEarly(siteRows), SNR.purpleLate(siteRows)];
[bestSNR, ~] = max(SNRmat, [], 2, 'omitnan');

keepSiteIdx = find(isfinite(bestSNR) & (bestSNR > SNRthr));
fprintf('Keeping %d / %d IT sites (bestSNR > %.2f)\n', numel(keepSiteIdx), nIT, SNRthr);
assert(~isempty(keepSiteIdx), ...
    'No IT sites passed the SNR threshold %.2f.', SNRthr);

ColorTune = compute_color_tuning_balanced_sites(R3, Tall_IT, siteRows, keepSiteIdx, 'Verbose', true);
ColorTune.thr = SNRthr;
ColorTune.bestSNR = bestSNR;
ColorTune.RFrange = RFrange;
ColorTune.monkeySuffix = monkeySuffix;

ciEarlyAll = ColorTune.early.colorIndex(ColorTune.keepSites);
pEarlyAll = ColorTune.early.p(ColorTune.keepSites);
keepEarly = isfinite(ciEarlyAll);
ciEarly = ciEarlyAll(keepEarly);
sigEarlyMask = keepEarly & isfinite(pEarlyAll) & (pEarlyAll < alpha);
ciEarlySig = ciEarlyAll(sigEarlyMask);

ciLateAll = ColorTune.late.colorIndex(ColorTune.keepSites);
pLateAll = ColorTune.late.p(ColorTune.keepSites);
keepLate = isfinite(ciLateAll);
ciLate = ciLateAll(keepLate);
sigLateMask = keepLate & isfinite(pLateAll) & (pLateAll < alpha);
ciLateSig = ciLateAll(sigLateMask);

fprintf('Significant IT color tuning at RF center by pooled test (p < %.2f): early=%d late=%d\n', ...
    alpha, nnz(sigEarlyMask), nnz(sigLateMask));

figName = sprintf('RF.m CENTER | IT color tuning (%s)', char(monkeySuffix));
figTitle = sprintf('RF.m CENTER | IT color tuning (%s), SNR>%.2f', char(monkeySuffix), SNRthr);
figNum = 101;
figure(figNum); clf;
set(gcf, 'Color', 'w', 'Name', figName, 'NumberTitle', 'off', 'Tag', 'IT_RFm_center_color_tuning');
fprintf('Opened figure %d: %s\n', figNum, figTitle);
useTiled = exist('tiledlayout', 'file') == 2;
if useTiled
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
end

if useTiled, nexttile; else, subplot(1, 2, 1); end
hEarly = histogram(ciEarly, 30, 'FaceColor', [0.82 0.82 0.82], 'EdgeColor', 'none');
hold on;
histogram(ciEarlySig, 'BinEdges', hEarly.BinEdges, 'FaceColor', [0.80 0.10 0.10], 'EdgeColor', 'none');
xline(0, 'k-');
xlabel('Color index (yellow - purple)');
ylabel('N sites');
title('Early');
legend('All sites', sprintf('Significant (p<%.2f)', alpha));
grid on;

if useTiled, nexttile; else, subplot(1, 2, 2); end
hLate = histogram(ciLate, 30, 'FaceColor', [0.82 0.82 0.82], 'EdgeColor', 'none');
hold on;
histogram(ciLateSig, 'BinEdges', hLate.BinEdges, 'FaceColor', [0.80 0.10 0.10], 'EdgeColor', 'none');
xline(0, 'k-');
xlabel('Color index (yellow - purple)');
ylabel('N sites');
title('Late');
legend('All sites', sprintf('Significant (p<%.2f)', alpha));
grid on;

if exist('sgtitle', 'file') == 2
    sgtitle(figTitle);
else
    annotation('textbox', [0.12 0.955 0.76 0.04], 'String', figTitle, ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end

if SaveResults
    save(outFile, 'ColorTune', '-v7.3');
    fprintf('Saved IT color tuning results to %s\n', outFile);
end
