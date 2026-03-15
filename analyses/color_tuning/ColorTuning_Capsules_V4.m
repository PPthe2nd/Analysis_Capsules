% Color tuning with complementary balancing for V4 sites in the line task.
%
% This mirrors the V1 ColorTuning_Capsules logic, but uses the V4 geometry
% tables and computes the histogram directly for the V4 onset-driven set.

Monkey = 1; % 1 = Nilson, 2 = Figaro
SNRthr = 0.7; % minimum onset SNR for inclusion in the color-tuning computation
alpha = 0.05; % significance threshold for the histogram overlay
SaveResults = true;

cfg = config();

if Monkey == 1
    monkeySuffix = "N";
    tallFile = 'Tall_V4_lines_N.mat';
    resp3binFile = 'SNR_capsules_N_d12.mat';
elseif Monkey == 2
    monkeySuffix = "F";
    tallFile = 'Tall_V4_lines_F.mat';
    resp3binFile = 'SNR_capsules_F_d12.mat';
else
    error('ColorTuning_Capsules_V4:InvalidMonkey', ...
        'Monkey must be 1 (Nilson) or 2 (Figaro).');
end

tallPath = fullfile(cfg.matDir, tallFile);
resp3binPath = fullfile(cfg.matDir, resp3binFile);
outFile = fullfile(cfg.matDir, sprintf('ColorTune_balanced_V4_%s.mat', char(monkeySuffix)));

assert(exist(tallPath, 'file') == 2, ...
    'Missing %s. Run Line_Stimuli_V4.m first.', tallPath);
assert(exist(resp3binPath, 'file') == 2, ...
    'Missing %s. Create the 3-bin response summary first.', resp3binPath);

Sgeo = load(tallPath);
assert(isfield(Sgeo, 'Tall_V4') && isstruct(Sgeo.Tall_V4), ...
    '%s must contain struct Tall_V4.', tallFile);
assert(isfield(Sgeo, 'RFrange') && ~isempty(Sgeo.RFrange), ...
    '%s must contain RFrange.', tallFile);

Tall_V4 = Sgeo.Tall_V4;
RFrange = Sgeo.RFrange(:);
nV4 = numel(RFrange);
siteRows = (1:nV4).';

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

SNR = compute_snr_per_color_sites(R3, Tall_V4, siteRows, 'Verbose', true);
SNRmat = [SNR.yellowEarly(siteRows), SNR.yellowLate(siteRows), ...
          SNR.purpleEarly(siteRows), SNR.purpleLate(siteRows)];
[bestSNR, ~] = max(SNRmat, [], 2, 'omitnan');

keepSiteIdx = find(isfinite(bestSNR) & (bestSNR > SNRthr));
fprintf('Keeping %d / %d V4 sites (bestSNR > %.2f)\n', numel(keepSiteIdx), nV4, SNRthr);
assert(~isempty(keepSiteIdx), ...
    'No V4 sites passed the SNR threshold %.2f.', SNRthr);

ColorTune = compute_color_tuning_balanced_sites(R3, Tall_V4, siteRows, keepSiteIdx, 'Verbose', true);
ColorTune.thr = SNRthr;
ColorTune.bestSNR = bestSNR;
ColorTune.RFrange = RFrange;
ColorTune.monkeySuffix = monkeySuffix;

ci = ColorTune.early.colorIndex(ColorTune.keepSites);
pvals = ColorTune.early.p(ColorTune.keepSites);
valid = isfinite(ci) & isfinite(pvals);
ci = ci(valid);
pvals = pvals(valid);
sig = pvals < alpha;

figure; hold on;
histogram(ci, 30, ...
    'FaceColor', [0.8 0.8 0.8], ...
    'EdgeColor', 'none');
histogram(ci(sig), 30, ...
    'FaceColor', [0.8 0 0], ...
    'EdgeColor', 'none');
xlabel('Color Index (yellow - purple)');
ylabel('Number of sites');
title(sprintf('V4 Color tuning (40-240ms, %s), SNR>%.2f', char(monkeySuffix), SNRthr));
legend('All sites', sprintf('Significant (p<%.2f)', alpha));
grid on;

if SaveResults
    save(outFile, 'ColorTune', '-v7.3');
    fprintf('Saved V4 color tuning results to %s\n', outFile);
end
