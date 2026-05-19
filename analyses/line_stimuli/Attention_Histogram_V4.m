% Attention_Histogram_V4
% Histogram of V4 attention indices in the line task, matching the V1 plot:
% - gray: all valid sites
% - red: sites with pTD < pThresh
%
% Attention index definition:
%   (T - D) / abs(mean(T, D))

Monkey = 1; % 1 = Nilson, 2 = Figaro
timeIdx3 = 3; % 3-bin dataset: 300-500 ms
pThresh = 0.05; % significance threshold for the red overlay
idxClip = 2; % clip only for display, matching the V1 histogram style
excludeOverlap = true;
nHistBins = 30;
epsDen = 1e-6; % denominator floor used by attention_modulation_V1_3bin

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
    error('Attention_Histogram_V4:InvalidMonkey', ...
        'Monkey must be 1 (Nilson) or 2 (Figaro).');
end

tallPath = fullfile(cfg.matDir, tallFile);
resp3binPath = fullfile(cfg.matDir, resp3binFile);

assert(exist(tallPath, 'file') == 2, ...
    'Missing %s. Run Line_Stimuli_V4.m first.', tallPath);
assert(exist(resp3binPath, 'file') == 2, ...
    'Missing %s. Create the 3-bin response summary first.', resp3binPath);

Sgeo = load(tallPath);
assert(isfield(Sgeo, 'Tall_V4') && isstruct(Sgeo.Tall_V4), ...
    '%s must contain struct Tall_V4.', tallFile);
assert(isfield(Sgeo, 'RFrange') && ~isempty(Sgeo.RFrange), ...
    '%s must contain non-empty RFrange.', tallFile);

Tall_V4 = Sgeo.Tall_V4;
RFrange = Sgeo.RFrange(:);
nV4 = numel(RFrange);
siteRows = (1:nV4).';

R3_full = load_capsules_struct_exclusion_aware(resp3binPath, monkeySuffix, 'cfg', cfg);
R3dec = R3_full;
R3dec.meanAct = R3_full.meanAct(RFrange, :, :);
R3dec.meanSqAct = R3_full.meanSqAct(RFrange, :, :);
if ismatrix(R3_full.nTrials) && size(R3_full.nTrials,1) >= max(RFrange)
    R3dec.nTrials = R3_full.nTrials(RFrange, :);
else
    R3dec.nTrials = R3_full.nTrials;
end

SNRnorm = compute_snr_per_color_sites(R3dec, Tall_V4, siteRows, 'Verbose', false);
optsTD = struct('v1Sites', siteRows, 'timeIdx', timeIdx3, ...
    'excludeOverlap', excludeOverlap, 'verbose', false, 'epsDen', epsDen);
OUT = attention_modulation_V1_3bin(R3dec, Tall_V4, SNRnorm, optsTD);

% Verify that OUT.index uses the requested attention-index definition.
denCheck = abs((OUT.muT + OUT.muD) / 2);
denCheck(denCheck < epsDen) = epsDen;
idxCheck = (OUT.muT - OUT.muD) ./ denCheck;
okCheck = isfinite(OUT.index) & isfinite(idxCheck);
if any(okCheck)
    maxErr = max(abs(OUT.index(okCheck) - idxCheck(okCheck)));
    assert(maxErr < 1e-10, ...
        'Attention index mismatch: max abs difference %.3g.', maxErr);
end

idxAll = OUT.index;
isSig = isfinite(OUT.pValueTD) & (OUT.pValueTD < pThresh);

vAllRaw = idxAll(isfinite(idxAll));
vSigRaw = idxAll(isSig & isfinite(idxAll));
vAll = max(min(vAllRaw, idxClip), -idxClip);
vSig = max(min(vSigRaw, idxClip), -idxClip);

fprintf(['V4 attention histogram (%s): valid=%d / %d, significant=%d, ' ...
         'display clipped at +/-%.1f\n'], ...
    char(monkeySuffix), numel(vAll), nV4, numel(vSig), idxClip);
fprintf('Clipped values for display only: all=%d, significant=%d\n', ...
    nnz(abs(vAllRaw) > idxClip), nnz(abs(vSigRaw) > idxClip));

figure('Color', 'w');
hold on;
histogram(vAll, nHistBins, ...
    'FaceColor', [0.80 0.80 0.80], ...
    'EdgeColor', 'none');
histogram(vSig, nHistBins, ...
    'FaceColor', [0.85 0.20 0.20], ...
    'EdgeColor', 'none');
xlabel('Attention index');
ylabel('Number of sites');
title(sprintf('V4 Attention index (%s, all vs significant, pTD < %.3f)', ...
    char(monkeySuffix), pThresh));
legend(sprintf('All sites (N=%d)', numel(vAll)), ...
       sprintf('Significant sites (N=%d)', numel(vSig)), ...
       'Location', 'best');
grid on;

AttentionV4Hist = struct();
AttentionV4Hist.Monkey = Monkey;
AttentionV4Hist.monkeySuffix = monkeySuffix;
AttentionV4Hist.ATT = OUT;
AttentionV4Hist.idxAllRaw = vAllRaw;
AttentionV4Hist.idxSigRaw = vSigRaw;
AttentionV4Hist.idxClip = idxClip;
AttentionV4Hist.pThresh = pThresh;
