% remember that arc_frac_edge is 0 close to the distractor and 1 close to
% the target

Monkey = 1; % 1 for Nilson, 2 for Figaro
TabFile = "ObjAtt_lines_monkeyN_20220201_B1";
cfg = config();
LogDir = cfg.logsDir;
hasSessionExclusions = ~isempty(site_session_exclusions("N"));

load(fullfile(cfg.matDir, "Tall_V1_lines_N.mat")); % load the RF information for each of 384 stimuli

% analyze only the correct trials; if all trials need to be analyzed, set
% onlyCorrect to false

% timeWindows_SNR = [-200 0; 40 240; 300 500];
% [meanActSNR, meanSqActSNR, nTrials, stimList] = ...
%    avg_byStim(m1, m2, timeWindows_SNR, 'days', [1 2]);  % only days 1 and 2 here
% R.meanAct   = meanActSNR;R.meanSqAct = meanSqActSNR;R.nTrials   = nTrials;
% R.stimList  = stimList;R.timeWindows = timeWindows_SNR;R.tb = tb;
% R.file_m1 = m1;R.file_m2 = m2;
% save the result in a mat file:
% save('SNR_capsules_N_d12.mat', 'R', '-v7.3');

R_snr = load_capsules_struct_exclusion_aware(fullfile(cfg.matDir, 'SNR_capsules_N_d12.mat'), "N", 'cfg', cfg);

% The logic of day 1 and 2:
% Stim 1 purple is target; Stim 6 same stimulus, but yellow is target; Stim 2 purple is target; Stim 5 same stimulus but yellow is target
% 1 and 5 have the target curve at the same location, so do 2 and 6

% Stim 3 purple is target; Stim 8 same stimulus, but yellow is target; Stim 4 purple is target; Stim 7 same stimulus but purple is target
% 3 and 7 have the target curve at the same location, so do 4 and 8

% ComputeSNR_perColor
% computes the SNR of all channels in every window - the result is in SNR_V1_byColor_byWindow.mat



%% --> ComputeSNR_perColor; % determines the SNR in the two time windows.
if hasSessionExclusions
    fprintf('Session exclusions are active for monkey N; recomputing V1 SNR from exclusion-aware 3-bin responses.\n');
    SNR = compute_snr_per_color_sites(R_snr, Tall_V1, (1:512).', 'Verbose', false);
else
    load(fullfile(cfg.matDir, 'SNR_V1_byColor_byWindow.mat'));
end
SNRmat = [SNR.yellowEarly(1:512), SNR.yellowLate(1:512), ...
    SNR.purpleEarly(1:512), SNR.purpleLate(1:512)];
[bestSNR, bestIdx] = max(SNRmat, [], 2, 'omitnan');
bestSNR = bestSNR(:);

x = bestSNR(isfinite(bestSNR));
figure;
histogram(x, 30);   % 30 bins 
xlabel('Max SNR');
ylabel('N sites');
title(sprintf('Max SNR per site (N = %d)', numel(x)));
grid on;

% --> Compare_Paolos_SNR_Color_SNR; % this script can be run to compare these results
% to the SNR defined in m1.SNR (which has 1024 channels)

SNRthr = 0.7;
pTDthr = 0.05;          % TD significance threshold (attention rescue)
NminMatched = 20;       % minimum matched T/D trial weight for rescue
pColorThr = 0.05;       % color significance threshold (color rescue)

% TD modulation from 3-bin data (same as Attention_Line_Stimuli logic)
R3 = R_snr;
optsTD = struct('timeIdx', 3, 'excludeOverlap', true, 'verbose', false);
OUTtd = attention_modulation_V1_3bin(R3, Tall_V1, SNR, optsTD);

% Color significance (purple vs yellow) from ColorTune
if exist(fullfile(cfg.matDir, 'ColorTune_balanced_V1.mat'), 'file') == 2 && ~hasSessionExclusions
    Sct = load(fullfile(cfg.matDir, 'ColorTune_balanced_V1.mat'));  % loads ColorTune
    ColorTune = Sct.ColorTune;
else
    keepSiteIdx = find(bestSNR > SNRthr);
    ColorTune = compute_color_tuning_balanced_sites(R3, Tall_V1, (1:512).', keepSiteIdx, 'Verbose', false);
    ColorTune.thr = SNRthr;
    ColorTune.bestSNR = bestSNR;
end
isColorSig = isfinite(ColorTune.early.p(1:512)) & (ColorTune.early.p(1:512) < pColorThr);

matchedN = OUTtd.wY + OUTtd.wP;
isMain = bestSNR > SNRthr;
isRescue = isfinite(OUTtd.pValueTD) & (OUTtd.pValueTD < pTDthr) & (matchedN >= NminMatched);
isKeep = isMain | isRescue | isColorSig;
keepSites = find(isKeep);   % indices 1..512

fprintf('Main SNR-selected sites (bestSNR > %.2f): %d / 512\n', SNRthr, nnz(isMain));
fprintf('TD rescue-only sites (pTD < %.3f, matchedN >= %d): %d\n', ...
    pTDthr, NminMatched, nnz(isRescue & ~isMain));
fprintf('Color rescue-only sites (pColor < %.3f): %d\n', ...
    pColorThr, nnz(isColorSig & ~isMain & ~isRescue));
fprintf('Total kept sites (union): %d / 512\n', numel(keepSites));

ExampleStimulus = 38;
h = plot_projected_RFs_on_example_stim(Tall_V1, ALLCOORDS, RTAB384, ExampleStimulus, ...
    'MarkerSize', 4, 'Alpha', 0.15);

h = plot_projected_RFs_on_example_stim(Tall_V1, ALLCOORDS, RTAB384, ExampleStimulus, ...
    'MarkerSize', 4, 'Alpha', 0.15,'SiteIdx', keepSites);

% look at distribution of RFs along the curve
x = [];  % histogram of GC values, x will collect along_GC values
for stim = 1:numel(Tall_V1)
    T = Tall_V1(stim).T;          % 512x18 table
    v = T.along_GC(keepSites, :); % keep rows; keep all columns (if along_GC is 512x18 inside table entries)
    x = [x; v(:)];
end

x = x(isfinite(x));  % remove NaN/Inf

figure;
histogram(x, 'BinMethod','fd');
xlabel('along\_GC');
ylabel('Count');
title(sprintf('along\\_GC for sites with bestSNR > %.2f (N=%d values)', SNRthr, numel(x)));
grid on;

%% Color tuning
% --> ColorTuning_Capsules; % determines the color tuning in the two time windows.
if ~exist('ColorTune','var')
    if exist(fullfile(cfg.matDir, 'ColorTune_balanced_V1.mat'), 'file') == 2 && ~hasSessionExclusions
        load(fullfile(cfg.matDir, 'ColorTune_balanced_V1.mat'));  % SNR threshold was 0.7
    else
        keepSiteIdx = find(bestSNR > SNRthr);
        ColorTune = compute_color_tuning_balanced_sites(R3, Tall_V1, (1:512).', keepSiteIdx, 'Verbose', false);
        ColorTune.thr = SNRthr;
        ColorTune.bestSNR = bestSNR;
    end
end

% histogram with colour tuning
alpha = 0.05; idx = keepSites;
ci = ColorTune.early.colorIndex(idx); pvals = ColorTune.early.p(idx);
valid = isfinite(ci) & isfinite(pvals); ci = ci(valid);
pvals = pvals(valid); sig = pvals < alpha;
figure; hold on;
histogram(ci, 30, ...
    'FaceColor', [0.8 0.8 0.8], ...
    'EdgeColor', 'none');
histogram(ci(sig), 30, ...
    'FaceColor', [0.8 0 0], ...
    'EdgeColor', 'none');
xlabel('Color Index (yellow - purple)'); ylabel('Number of sites');
title(sprintf('V1 Color tuning (40-240ms), SNR>%.2f', SNRthr));
legend('All sites','Significant (p<0.05)'); grid on;

% compute the entire matrix of responses per site per stimulus, day 1 and 2.
winSize = 10;                 % ms
edges = -200:winSize:500;     % window edges
% timeWindowsResp = [edges(1:end-1).'  edges(2:end).'];
% [meanActResp, meanSqActResp, nTrials, stimList] = ...
%    avg_byStim(m1, m2, timeWindowsResp, 'days', [1 2]);  % only days 1 and 2 here
% 
% R.meanAct   = meanActResp;R.meanSqAct = meanSqActResp;R.nTrials   = nTrials;
% R.stimList  = stimList;R.timeWindows = timeWindowsResp;R.tb = tb;
% R.file_m1 = m1;R.file_m2 = m2;
% save('Resp_capsules_N_d12.mat', 'R', '-v7.3');
R_resp = load_capsules_struct_exclusion_aware(fullfile(cfg.matDir, 'Resp_capsules_N_d12.mat'), "N", 'cfg', cfg);

% best vs worst color and grey with a minimal distance to a coloured pixel
ciThr = 0.35;                 % abs(colorIndex) threshold for "color tuned"
useColorIndexFrom = "early";   % "early" or "late"
minDistThr = 30;              % px: include ONLY gray stimuli with distance >= minDistThr in pixels
R = R_resp; %#ok<NASGU> % legacy workspace name expected by PSTH_colorPref_V1_grayDist
PSTH_colorPref_V1_grayDist
% PSTH_colorPref_V1;            % best vs worst color response for tuned
% sites, no grey
% first attempt on colour decoding in 10ms time bins, based on PSTHs. This 
% is a coarse, first attempt at this and it should be done better at the single trial level  
Color_decoding_V1; 

% Plot activity overlay on example stimulus, time bin 30
h = plot_projected_activity_on_example_stim(Tall_V1, ALLCOORDS, RTAB384, ExampleStimulus, R_resp, SNR, ...
    'TimeBin', 4, 'UseOnlyV1', true, 'OnlyOnObjects', false, 'MarkerSize', 5, 'SiteIdx', keepSites);
h = plot_projected_activity_on_example_stim(Tall_V1, ALLCOORDS, RTAB384, ExampleStimulus, R_resp, SNR, ...
    'TimeBin', 24, 'UseOnlyV1', true, 'OnlyOnObjects', false, 'MarkerSize', 5, 'SiteIdx', keepSites);




% Plot activity overlay on example stimulus, time bin 30
% make_activity_movie(Tall, ALLCOORDS, RTAB384, 1, R, SNR, 'V1_activity_movie.mp4');

% make_activity_movie_wrapper_safe('V1_activity_movie_tmp2.mp4', ...
%     Tall_V1, ALLCOORDS, RTAB384, ExampleStimulus, R_resp, SNR, 'UseOnlyV1', true,'OnlyOnObjects', false, ...
%     'SiteIdx', keepSites, 'MarkerSize', 5, 'SiteIdx', keepSites);
