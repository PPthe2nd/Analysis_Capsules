% Render a visually driven V4 activity movie for the line-task stimulus set.
%
% Uses the same projection and coloring path as the V1 movie wrapper, but
% localizes the response and normalization structs to V4 rows and reuses the
% current first-pass V4 selection logic.

Monkey = 1;              % 1 = Nilson, 2 = Figaro
ExampleStimulus = 38;    % stimulus that defines the canonical projection frame for the movie
SNRthr = 0.7;            % minimum onset SNR for the primary visually driven site set
pTDthr = 0.05;           % p-value threshold for attention-modulation rescue
NminMatched = 20;        % minimum matched T/D trial weight after color balancing for attention rescue
MinObjectStim = 1;       % site must fall on target/distractor in at least this many stimuli
UseExpandedSelection = true; % true: onset OR attention sites, false: onset-only sites
OnlyOnObjects = false;   % true: render only target/distractor points, false: include projected background RFs too
MarkerSize = 5;          % scatter marker size in the movie frames
ForceRender = false;     % if false, skip rendering when the output movie already exists

cfg = config();

if Monkey == 1
    monkeySuffix = "N";
    tallFile = 'Tall_V4_lines_N.mat';
    resp3binFile = 'SNR_capsules_N_d12.mat';
    respMovieFile = 'Resp_capsules_N_d12.mat';
elseif Monkey == 2
    monkeySuffix = "F";
    tallFile = 'Tall_V4_lines_F.mat';
    resp3binFile = 'SNR_capsules_F_d12.mat';
    respMovieFile = 'Resp_capsules_F_d12.mat';
else
    error('Make_V4_activity_movie:InvalidMonkey', ...
        'Monkey must be 1 (Nilson) or 2 (Figaro).');
end

outFile = fullfile(cfg.resultsDir, ...
    sprintf('V4_activity_movie_%s_stim%03d.mp4', char(monkeySuffix), ExampleStimulus));

if exist(outFile, 'file') == 2 && ~ForceRender
    fprintf('Skipping V4 movie render because output already exists:\n%s\n', outFile);
    fprintf('Set ForceRender = true to re-render the movie.\n');
    return;
end

tallPath = fullfile(cfg.matDir, tallFile);
resp3binPath = fullfile(cfg.matDir, resp3binFile);
respMoviePath = fullfile(cfg.matDir, respMovieFile);

assert(exist(tallPath, 'file') == 2, ...
    'Missing %s. Run Line_Stimuli_V4.m first.', tallPath);
assert(exist(resp3binPath, 'file') == 2, ...
    'Missing %s. Create the 3-bin response summary first.', resp3binPath);
assert(exist(respMoviePath, 'file') == 2, ...
    'Missing %s. Create the high-resolution response summary first.', respMoviePath);

Sgeo = load(tallPath);
assert(isfield(Sgeo, 'Tall_V4') && isstruct(Sgeo.Tall_V4), ...
    '%s must contain struct Tall_V4.', tallFile);
assert(isfield(Sgeo, 'ALLCOORDS') && isfield(Sgeo, 'RTAB384') && isfield(Sgeo, 'RFrange'), ...
    '%s must contain Tall_V4, ALLCOORDS, RTAB384, and RFrange.', tallFile);

Tall_V4 = Sgeo.Tall_V4;
ALLCOORDS = Sgeo.ALLCOORDS;
RTAB384 = Sgeo.RTAB384;
RFrange = Sgeo.RFrange(:);

nV4 = numel(RFrange);
siteRows = (1:nV4).';

% Geometry summary needed for the same inclusion logic used in Analyse_Line_Stimuli_V4.
nTargetStim = zeros(nV4,1);
nDistrStim  = zeros(nV4,1);

for stim = 1:numel(Tall_V4)
    T = Tall_V4(stim).T;
    assign = string(T.assignment(siteRows));
    nTargetStim = nTargetStim + (assign == "target");
    nDistrStim  = nDistrStim  + (assign == "distractor");
end

nObjectStim = nTargetStim + nDistrStim;
hasObjectRF = nObjectStim >= MinObjectStim;

% 3-bin responses for site screening and normalization.
R3_full = load_capsules_struct_exclusion_aware(resp3binPath, monkeySuffix, 'cfg', cfg);
R3 = R3_full;
R3.meanAct = R3_full.meanAct(RFrange, :, :);
R3.meanSqAct = R3_full.meanSqAct(RFrange, :, :);
if ismatrix(R3_full.nTrials) && size(R3_full.nTrials,1) >= max(RFrange)
    R3.nTrials = R3_full.nTrials(RFrange, :);
else
    R3.nTrials = R3_full.nTrials;
end

SNR = compute_snr_per_color_sites(R3, Tall_V4, siteRows, 'Verbose', false);
SNRmat = [SNR.yellowEarly(siteRows), SNR.yellowLate(siteRows), ...
          SNR.purpleEarly(siteRows), SNR.purpleLate(siteRows)];
[bestSNR, ~] = max(SNRmat, [], 2, 'omitnan');

optsTD = struct('v1Sites', siteRows, 'timeIdx', 3, 'excludeOverlap', true, 'verbose', false);
OUTtd = attention_modulation_V1_3bin(R3, Tall_V4, SNR, optsTD);

matchedN = OUTtd.wY + OUTtd.wP;
isOnsetDriven = hasObjectRF & isfinite(bestSNR) & (bestSNR > SNRthr);
isAttentionMod = hasObjectRF & isfinite(OUTtd.pValueTD) & ...
    (OUTtd.pValueTD < pTDthr) & (matchedN >= NminMatched);

keepSitesPrimary = find(isOnsetDriven);
keepSitesExpanded = find(isOnsetDriven | isAttentionMod);

if UseExpandedSelection
    siteIdx = keepSitesExpanded;
    selectionLabel = 'onset OR attention';
else
    siteIdx = keepSitesPrimary;
    selectionLabel = 'onset only';
end

assert(~isempty(siteIdx), ...
    'No V4 sites passed the requested selection (%s).', selectionLabel);

fprintf('Rendering V4 movie for monkey %s\n', char(monkeySuffix));
fprintf('Example stimulus: %d\n', ExampleStimulus);
fprintf('Selection: %s (%d sites)\n', selectionLabel, numel(siteIdx));
fprintf('OnlyOnObjects: %d\n', OnlyOnObjects);
fprintf('Output: %s\n', outFile);

% High-resolution response struct for movie frames.
R_resp_full = load_capsules_struct_exclusion_aware(respMoviePath, monkeySuffix, 'cfg', cfg);
R_resp = R_resp_full;
R_resp.meanAct = R_resp_full.meanAct(RFrange, :, :);
R_resp.meanSqAct = R_resp_full.meanSqAct(RFrange, :, :);
if ismatrix(R_resp_full.nTrials) && size(R_resp_full.nTrials,1) >= max(RFrange)
    R_resp.nTrials = R_resp_full.nTrials(RFrange, :);
else
    R_resp.nTrials = R_resp_full.nTrials;
end

make_activity_movie_wrapper_safe(outFile, ...
    Tall_V4, ALLCOORDS, RTAB384, ExampleStimulus, R_resp, SNR, ...
    'UseOnlyV1', false, ...
    'OnlyOnObjects', OnlyOnObjects, ...
    'SiteIdx', siteIdx, ...
    'MarkerSize', MarkerSize);
