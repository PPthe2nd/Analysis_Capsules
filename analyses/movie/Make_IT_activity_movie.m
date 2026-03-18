% Render a visually driven IT activity movie for the line-task stimulus set.
%
% Uses the same projection and coloring path as the V1/V4 movie wrapper, but
% localizes the response and normalization structs to IT rows and uses a
% quartet-based responsiveness mask by default.

Monkey = 1;              % 1 = Nilson, 2 = Figaro
ExampleStimulus = 38;    % stimulus that defines the canonical projection frame for the movie
MinObjectStim = 1;       % site must fall on target/distractor in at least this many stimuli
RespThr = 0.7;           % quartet-responsiveness threshold in spontaneous-SD units
TopKQuartets = 5;        % responsiveness score uses mean abs(response) of the top-k quartets
UseQuartetResponsive = true; % true: use quartet-responsive mask, false: use all object-assigned IT sites
OnlyOnObjects = false;   % true: render only target/distractor points, false: include projected background RFs too
MarkerSize = 5;          % scatter marker size in the movie frames
ForceRender = false;     % if false, skip rendering when the output movie already exists

cfg = config();

if Monkey == 1
    monkeySuffix = "N";
    tallFile = 'Tall_IT_lines_N.mat';
    resp3binFile = 'SNR_capsules_N_d12.mat';
    respMovieFile = 'Resp_capsules_N_d12.mat';
elseif Monkey == 2
    monkeySuffix = "F";
    tallFile = 'Tall_IT_lines_F.mat';
    resp3binFile = 'SNR_capsules_F_d12.mat';
    respMovieFile = 'Resp_capsules_F_d12.mat';
else
    error('Make_IT_activity_movie:InvalidMonkey', ...
        'Monkey must be 1 (Nilson) or 2 (Figaro).');
end

if UseQuartetResponsive
    selectionSuffix = sprintf('top%d_%.2f', TopKQuartets, RespThr);
else
    selectionSuffix = 'allobj';
end
selectionSuffix = strrep(selectionSuffix, '.', 'p');

outFile = fullfile(cfg.resultsDir, ...
    sprintf('IT_activity_movie_%s_stim%03d_%s.mp4', ...
    char(monkeySuffix), ExampleStimulus, selectionSuffix));

if exist(outFile, 'file') == 2 && ~ForceRender
    fprintf('Skipping IT movie render because output already exists:\n%s\n', outFile);
    fprintf('Set ForceRender = true to re-render the movie.\n');
    return;
end

tallPath = fullfile(cfg.matDir, tallFile);
resp3binPath = fullfile(cfg.matDir, resp3binFile);
respMoviePath = fullfile(cfg.matDir, respMovieFile);

assert(exist(tallPath, 'file') == 2, ...
    'Missing %s. Run Line_Stimuli_IT.m first.', tallPath);
assert(exist(resp3binPath, 'file') == 2, ...
    'Missing %s. Create the 3-bin response summary first.', resp3binPath);
assert(exist(respMoviePath, 'file') == 2, ...
    'Missing %s. Create the high-resolution response summary first.', respMoviePath);

Sgeo = load(tallPath);
assert(isfield(Sgeo, 'Tall_IT') && isstruct(Sgeo.Tall_IT), ...
    '%s must contain struct Tall_IT.', tallFile);
assert(isfield(Sgeo, 'ALLCOORDS') && isfield(Sgeo, 'RTAB384') && isfield(Sgeo, 'RFrange'), ...
    '%s must contain Tall_IT, ALLCOORDS, RTAB384, and RFrange.', tallFile);

Tall_IT = Sgeo.Tall_IT;
ALLCOORDS = Sgeo.ALLCOORDS;
RTAB384 = Sgeo.RTAB384;
RFrange = Sgeo.RFrange(:);

nIT = numel(RFrange);
siteRows = (1:nIT).';

stimNums = arrayfun(@(x) x.stimNum, Tall_IT(:));
[stimNumsSorted, ordStim] = sort(stimNums(:));
assert(all(stimNumsSorted(:).' == 1:numel(Tall_IT)), ...
    'Tall_IT.stimNum must cover 1..%d exactly.', numel(Tall_IT));
Tall_IT = Tall_IT(ordStim);
nStim = numel(Tall_IT);

% Geometry summary needed for site-selection logic.
nTargetStim = zeros(nIT,1);
nDistrStim  = zeros(nIT,1);
for stim = 1:nStim
    T = Tall_IT(stim).T;
    assign = string(T.assignment(siteRows));
    nTargetStim = nTargetStim + (assign == "target");
    nDistrStim  = nDistrStim  + (assign == "distractor");
end
nObjectStim = nTargetStim + nDistrStim;
hasObjectRF = nObjectStim >= MinObjectStim;

% 3-bin responses for site screening and normalization.
S3 = load(resp3binPath);
assert(isfield(S3, 'R') && isstruct(S3.R), ...
    '%s must contain struct R.', resp3binFile);
R3_full = S3.R;
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

% Build quartet-based responsiveness scores.
nQuartets = floor(nStim / 8) * 2;
quartetMembers = zeros(nQuartets, 4);
q = 0;
for base = 0:8:(nStim - 8)
    q = q + 1;
    quartetMembers(q,:) = base + [1 2 5 6];
    q = q + 1;
    quartetMembers(q,:) = base + [3 4 7 8];
end

if isvector(R3.nTrials)
    nTrialsAll = double(R3.nTrials(:));
    perSiteTrials = false;
else
    perSiteTrials = true;
end

muQuartetEarly = nan(nIT, nQuartets);
muQuartetLate  = nan(nIT, nQuartets);

for iSite = 1:nIT
    if perSiteTrials
        nTrSite = double(R3.nTrials(iSite, :)).';
    else
        nTrSite = nTrialsAll;
    end

    rEarly = squeeze(R3.meanAct(iSite,:,2)).';
    rLate  = squeeze(R3.meanAct(iSite,:,3)).';

    for qIdx = 1:nQuartets
        stimQ = quartetMembers(qIdx,:);
        nTrQ = nTrSite(stimQ);

        rEarlyQ = rEarly(stimQ);
        idxEarly = isfinite(rEarlyQ) & isfinite(nTrQ) & (nTrQ > 0);
        if any(idxEarly)
            muQuartetEarly(iSite, qIdx) = sum(nTrQ(idxEarly) .* rEarlyQ(idxEarly)) / sum(nTrQ(idxEarly));
        end

        rLateQ = rLate(stimQ);
        idxLate = isfinite(rLateQ) & isfinite(nTrQ) & (nTrQ > 0);
        if any(idxLate)
            muQuartetLate(iSite, qIdx) = sum(nTrQ(idxLate) .* rLateQ(idxLate)) / sum(nTrQ(idxLate));
        end
    end
end

signedQuartetEarly = bsxfun(@rdivide, bsxfun(@minus, muQuartetEarly, muSpont), sdSpont);
signedQuartetLate  = bsxfun(@rdivide, bsxfun(@minus, muQuartetLate,  muSpont), sdSpont);
badNoise = ~isfinite(sdSpont) | (sdSpont <= 0);
signedQuartetEarly(badNoise, :) = NaN;
signedQuartetLate(badNoise, :) = NaN;

topKAbsQuartetEarly = nan(nIT,1);
topKAbsQuartetLate = nan(nIT,1);
for iSite = 1:nIT
    valsEarly = abs(signedQuartetEarly(iSite, :));
    valsEarly = valsEarly(isfinite(valsEarly));
    if ~isempty(valsEarly)
        valsEarly = sort(valsEarly, 'descend');
        kEarly = min(TopKQuartets, numel(valsEarly));
        topKAbsQuartetEarly(iSite) = mean(valsEarly(1:kEarly));
    end

    valsLate = abs(signedQuartetLate(iSite, :));
    valsLate = valsLate(isfinite(valsLate));
    if ~isempty(valsLate)
        valsLate = sort(valsLate, 'descend');
        kLate = min(TopKQuartets, numel(valsLate));
        topKAbsQuartetLate(iSite) = mean(valsLate(1:kLate));
    end
end

quartetRespScore = max(topKAbsQuartetEarly, topKAbsQuartetLate);
isResponsiveQuartet = hasObjectRF & isfinite(quartetRespScore) & (quartetRespScore > RespThr);

if UseQuartetResponsive
    siteIdx = find(isResponsiveQuartet);
    selectionLabel = sprintf('top-%d quartet-responsive (>%0.2f)', TopKQuartets, RespThr);
else
    siteIdx = find(hasObjectRF);
    selectionLabel = sprintf('all object-assigned (>= %d object stim)', MinObjectStim);
end

assert(~isempty(siteIdx), ...
    'No IT sites passed the requested selection (%s).', selectionLabel);

fprintf('Rendering IT movie for monkey %s\n', char(monkeySuffix));
fprintf('Example stimulus: %d\n', ExampleStimulus);
fprintf('Selection: %s (%d sites)\n', selectionLabel, numel(siteIdx));
fprintf('OnlyOnObjects: %d\n', OnlyOnObjects);
fprintf('Output: %s\n', outFile);

% High-resolution response struct for movie frames.
Sresp = load(respMoviePath);
assert(isfield(Sresp, 'R') && isstruct(Sresp.R), ...
    '%s must contain struct R.', respMovieFile);
R_resp_full = Sresp.R;
R_resp = R_resp_full;
R_resp.meanAct = R_resp_full.meanAct(RFrange, :, :);
R_resp.meanSqAct = R_resp_full.meanSqAct(RFrange, :, :);
if ismatrix(R_resp_full.nTrials) && size(R_resp_full.nTrials,1) >= max(RFrange)
    R_resp.nTrials = R_resp_full.nTrials(RFrange, :);
else
    R_resp.nTrials = R_resp_full.nTrials;
end

make_activity_movie_wrapper_safe(outFile, ...
    Tall_IT, ALLCOORDS, RTAB384, ExampleStimulus, R_resp, SNR, ...
    'UseOnlyV1', false, ...
    'OnlyOnObjects', OnlyOnObjects, ...
    'SiteIdx', siteIdx, ...
    'MarkerSize', MarkerSize);
