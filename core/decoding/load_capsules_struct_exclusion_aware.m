function R = load_capsules_struct_exclusion_aware(pooledPath, monkeySuffix, varargin)
% LOAD_CAPSULES_STRUCT_EXCLUSION_AWARE
% Load a pooled response summary (e.g. SNR_capsules_*_d12.mat or
% Resp_capsules_*_d12.mat) and, if site-by-session exclusions exist,
% recombine day-specific summaries so that bad sessions are dropped while
% preserving the good session for each site.
%
% The returned struct R matches the usual fields, but nTrials is returned
% as an [NChans x NStim] matrix because the effective per-stimulus trial
% count becomes site-specific once session exclusions are applied.

p = inputParser;
p.addRequired('pooledPath', @(x) ischar(x) || isstring(x));
p.addRequired('monkeySuffix', @(x) ischar(x) || isstring(x));
p.addParameter('cfg', [], @(x) isempty(x) || isstruct(x));
p.addParameter('sessions', [1 2], @(x) isnumeric(x) && isvector(x) && numel(x) >= 1);
p.addParameter('onlyCorrect', true, @(x)islogical(x) && isscalar(x));
p.addParameter('correctCol', 9, @(x)isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('correctVal', 1, @(x)isnumeric(x) && isscalar(x));
p.addParameter('dayCol', 11, @(x)isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('chunkTrials', 200, @(x)isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('expectedMaxStim', 384, @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>=1));
p.addParameter('useCache', true, @(x)islogical(x) && isscalar(x));
p.parse(pooledPath, monkeySuffix, varargin{:});
opt = p.Results;

pooledPath = char(pooledPath);
monkeySuffix = string(monkeySuffix);

if isempty(opt.cfg)
    cfg = config();
else
    cfg = opt.cfg;
end

assert(exist(pooledPath, 'file') == 2, 'Missing %s.', pooledPath);
S = load(pooledPath, 'R');
assert(isfield(S, 'R') && isstruct(S.R), '%s must contain struct R.', pooledPath);
Rpooled = S.R;

Texcl = site_session_exclusions(monkeySuffix);
if isempty(Texcl)
    R = Rpooled;
    return;
end

[~, daySplitCachePath] = local_cache_info(cfg, pooledPath, monkeySuffix);
[Split, ~, ~] = load_capsules_day_split(pooledPath, monkeySuffix, ...
    'cfg', cfg, ...
    'sessions', opt.sessions, ...
    'onlyCorrect', opt.onlyCorrect, ...
    'correctCol', opt.correctCol, ...
    'correctVal', opt.correctVal, ...
    'dayCol', opt.dayCol, ...
    'chunkTrials', opt.chunkTrials, ...
    'expectedMaxStim', opt.expectedMaxStim, ...
    'useCache', opt.useCache);

R = local_combine_splits_with_exclusions(Rpooled, Split, Texcl, opt.sessions);
R.sessionExclusions = Texcl;
R.sessionIdsUsed = opt.sessions(:);
R.sourcePooledPath = pooledPath;
R.sourceDaySplitCache = daySplitCachePath;
end

function [kindName, cachePath] = local_cache_info(cfg, pooledPath, monkeySuffix)
[~, baseName, ~] = fileparts(pooledPath);
if contains(baseName, 'SNR_capsules')
    kindName = 'SNR';
elseif contains(baseName, 'Resp_capsules')
    kindName = 'Resp';
else
    error('Unsupported capsules file: %s', pooledPath);
end
cacheName = sprintf('%s_daySplit_%s.mat', baseName, char(monkeySuffix));
cachePath = fullfile(cfg.resultsDir, cacheName);
end

function R = local_combine_splits_with_exclusions(Rpooled, Split, Texcl, sessions)
nCh = size(Rpooled.meanAct, 1);
nStim = size(Rpooled.meanAct, 2);
nWin = size(Rpooled.meanAct, 3);
nSess = numel(sessions);
assert(numel(Split) == nSess, 'Split/session mismatch.');

keepBySess = true(nCh, nSess);
for i = 1:height(Texcl)
    sessIdx = find(sessions == Texcl.day(i), 1, 'first');
    if isempty(sessIdx)
        continue;
    end
    ch = double(Texcl.siteGlobal(i));
    if ch >= 1 && ch <= nCh
        keepBySess(ch, sessIdx) = false;
    end
end

sumAct = zeros(nCh, nStim, nWin, 'double');
sumSqAct = zeros(nCh, nStim, nWin, 'double');
nTrialsMat = zeros(nCh, nStim, 'double');

for iSess = 1:nSess
    meanAct = double(Split(iSess).meanAct);
    meanSqAct = double(Split(iSess).meanSqAct);
    nTrials = double(Split(iSess).nTrials(:)');
    assert(size(meanAct,1) == nCh && size(meanAct,2) == nStim && size(meanAct,3) == nWin, ...
        'Day split meanAct size mismatch.');
    keep = keepBySess(:, iSess);
    nSiteStim = repmat(nTrials, nCh, 1);
    nSiteStim(~keep, :) = 0;

    sumActSess = bsxfun(@times, meanAct, reshape(nSiteStim, [nCh nStim 1]));
    sumSqSess = bsxfun(@times, meanSqAct, reshape(nSiteStim, [nCh nStim 1]));
    sumActSess(~isfinite(sumActSess)) = 0;
    sumSqSess(~isfinite(sumSqSess)) = 0;

    sumAct = sumAct + sumActSess;
    sumSqAct = sumSqAct + sumSqSess;
    nTrialsMat = nTrialsMat + nSiteStim;
end

R = Rpooled;
R.meanAct = sumAct ./ reshape(max(nTrialsMat, 1), [nCh nStim 1]);
R.meanSqAct = sumSqAct ./ reshape(max(nTrialsMat, 1), [nCh nStim 1]);
zeroMask3 = repmat(nTrialsMat == 0, [1 1 nWin]);
R.meanAct(zeroMask3) = NaN;
R.meanSqAct(zeroMask3) = NaN;
R.nTrials = nTrialsMat;
R.siteSessionKeepMask = keepBySess;
end
