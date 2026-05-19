function [Split, Rpooled, cachePath] = load_capsules_day_split(pooledPath, monkeySuffix, varargin)
% LOAD_CAPSULES_DAY_SPLIT
% Load or build day-specific stimulus averages corresponding to a pooled
% capsules response file (e.g. SNR_capsules_*_d12.mat or Resp_capsules_*_d12.mat).
%
% Returns:
%   Split    : struct array with fields sessionId, meanAct, meanSqAct,
%              nTrials, stimList
%   Rpooled  : pooled struct R loaded from pooledPath
%   cachePath: path of the cached day-split file in resultsDir

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

[~, baseName, ~] = fileparts(pooledPath);
cacheName = sprintf('%s_daySplit_%s.mat', baseName, char(monkeySuffix));
cachePath = fullfile(cfg.resultsDir, cacheName);
wantSessions = opt.sessions(:).';

if opt.useCache && exist(cachePath, 'file') == 2
    S = load(cachePath, 'Split', 'meta');
    if isfield(S, 'Split') && isfield(S, 'meta') ...
            && isequal(double(S.meta.timeWindows), double(Rpooled.timeWindows)) ...
            && isequal(double(S.meta.sessions(:).'), double(wantSessions)) ...
            && isequal(double(S.meta.tb(:)), double(Rpooled.tb(:)))
        Split = S.Split;
        return;
    end
end

if monkeySuffix == "N"
    monkeyFolder = 'Mr Nilson';
elseif monkeySuffix == "F"
    monkeyFolder = 'Figaro';
else
    error('Unknown monkey suffix %s.', char(monkeySuffix));
end

trialRespPath = fullfile(cfg.dataRoot, monkeyFolder, 'ObjAtt_lines_normMUA.mat');
trialInfoPath = fullfile(cfg.dataRoot, monkeyFolder, 'ObjAtt_lines_MUA_trials.mat');
assert(exist(trialRespPath, 'file') == 2, 'Missing %s.', trialRespPath);
assert(exist(trialInfoPath, 'file') == 2, 'Missing %s.', trialInfoPath);

m1 = matfile(trialRespPath);
m2 = matfile(trialInfoPath);
Split = repmat(struct('sessionId', [], 'meanAct', [], 'meanSqAct', [], 'nTrials', [], 'stimList', []), ...
    numel(wantSessions), 1);

for iSess = 1:numel(wantSessions)
    sess = wantSessions(iSess);
    fprintf('Building day split for %s, monkey %s, day %d...\n', baseName, char(monkeySuffix), sess);
    [meanAct, meanSqAct, nTrials, stimList] = avg_byStim( ...
        m1, m2, double(Rpooled.timeWindows), ...
        'onlyCorrect', opt.onlyCorrect, ...
        'correctCol', opt.correctCol, ...
        'correctVal', opt.correctVal, ...
        'days', sess, ...
        'dayCol', opt.dayCol, ...
        'chunkTrials', opt.chunkTrials, ...
        'expectedMaxStim', opt.expectedMaxStim, ...
        'verbose', true);
    Split(iSess).sessionId = sess;
    Split(iSess).meanAct = meanAct;
    Split(iSess).meanSqAct = meanSqAct;
    Split(iSess).nTrials = nTrials;
    Split(iSess).stimList = stimList;
end

if opt.useCache
    meta = struct();
    meta.timeWindows = double(Rpooled.timeWindows);
    meta.tb = double(Rpooled.tb(:));
    meta.sessions = wantSessions;
    save(cachePath, 'Split', 'meta', '-v7.3');
    fprintf('Saved day-split cache to %s\n', cachePath);
end
end
