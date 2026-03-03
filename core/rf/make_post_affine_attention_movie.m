function M = make_post_affine_attention_movie(outMovie, stimID_example, ALLCOORDS, RTAB384, OUT_postAffine, varargin)
% MAKE_POST_AFFINE_ATTENTION_MOVIE
% Build movie frames from post-affine bins using stream-aware KNN smoothing
% and prep-calibrated thresholds.

p = inputParser;
p.addParameter('K', 30, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('siteWeightsByIndex', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.addParameter('alphaFullAt', 0.2, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('colorRedAt', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
p.addParameter('cMaxFixed', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
p.addParameter('markerSize', 8, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('alpha', 1, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
p.addParameter('bgColor', [0.5 0.5 0.5], @(x) isnumeric(x) && numel(x)==3);
p.addParameter('cLow', [0.50 0.50 0.50], @(x) isnumeric(x) && numel(x)==3);
p.addParameter('cHigh', [0.85 0.05 0.05], @(x) isnumeric(x) && numel(x)==3);
p.addParameter('hotScale', false, @(x) islogical(x) && isscalar(x));
p.addParameter('colorHotMaxFactor', 8.0, @(x) isnumeric(x) && isscalar(x) && x > 1);
p.addParameter('hardAlphaCutoff', false, @(x) islogical(x) && isscalar(x));
p.addParameter('frameRate', 10, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('quality', 95, @(x) isnumeric(x) && isscalar(x) && x >= 1 && x <= 100);
p.addParameter('enforceK', false, @(x) islogical(x) && isscalar(x));
p.addParameter('verbose', true, @(x) islogical(x) && isscalar(x));
p.addParameter('timeLabelRef', 'start', @(x) ischar(x) || isstring(x));
p.addParameter('stimOnsetMs', 0, @(x) isnumeric(x) && isscalar(x) && isfinite(x));
p.parse(varargin{:});
opt = p.Results;

assert(isstruct(OUT_postAffine) && isfield(OUT_postAffine,'bins') && ~isempty(OUT_postAffine.bins), ...
    'OUT_postAffine must contain non-empty bins.');
assert(isfield(OUT_postAffine.bins, 'stream'), ...
    'OUT_postAffine.bins.stream missing. Rebuild post-affine export first.');

bins = OUT_postAffine.bins;
nFrames = numel(bins);

if isfield(OUT_postAffine, 'meta') && isfield(OUT_postAffine.meta, 'timeWindows')
    twAll = double(OUT_postAffine.meta.timeWindows);
else
    twAll = nan(nFrames,2);
    for tb = 1:nFrames
        twAll(tb,:) = double(bins(tb).timeWindow);
    end
end

if isempty(opt.colorRedAt)
    redAtUse = opt.alphaFullAt;
else
    redAtUse = opt.colorRedAt;
end

if opt.verbose
    fprintf('Writing post-affine movie: %s\n', outMovie);
    fprintf(['  nFrames=%d | K=%d | alphaFullAt=%.6g | colorRedAt=%.6g | ' ...
             'cMaxFixed=%s | stimOnsetMs=%.3g | hardAlphaCutoff=%d\n'], ...
        nFrames, round(opt.K), opt.alphaFullAt, redAtUse, mat2str(opt.cMaxFixed), ...
        opt.stimOnsetMs, opt.hardAlphaCutoff);
end

vw = VideoWriter(outMovie, 'MPEG-4');
vw.FrameRate = opt.frameRate;
vw.Quality = round(opt.quality);
open(vw);

fracAbove = nan(nFrames,1);
for tb = 1:nFrames
    if size(twAll,1) >= tb
        tw = twAll(tb,:);
    else
        tw = [NaN NaN];
    end

    showStimulus = true;
    if all(isfinite(tw))
        % Draw stimulus from the first bin whose window starts at/after onset.
        showStimulus = (tw(1) >= opt.stimOnsetMs);
    end

    h = plot_post_affine_knn_frame( ...
        stimID_example, ALLCOORDS, RTAB384, bins(tb), ...
        'K', round(opt.K), ...
        'siteWeightsByIndex', opt.siteWeightsByIndex, ...
        'alphaFullAt', opt.alphaFullAt, ...
        'colorRedAt', redAtUse, ...
        'cMaxFixed', opt.cMaxFixed, ...
        'markerSize', opt.markerSize, ...
        'alpha', opt.alpha, ...
        'bgColor', opt.bgColor, ...
        'cLow', opt.cLow, ...
        'cHigh', opt.cHigh, ...
        'hotScale', opt.hotScale, ...
        'colorHotMaxFactor', opt.colorHotMaxFactor, ...
        'hardAlphaCutoff', opt.hardAlphaCutoff, ...
        'timeWindow', tw, ...
        'timeLabelRef', opt.timeLabelRef, ...
        'showStimulus', showStimulus, ...
        'enforceK', opt.enforceK);

    fracAbove(tb) = h.fracAboveThreshold;
    fr = getframe(h.fig);
    writeVideo(vw, fr);
    close(h.fig);

    if opt.verbose && (tb == 1 || tb == nFrames || mod(tb,10) == 0)
        fprintf('  frame %2d/%2d | >thr %.2f%%\n', tb, nFrames, 100*fracAbove(tb));
    end
end

close(vw);

M = struct();
M.outMovie = outMovie;
M.nFrames = nFrames;
M.K = round(opt.K);
M.alphaFullAt = opt.alphaFullAt;
M.colorRedAt = redAtUse;
M.cMaxFixed = opt.cMaxFixed;
M.stimOnsetMs = opt.stimOnsetMs;
M.hardAlphaCutoff = opt.hardAlphaCutoff;
M.fracAboveByFrame = fracAbove;
M.meanFracAbove = mean(fracAbove, 'omitnan');

if opt.verbose
    fprintf('Movie done: %s\n', outMovie);
    fprintf('Mean frame exceedance >thr: %.2f%%\n', 100*M.meanFracAbove);
end

end
