function M = make_grouped_alonggc_polygon_movie(outMovie, stimID_example, ALLCOORDS, RTAB384, G, varargin)
% MAKE_GROUPED_ALONGGC_POLYGON_MOVIE
% Render grouped along_GC polygon frames across all time bins.

p = inputParser;
p.addParameter('alphaFullAt', 0.2, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('colorRedAt', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
p.addParameter('cMaxFixed', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
p.addParameter('alpha', 1, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
p.addParameter('bgColor', [0.5 0.5 0.5], @(x) isnumeric(x) && numel(x)==3);
p.addParameter('cLow', [0.50 0.50 0.50], @(x) isnumeric(x) && numel(x)==3);
p.addParameter('cHigh', [0.85 0.05 0.05], @(x) isnumeric(x) && numel(x)==3);
p.addParameter('hotScale', false, @(x) islogical(x) && isscalar(x));
p.addParameter('colorHotMaxFactor', 8.0, @(x) isnumeric(x) && isscalar(x) && x > 1);
p.addParameter('hardAlphaCutoff', false, @(x) islogical(x) && isscalar(x));
p.addParameter('timeLabelRef', 'start', @(x) ischar(x) || isstring(x));
p.addParameter('stimOnsetMs', 0, @(x) isnumeric(x) && isscalar(x) && isfinite(x));
p.addParameter('frameRate', 10, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('quality', 95, @(x) isnumeric(x) && isscalar(x) && x >= 1 && x <= 100);
p.addParameter('verbose', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;

assert(isstruct(G) && isfield(G, 'groupMeanSigned') && isfield(G, 'timeWindows'), ...
    'G must be grouped output from build_grouped_alonggc_polygons_allbins.');
nFrames = size(G.groupMeanSigned, 2);

if isempty(opt.colorRedAt)
    redAtUse = opt.alphaFullAt;
else
    redAtUse = opt.colorRedAt;
end

if opt.verbose
    fprintf('Writing grouped polygon movie: %s\n', outMovie);
    fprintf(['  nFrames=%d | alphaFullAt=%.6g | colorRedAt=%.6g | cMaxFixed=%s | ' ...
             'stimOnsetMs=%.3g | hardAlphaCutoff=%d\n'], ...
        nFrames, opt.alphaFullAt, redAtUse, mat2str(opt.cMaxFixed), ...
        opt.stimOnsetMs, opt.hardAlphaCutoff);
end

vw = VideoWriter(outMovie, 'MPEG-4');
vw.FrameRate = opt.frameRate;
vw.Quality = round(opt.quality);
open(vw);

fracAbove = nan(nFrames,1);
for tb = 1:nFrames
    tw = double(G.timeWindows(tb,:));
    showStimulus = true;
    if all(isfinite(tw))
        showStimulus = (tw(1) >= opt.stimOnsetMs);
    end

    h = plot_grouped_alonggc_polygon_frame( ...
        stimID_example, ALLCOORDS, RTAB384, G, tb, ...
        'alphaFullAt', opt.alphaFullAt, ...
        'colorRedAt', redAtUse, ...
        'cMaxFixed', opt.cMaxFixed, ...
        'alpha', opt.alpha, ...
        'bgColor', opt.bgColor, ...
        'cLow', opt.cLow, ...
        'cHigh', opt.cHigh, ...
        'hotScale', opt.hotScale, ...
        'colorHotMaxFactor', opt.colorHotMaxFactor, ...
        'hardAlphaCutoff', opt.hardAlphaCutoff, ...
        'timeLabelRef', opt.timeLabelRef, ...
        'showStimulus', showStimulus);

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
M.alphaFullAt = opt.alphaFullAt;
M.colorRedAt = redAtUse;
M.cMaxFixed = opt.cMaxFixed;
M.stimOnsetMs = opt.stimOnsetMs;
M.hardAlphaCutoff = opt.hardAlphaCutoff;
M.fracAboveByFrame = fracAbove;
M.meanFracAbove = mean(fracAbove, 'omitnan');

if opt.verbose
    fprintf('Grouped movie done: %s\n', outMovie);
    fprintf('Mean frame exceedance >thr: %.2f%%\n', 100*M.meanFracAbove);
end

end
