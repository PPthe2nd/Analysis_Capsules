function h = plot_grouped_alonggc_polygon_frame(stimID_example, ALLCOORDS, RTAB384, G, timeIdx, varargin)
% PLOT_GROUPED_ALONGGC_POLYGON_FRAME
% Plot one frame where each along_GC group is rendered as a polygon patch.
% Group routing:
%   mean signed value >= 0 -> target polygon (stream 1)
%   mean signed value < 0  -> distractor polygon (stream 2)

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
p.addParameter('showStimulus', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;

assert(isstruct(G) && isfield(G, 'groupMeanSigned') && isfield(G, 'groupPolygons') && isfield(G, 'timeWindows'), ...
    'G must be grouped output from build_grouped_alonggc_polygons_allbins.');

nGroups = size(G.groupMeanSigned, 1);
nBins = size(G.groupMeanSigned, 2);
assert(timeIdx >= 1 && timeIdx <= nBins, 'timeIdx=%d out of range [1..%d].', timeIdx, nBins);
assert(size(G.groupPolygons,1) == nGroups && size(G.groupPolygons,2) == 2, ...
    'G.groupPolygons must be [nGroups x 2].');

vSigned = double(G.groupMeanSigned(:, timeIdx));
absV = abs(vSigned);
streamSel = ones(nGroups,1);
streamSel(vSigned < 0) = 2;
streamSel(~isfinite(vSigned)) = 0;

if isempty(opt.cMaxFixed)
    cMax = prctile(absV, 95);
    if ~isfinite(cMax) || cMax <= 0
        cMax = max(absV);
    end
    if ~isfinite(cMax) || cMax <= 0
        cMax = 1;
    end
else
    cMax = opt.cMaxFixed;
end

if isempty(opt.colorRedAt)
    redAt = opt.alphaFullAt;
else
    redAt = opt.colorRedAt;
end

r = max(0, absV ./ redAt);
if opt.hotScale
    hotCap = max(1.01, opt.colorHotMaxFactor);
    tColor = min(r, hotCap);
    tk = [0.00, 1.00, min(2.0, hotCap), hotCap];
    ck = [opt.cLow;
          [0.85 0.05 0.05];
          [1.00 0.90 0.20];
          [1.00 1.00 1.00]];
    if hotCap <= 2
        tk = [0.00, 1.00, hotCap];
        ck = [opt.cLow;
              [0.85 0.05 0.05];
              [1.00 1.00 1.00]];
    end
    C = zeros(nGroups, 3);
    for j = 1:3
        C(:,j) = interp1(tk, ck(:,j), tColor, 'linear', 'extrap');
    end
else
    tr = min(max(r, 0), 1);
    C = (1-tr).*opt.cLow + tr.*opt.cHigh;
end

if opt.hardAlphaCutoff
    alphaScale = double(absV >= opt.alphaFullAt);
else
    alphaScale = min(max(absV ./ opt.alphaFullAt, 0), 1);
end
alphaPoint = min(1, max(0, opt.alpha * alphaScale));

W = 1024; H = 768;
if opt.showStimulus
    img = render_stim_from_ALLCOORDS(ALLCOORDS, RTAB384, stimID_example);
    [H, W, ~] = size(img);
else
    bg = reshape(uint8(255 * max(0, min(1, opt.bgColor(:)'))), 1, 1, 3);
    img = repmat(bg, H, W);
end

fig = figure('Color', opt.bgColor);
ax = axes('Position', [0 0 1 1]); hold(ax, 'on');
imshow(img, 'Parent', ax, 'InitialMagnification', 'fit');
set(ax, 'Position', [0 0 1 1], 'Color', opt.bgColor);
axis(ax, 'ij');

% Draw weak groups first, strong groups last.
[~, ord] = sort(absV, 'ascend');
hPatch = gobjects(0);
for ii = 1:numel(ord)
    g = ord(ii);
    if streamSel(g) == 0 || ~isfinite(alphaPoint(g)) || alphaPoint(g) <= 0
        continue;
    end
    P = G.groupPolygons(g, streamSel(g));
    if isempty(P.x) || isempty(P.y)
        continue;
    end
    if numel(P.x) < 3
        continue;
    end

    col = C(g,:);
    a = alphaPoint(g);
    hp = patch(ax, P.x(:), P.y(:), col, ...
        'FaceColor', col, ...
        'FaceAlpha', a, ...
        'EdgeColor', col, ...
        'LineWidth', 1.2);
    if isprop(hp, 'EdgeAlpha')
        hp.EdgeAlpha = min(1, max(0, a));
    end
    hPatch(end+1,1) = hp; %#ok<AGROW>
end

xlim(ax, [1 W]);
ylim(ax, [1 H]);
axis(ax, 'equal');
set(ax, 'YDir', 'reverse');
hFrame = rectangle(ax, 'Position', [0.5 0.5 W H], 'EdgeColor', [0.85 0.85 0.85], 'LineWidth', 1);
uistack(hFrame, 'top');

hTime = [];
tw = double(G.timeWindows(timeIdx,:));
if all(isfinite(tw))
    ref = lower(string(opt.timeLabelRef));
    switch ref
        case "start"
            tLbl = tw(1);
        case "end"
            tLbl = tw(2);
        otherwise
            tLbl = mean(tw);
    end
    hTime = text(ax, 14, 14, sprintf('%d ms', round(tLbl)), ...
        'Color', [1 1 1], 'FontSize', 14, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
end

h = struct();
h.fig = fig;
h.ax = ax;
h.hPatch = hPatch;
h.hTime = hTime;
h.groupSigned = vSigned;
h.groupAbs = absV;
h.groupStream = streamSel;
h.alphaPoint = alphaPoint;
h.threshold = opt.alphaFullAt;
h.fracAboveThreshold = mean(absV > opt.alphaFullAt, 'omitnan');
h.nTargetGroups = nnz(streamSel == 1 & isfinite(vSigned));
h.nDistrGroups = nnz(streamSel == 2 & isfinite(vSigned));
h.showStimulus = opt.showStimulus;
h.cMax = cMax;

end
