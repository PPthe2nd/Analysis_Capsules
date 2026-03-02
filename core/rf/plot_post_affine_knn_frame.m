function h = plot_post_affine_knn_frame(stimID_example, ALLCOORDS, RTAB384, binData, varargin)
% PLOT_POST_AFFINE_KNN_FRAME
% Plot one post-affine frame using the same per-bin KNN value pipeline as prep.
%
% The value shown per point is abs(mean of signed delta over K nearest points
% in post-affine pixel space).

p = inputParser;
p.addParameter('K', 30, @(x) isnumeric(x) && isscalar(x) && x >= 1);
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
p.addParameter('enforceK', true, @(x) islogical(x) && isscalar(x));
p.addParameter('timeWindow', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
p.addParameter('timeLabelRef', 'center', @(x) ischar(x) || isstring(x));
p.addParameter('showStimulus', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;

[x, y, vSigned, stream] = smooth_post_affine_bin_knn(binData, round(opt.K), opt.enforceK);
V = route_by_stream(vSigned, stream);

if isempty(opt.cMaxFixed)
    cMax = prctile(V, 95);
    if ~isfinite(cMax) || cMax <= 0
        cMax = max(V);
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

r = max(0, V ./ redAt);
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
    C = zeros(numel(tColor), 3);
    for j = 1:3
        C(:,j) = interp1(tk, ck(:,j), tColor, 'linear', 'extrap');
    end
else
    tr = min(max(r, 0), 1);
    C = (1-tr).*opt.cLow + tr.*opt.cHigh;
end

if opt.hardAlphaCutoff
    alphaScale = double(V >= opt.alphaFullAt);
else
    alphaScale = min(max(V ./ opt.alphaFullAt, 0), 1);
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

if ~isempty(x)
    hSc = scatter(ax, x, y, opt.markerSize, C, 'filled');
    hSc.MarkerEdgeColor = 'none';
    try
        hSc.MarkerFaceAlpha = 'flat';
        hSc.AlphaData = alphaPoint;
        hSc.AlphaDataMapping = 'none';
        if isprop(hSc, 'MarkerEdgeAlpha')
            hSc.MarkerEdgeAlpha = 'flat';
        end
    catch
        delete(hSc);
        scatter_alpha_fallback(ax, x, y, C, alphaPoint, opt.markerSize);
    end
end

xlim(ax, [1 W]);
ylim(ax, [1 H]);
axis(ax, 'equal');
set(ax, 'YDir', 'reverse');

hFrame = rectangle(ax, 'Position', [0.5 0.5 W H], 'EdgeColor', [0.85 0.85 0.85], 'LineWidth', 1);
uistack(hFrame, 'top');

hTime = [];
if ~isempty(opt.timeWindow)
    tw = double(opt.timeWindow(:)');
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
end

h = struct();
h.fig = fig;
h.ax = ax;
h.hTime = hTime;
h.x = x;
h.y = y;
h.stream = stream;
h.VSigned = vSigned;
h.VAbs = V;
h.cMax = cMax;
h.threshold = opt.alphaFullAt;
h.fracAboveThreshold = mean(V > opt.alphaFullAt);
h.showStimulus = opt.showStimulus;
h.hardAlphaCutoff = opt.hardAlphaCutoff;

end

function scatter_alpha_fallback(ax, X, Y, C, alphaPoint, markerSize)
if isempty(X)
    return;
end

alphaQ = min(1, max(0, round(alphaPoint*10)/10));
[aq, ~, g] = unique(alphaQ);
for i = 1:numel(aq)
    a = aq(i);
    idx = (g == i);
    if ~any(idx)
        continue;
    end
    h = scatter(ax, X(idx), Y(idx), markerSize, C(idx,:), 'filled');
    h.MarkerEdgeColor = 'none';
    h.MarkerFaceAlpha = a;
    if isprop(h, 'MarkerEdgeAlpha')
        h.MarkerEdgeAlpha = a;
    end
end
end

function V = route_by_stream(vSigned, stream)
% Match target/distractor routing semantics:
% stream 1 (target-assigned):     max(0, +delta)
% stream 2 (distractor-assigned): max(0, -delta)
% unknown stream:                 abs(delta)
V = zeros(size(vSigned));
isT = (stream == 1);
isD = (stream == 2);
V(isT) = max(0,  vSigned(isT));
V(isD) = max(0, -vSigned(isD));
isOther = ~(isT | isD);
V(isOther) = abs(vSigned(isOther));
end
