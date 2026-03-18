% Analysing the RFs in the two monkeys

% clear all
% close all

Monkey=1; % (1=Nilson; 2= Figaro)
cfg = config();
cd(cfg.repoRoot);

if Monkey  == 1     %     --- Region indices Nilson ---
    V1 = 1:512;
    V4 = 513:768;
    IT = 769:1024;
    rfFile = resolve_rf_file(cfg, 'THINGS_RF1s_N.mat', 'Mr Nilson');
    load(rfFile);
    MonkeyName="N";
else%                   --- Region indices Figaro ---
    V1 = 1:512;
    IT = 513:832;
    V4 = 833:1024;
    rfFile = resolve_rf_file(cfg, 'THINGS_RF1s_F.mat', 'Figaro');
    load(rfFile);
    MonkeyName="F";
end


N = numel(all_centrex);
fprintf('Loaded %s: N=%d (centrex)\n', MonkeyName, N);

% Clip indices so they never exceed the data length
V1 = V1(V1 <= N);
V4 = V4(V4 <= N);
IT = IT(IT <= N);

fprintf('Using indices: V1=%d, V4=%d, IT=%d\n', numel(V1), numel(V4), numel(IT));


%% RFs.m — Plot RF rectangles by area (V1/V4/IT), clean figure (no grid, no ecc rings)
% Assumes variables exist in workspace:
%   all_centrex, all_centrey, all_szx, all_szy
%
% Regions:
%   V1 = 1:512; V4 = 513:768; IT = 769:1024;


k = 0.2;                  % scale factor for RF box size (set to 2 for 2*sigma, etc.)
pad = 1;                  % padding (deg) around plot limits

% Colors (RGB)
cV1 = [0.20 0.45 0.90];   % V1 blue
cV4 = [0.90 0.45 0.20];   % V4 orange
cIT = [0.60 0.60 0.60];   % IT light gray (light)

% Rectangle styling
lw_V1 = 0.7;   edgeA_V1 = 0.30;  faceA_V1 = 0.03;
lw_V4 = 0.7;   edgeA_V4 = 0.25;  faceA_V4 = 0.025;
lw_IT = 0.4;   edgeA_IT = 0.12;  faceA_IT = 0.01;

% Center marker styling (make visible)
sV1 = 14;  aV1 = 0.25;
sV4 = 18;  aV4 = 0.30;
% (skip IT centers by default)

%% --- Data vectors ---
x  = all_centrex(:);
y  = all_centrey(:);
sx = all_szx(:);
sy = all_szy(:);

w = k*sx;
h = k*sy;

valid = isfinite(x) & isfinite(y) & isfinite(w) & isfinite(h) & w>0 & h>0;

fprintf('Valid RFs: %d / %d\n', sum(valid), numel(x));
fprintf('x range: [%g %g]\n', min(x(valid)), max(x(valid)));
fprintf('y range: [%g %g]\n', min(y(valid)), max(y(valid)));
fprintf('w range: [%g %g]\n', min(w(valid)), max(w(valid)));
fprintf('h range: [%g %g]\n', min(h(valid)), max(h(valid)));

% Keep indexing intact; mark invalid entries
bad = ~isfinite(x) | ~isfinite(y) | ~isfinite(w) | ~isfinite(h) | w<=0 | h<=0;

% NEW: reject absurd center values (likely failed fits)
bad = bad | abs(x) > 200 | abs(y) > 200;

x(bad) = NaN; y(bad) = NaN; w(bad) = NaN; h(bad) = NaN;



%% --- Figure setup ---
figure('Color','w'); hold on;
axis equal;
set(gca,'Box','off','TickDir','out','XGrid','off','YGrid','off','Layer','top');
xlabel('x (deg)'); ylabel('y (deg)');
title("RF locations by area (rectangles) "+MonkeyName);

%% --- Draw rectangles (back to front) ---
draw_rf_rects(IT, x, y, w, h, cIT, edgeA_IT, faceA_IT, lw_IT);
draw_rf_rects(V4, x, y, w, h, cV4, edgeA_V4, faceA_V4, lw_V4);
draw_rf_rects(V1, x, y, w, h, cV1, edgeA_V1, faceA_V1, lw_V1);

%% --- Centers (more visible) ---
scatter(x(V1), y(V1), sV1, cV1, 'filled', ...
    'MarkerFaceAlpha', aV1, 'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', aV1);

scatter(x(V4), y(V4), sV4, cV4, 'filled', ...
    'MarkerFaceAlpha', aV4, 'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', aV4);

% Optional fixation cross at (0,0)
plot([-2 2], [0 0], 'k-', 'LineWidth', 1, 'HandleVisibility','off');
plot([0 0], [-2 2], 'k-', 'LineWidth', 1, 'HandleVisibility','off');

%% --- Limits ---
valid = isfinite(x) & isfinite(y) & isfinite(w) & isfinite(h) & w>0 & h>0;

xmin = min(x(valid) - w(valid)/2) - pad;
xmax = max(x(valid) + w(valid)/2) + pad;
ymin = min(y(valid) - h(valid)/2) - pad;
ymax = max(y(valid) + h(valid)/2) + pad;

xlim([xmin xmax]); ylim([ymin ymax]);

%% --- Legend (use dummy handles; avoids warnings) ---
hIT = plot(nan, nan, '-', 'Color', cIT, 'LineWidth', 2);
hV4 = plot(nan, nan, '-', 'Color', cV4, 'LineWidth', 2);
hV1 = plot(nan, nan, '-', 'Color', cV1, 'LineWidth', 2);
legend([hIT hV4 hV1], {'IT','V4','V1'}, 'Location','southeast');

%% ========================= Helper function (place at end) =========================
function draw_rf_rects(idx, x, y, w, h, col, edgeAlpha, faceAlpha, lw)
% Draw axis-aligned RF rectangles centered at (x,y) with width w and height h
% Skips NaNs to preserve original indexing.

for ii = 1:numel(idx)
    i = idx(ii);

    if i < 1 || i > numel(x)
        continue
    end
    if isnan(x(i)) || isnan(y(i)) || isnan(w(i)) || isnan(h(i))
        continue
    end

    rectangle('Position', [x(i)-w(i)/2, y(i)-h(i)/2, w(i), h(i)], ...
        'EdgeColor', [col edgeAlpha], ...
        'FaceColor', [col faceAlpha], ...
        'LineWidth', lw, ...
        'HandleVisibility','off');   % prevent rectangles from flooding legend
end
end

function rfPath = resolve_rf_file(cfg, fileName, monkeySubdir)
% Prefer repo-local data_mat, but allow legacy location under dataRoot.
candidates = { ...
    fullfile(cfg.matDir, fileName), ...
    fullfile(cfg.dataRoot, monkeySubdir, fileName) ...
};

rfPath = '';
for i = 1:numel(candidates)
    if exist(candidates{i}, 'file') == 2
        rfPath = candidates{i};
        return;
    end
end

error('RFs:MissingRFFile', ...
    ['Could not find %s. Checked:\n - %s\n - %s'], ...
    fileName, candidates{1}, candidates{2});
end
