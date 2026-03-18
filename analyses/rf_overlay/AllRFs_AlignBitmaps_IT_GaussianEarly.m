function RES = AllRFs_AlignBitmaps_IT_GaussianEarly(Monkey, Opts)
% Align early Gaussian IT RF centers onto one reference bitmap.

if nargin < 1 || isempty(Monkey)
    Monkey = 1; % 1 = Nilson, 2 = Figaro
end
if nargin < 2 || isempty(Opts)
    Opts = struct();
end

Opts = normalize_opts_local(Opts);
cfg = config();

if Monkey == 1
    monkeySuffix = "N";
elseif Monkey == 2
    monkeySuffix = "F";
else
    error('AllRFs_AlignBitmaps_IT_GaussianEarly:InvalidMonkey', ...
        'Monkey must be 1 (Nilson) or 2 (Figaro).');
end

gaussPath = fullfile(cfg.matDir, sprintf('GaussianOccupancy_Tuning_IT_%s.mat', char(monkeySuffix)));
assert(exist(gaussPath, 'file') == 2, ...
    'Missing %s. Run GaussianOccupancy_Tuning_IT.m first.', gaussPath);

Sg = load(gaussPath, 'OUT');
assert(isfield(Sg, 'OUT') && isstruct(Sg.OUT), '%s must contain struct OUT.', gaussPath);
OUT = Sg.OUT;
assert(isfield(OUT, 'FitSpatialEarly') && isfield(OUT, 'RFrange'), ...
    'Gaussian OUT is missing FitSpatialEarly or RFrange.');

FitEarly = OUT.FitSpatialEarly;
RFrange = OUT.RFrange(:);
gxEarly = [FitEarly.centerX].';
gyEarly = [FitEarly.centerY].';
sigmaEarly = [FitEarly.sigmaPx].';
veEarly = [FitEarly.r2TrainPct].';
pEarly = [FitEarly.pValueApprox].';

keepMask = isfinite(gxEarly) & isfinite(gyEarly) & ...
    isfinite(sigmaEarly) & (sigmaEarly > 0) & ...
    isfinite(veEarly) & (veEarly >= Opts.MinSpatialVE) & ...
    (sigmaEarly <= Opts.MaxSigmaPx);
if Opts.RequireSpatialSig
    keepMask = keepMask & isfinite(pEarly) & (pEarly < Opts.Palpha);
end

assert(any(keepMask), ...
    'No early Gaussian IT RFs passed the requested filters.');

x_rf = gxEarly(keepMask);
y_rf = gyEarly(keepMask);
siteGlobal = RFrange(keepMask);
nSites = numel(siteGlobal);

bitmapNums = choose_bitmap_numbers_local(cfg.stimDir, Opts);
assert(any(bitmapNums == Opts.ReferenceBitmap), ...
    'Reference bitmap %d is not in the selected bitmap list.', Opts.ReferenceBitmap);

refFile = fullfile(cfg.stimDir, sprintf('%03d.bmp', Opts.ReferenceBitmap));
assert(exist(refFile, 'file') == 2, 'Missing reference bitmap: %s', refFile);
Iref = imread(refFile);
[Href, Wref, ~] = size(Iref);
Fref = extract_capsule_axes(Iref);

nBmp = numel(bitmapNums);
All_RFs_x = nan(nSites, nBmp);
All_RFs_y = nan(nSites, nBmp);

perBitmap = repmat(struct( ...
    'bmp', NaN, 'ok', false, 'mirror', "", 'rotDeg', NaN, ...
    'errY', NaN, 'errP', NaN, 'nInBounds', 0, 'message', ""), nBmp, 1);

fprintf('Early Gaussian IT RF affine overlay (%s)\n', char(monkeySuffix));
fprintf('Using %d / %d early IT Gaussian centers\n', nSites, numel(RFrange));
fprintf('Filters: %s\n', char(describe_filters_local(Opts)));
fprintf('Reference bitmap: %03d | moving bitmaps: %d\n', Opts.ReferenceBitmap, nBmp);

for iBmp = 1:nBmp
    movingIdx = bitmapNums(iBmp);
    perBitmap(iBmp).bmp = movingIdx;

    try
        movFile = fullfile(cfg.stimDir, sprintf('%03d.bmp', movingIdx));
        assert(exist(movFile, 'file') == 2, 'Missing bitmap %03d.', movingIdx);
        Imov = imread(movFile);

        if movingIdx == Opts.ReferenceBitmap
            Fmov_used = Fref;
            outAlign = struct('ok', true, 'mirror', "none", 'rotDeg', 0, ...
                'errY', 0, 'errP', 0, 'tform', affine2d(eye(3)));
        else
            Fmov = extract_capsule_axes(Imov);
            [outAlign, Fmov_used] = robust_align_local(Fmov, Fref, Opts.rotTolDeg, Opts.centTolPx);
        end

        if Opts.UseDetectedStimCenter
            Cmov = Fmov_used.stimCenter_xy(:)';
        else
            Cmov = Opts.ScreenCenter(:)';
        end
        x_pix = Cmov(1) + x_rf(:);
        if Opts.flipY
            y_pix = Cmov(2) - y_rf(:);
        else
            y_pix = Cmov(2) + y_rf(:);
        end

        [x_ref, y_ref] = transformPointsForward(outAlign.tform, x_pix, y_pix);
        inBounds = isfinite(x_ref) & isfinite(y_ref) & ...
            x_ref >= 1 & x_ref <= Wref & y_ref >= 1 & y_ref <= Href;

        All_RFs_x(:, iBmp) = x_ref(:);
        All_RFs_y(:, iBmp) = y_ref(:);

        perBitmap(iBmp).ok = true;
        perBitmap(iBmp).mirror = mirror_label_local(outAlign.mirror);
        perBitmap(iBmp).rotDeg = outAlign.rotDeg;
        perBitmap(iBmp).errY = outAlign.errY;
        perBitmap(iBmp).errP = outAlign.errP;
        perBitmap(iBmp).nInBounds = nnz(inBounds);
        perBitmap(iBmp).message = "";

        if Opts.Verbose && (iBmp == 1 || mod(iBmp, Opts.ProgressEvery) == 0 || iBmp == nBmp)
            fprintf('  Bitmap %03d (%d / %d): in-bounds=%d | mirror=%s | rot=%.2f\n', ...
                movingIdx, iBmp, nBmp, nnz(inBounds), char(perBitmap(iBmp).mirror), perBitmap(iBmp).rotDeg);
        end
    catch ME
        perBitmap(iBmp).ok = false;
        perBitmap(iBmp).message = string(ME.message);
        if Opts.Verbose && Opts.ReportFailures
            fprintf('  Bitmap %03d failed: %s\n', movingIdx, ME.message);
        end
    end
end

okBmp = [perBitmap.ok].';
Xall = All_RFs_x(:, okBmp);
Yall = All_RFs_y(:, okBmp);
xPlot = Xall(:);
yPlot = Yall(:);
goodPlot = isfinite(xPlot) & isfinite(yPlot) & ...
    xPlot >= 1 & xPlot <= Wref & yPlot >= 1 & yPlot <= Href;
xPlot = xPlot(goodPlot);
yPlot = yPlot(goodPlot);

fprintf('Successful bitmap alignments: %d / %d\n', nnz(okBmp), nBmp);
fprintf('Total transformed IT Gaussian points plotted on %03d: %d\n', Opts.ReferenceBitmap, numel(xPlot));

figTitle = sprintf('Early Gaussian IT RFs on reference %03d (%s)', Opts.ReferenceBitmap, char(monkeySuffix));
figure(Opts.FigureNumber); clf;
set(gcf, 'Color', 'w', 'Name', figTitle, 'NumberTitle', 'off', ...
    'Tag', 'IT_GaussianEarly_affine_overlay');
ax = axes();
imshow(Iref, 'Parent', ax);
hold(ax, 'on');
axis(ax, 'image');

hRF = scatter(ax, xPlot, yPlot, Opts.MarkerSize, Opts.MarkerColor, 'filled', ...
    'MarkerFaceAlpha', Opts.MarkerAlpha, 'MarkerEdgeColor', 'k', ...
    'MarkerEdgeAlpha', min(1, Opts.MarkerAlpha + 0.10));

Cref = Fref.stimCenter_xy(:)';
plot(ax, Cref(1), Cref(2), 'wx', 'MarkerSize', 12, 'LineWidth', 5);
hCtr = plot(ax, Cref(1), Cref(2), 'kx', 'MarkerSize', 12, 'LineWidth', 2);
title(ax, sprintf('%s | sites=%d | bitmaps ok=%d/%d | points=%d', ...
    figTitle, nSites, nnz(okBmp), nBmp, numel(xPlot)));
legend(ax, [hRF, hCtr], {'Early Gaussian IT RFs', 'Stim center'}, 'Location', 'best');

RFmap = build_rfmap_table_local(bitmapNums, okBmp, All_RFs_x, All_RFs_y, siteGlobal);

RES = struct();
RES.monkeySuffix = monkeySuffix;
RES.gaussPath = gaussPath;
RES.referenceBitmap = Opts.ReferenceBitmap;
RES.bitmapNums = bitmapNums(:);
RES.siteGlobal = siteGlobal(:);
RES.x_rf = x_rf(:);
RES.y_rf = y_rf(:);
RES.keepMask = keepMask(:);
RES.filters = Opts;
RES.perBitmap = perBitmap;
RES.All_RFs_x = All_RFs_x;
RES.All_RFs_y = All_RFs_y;
RES.RFmap = RFmap;
RES.refStimCenter = Cref;

if Opts.SaveResult
    outFile = fullfile(cfg.matDir, sprintf('AllRFs_AlignBitmaps_IT_GaussianEarly_%s.mat', char(monkeySuffix)));
    save(outFile, 'RES', '-v7.3');
    fprintf('Saved affine-overlay result to %s\n', outFile);
    RES.outFile = outFile;
end
end

function RFmap = build_rfmap_table_local(bitmapNums, okBmp, All_RFs_x, All_RFs_y, siteGlobal)
bmpOk = bitmapNums(okBmp);
nSites = numel(siteGlobal);
nBmp = numel(bmpOk);

if nBmp == 0
    RFmap = table();
    return;
end

RFmap = table( ...
    repelem(bmpOk(:), nSites, 1), ...
    repmat((1:nSites).', nBmp, 1), ...
    repmat(siteGlobal(:), nBmp, 1), ...
    reshape(All_RFs_x(:, okBmp), [], 1), ...
    reshape(All_RFs_y(:, okBmp), [], 1), ...
    'VariableNames', {'bmp', 'rfId', 'siteGlobal', 'x_ref', 'y_ref'});
end

function bitmapNums = choose_bitmap_numbers_local(stimDir, Opts)
if ~isempty(Opts.BitmapNums)
    bitmapNums = unique(double(Opts.BitmapNums(:).'));
    return;
end

nums = 1:Opts.MaxBitmapNumberDefault;
if Opts.UseOddBitmapsOnly
    nums = nums(mod(nums, 2) == 1);
end
files = dir(fullfile(stimDir, '*.bmp'));
available = nan(numel(files), 1);
for i = 1:numel(files)
    tok = regexp(files(i).name, '^(\d+)\.bmp$', 'tokens', 'once');
    if ~isempty(tok)
        available(i) = str2double(tok{1});
    end
end
available = unique(available(isfinite(available)));
nums = nums(ismember(nums, available));
bitmapNums = nums(:).';
end

function [outAlign, Fmov_used] = robust_align_local(Fmov, Fref, rotTolDeg, centTolPx)
Fmov_used = Fmov;

try
    outAlign = align_centroids_about_center(Fmov_used, Fref, ...
        'rotTolDeg', rotTolDeg, 'centTolPx', centTolPx);
    return
catch ME
end

Fswap = Fmov_used;
tmp = Fswap.yellow;
Fswap.yellow = Fswap.purple;
Fswap.purple = tmp;

try
    outAlign = align_centroids_about_center(Fswap, Fref, ...
        'rotTolDeg', rotTolDeg, 'centTolPx', centTolPx);
    Fmov_used = Fswap;
    return
catch ME2
    error(['Alignment failed even after swapping Y/P in moving. ' ...
        'Original error: %s | Swap retry error: %s'], ME.message, ME2.message);
end
end

function txt = mirror_label_local(val)
if islogical(val)
    if val
        txt = "mirror";
    else
        txt = "none";
    end
else
    txt = string(val);
end
end

function Opts = normalize_opts_local(Opts)
defaults = struct();
defaults.ReferenceBitmap = 61;
defaults.BitmapNums = [];
defaults.UseOddBitmapsOnly = true;
defaults.MaxBitmapNumberDefault = 183;
defaults.flipY = true;
defaults.ScreenCenter = [512 384];
defaults.UseDetectedStimCenter = false;
defaults.rotTolDeg = 14;
defaults.centTolPx = 15;
defaults.RequireSpatialSig = true;
defaults.Palpha = 0.05;
defaults.MinSpatialVE = -inf;
defaults.MaxSigmaPx = inf;
defaults.MarkerSize = 14;
defaults.MarkerAlpha = 0.10;
defaults.MarkerColor = [0.60 0.60 0.60];
defaults.FigureNumber = 106;
defaults.ProgressEvery = 10;
defaults.Verbose = true;
defaults.ReportFailures = false;
defaults.SaveResult = false;

fn = fieldnames(defaults);
for i = 1:numel(fn)
    if ~isfield(Opts, fn{i}) || isempty(Opts.(fn{i}))
        Opts.(fn{i}) = defaults.(fn{i});
    end
end
end

function txt = describe_filters_local(Opts)
parts = {};
if Opts.RequireSpatialSig
    parts{end+1} = sprintf('early p < %.2f', Opts.Palpha); %#ok<AGROW>
else
    parts{end+1} = 'finite early Gaussian fits'; %#ok<AGROW>
end
if isfinite(Opts.MinSpatialVE)
    parts{end+1} = sprintf('VE >= %.1f%%', Opts.MinSpatialVE); %#ok<AGROW>
end
if isfinite(Opts.MaxSigmaPx)
    parts{end+1} = sprintf('sigma <= %.1f px', Opts.MaxSigmaPx); %#ok<AGROW>
end
if Opts.UseOddBitmapsOnly && isempty(Opts.BitmapNums)
    parts{end+1} = sprintf('odd-numbered bitmaps only up to %d', Opts.MaxBitmapNumberDefault); %#ok<AGROW>
end
if Opts.UseDetectedStimCenter
    parts{end+1} = 'using detected stimulus center'; %#ok<AGROW>
else
    parts{end+1} = sprintf('using fixed screen center [%g %g]', ...
        Opts.ScreenCenter(1), Opts.ScreenCenter(2)); %#ok<AGROW>
end
txt = strjoin(parts, ', ');
end
