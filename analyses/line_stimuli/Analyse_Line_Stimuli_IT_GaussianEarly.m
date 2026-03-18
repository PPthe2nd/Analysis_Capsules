function OUTG = Analyse_Line_Stimuli_IT_GaussianEarly(Monkey, ExampleStimulus, Opts)
% Plot all projected IT RFs on one stimulus using early Gaussian RF centers.

if nargin < 1 || isempty(Monkey)
    Monkey = 1; % 1 = Nilson, 2 = Figaro
end
if nargin < 2 || isempty(ExampleStimulus)
    ExampleStimulus = 38;
end
if nargin < 3 || isempty(Opts)
    Opts = struct();
end

Opts = normalize_opts_local(Opts);
cfg = config();

if Monkey == 1
    monkeySuffix = "N";
    gaussFile = 'GaussianOccupancy_Tuning_IT_N.mat';
    tallBaseFile = 'Tall_IT_lines_N.mat';
elseif Monkey == 2
    monkeySuffix = "F";
    gaussFile = 'GaussianOccupancy_Tuning_IT_F.mat';
    tallBaseFile = 'Tall_IT_lines_F.mat';
else
    error('Analyse_Line_Stimuli_IT_GaussianEarly:InvalidMonkey', ...
        'Monkey must be 1 (Nilson) or 2 (Figaro).');
end

gaussPath = fullfile(cfg.matDir, gaussFile);
tallBasePath = fullfile(cfg.matDir, tallBaseFile);
assert(exist(gaussPath, 'file') == 2, ...
    'Missing %s. Run GaussianOccupancy_Tuning_IT.m first.', gaussPath);
assert(exist(tallBasePath, 'file') == 2, ...
    'Missing %s. Run Line_Stimuli_IT.m first.', tallBasePath);

Sg = load(gaussPath, 'OUT');
Sgeo = load(tallBasePath, 'ALLCOORDS', 'RTAB384');
assert(isfield(Sg, 'OUT') && isstruct(Sg.OUT), '%s must contain struct OUT.', gaussPath);
assert(isfield(Sgeo, 'ALLCOORDS') && isfield(Sgeo, 'RTAB384'), ...
    '%s must contain ALLCOORDS and RTAB384.', tallBasePath);

OUT = Sg.OUT;
ALLCOORDS = Sgeo.ALLCOORDS;
RTAB384 = Sgeo.RTAB384;

assert(isfield(OUT, 'FitSpatialEarly') && isfield(OUT, 'RFrange'), ...
    '%s must contain FitSpatialEarly and RFrange.', gaussFile);
assert(ExampleStimulus >= 1 && ExampleStimulus <= size(RTAB384, 1), ...
    'ExampleStimulus must be between 1 and %d.', size(RTAB384, 1));

FitEarly = OUT.FitSpatialEarly;
RFrange = OUT.RFrange(:);
gxEarly = [FitEarly.centerX].';
gyEarly = [FitEarly.centerY].';
sigmaEarly = [FitEarly.sigmaPx].';
veEarly = [FitEarly.r2TrainPct].';
pEarly = [FitEarly.pValueApprox].';

keepMask = isfinite(gxEarly) & isfinite(gyEarly);
if Opts.RequireFiniteSigma
    keepMask = keepMask & isfinite(sigmaEarly) & (sigmaEarly > 0);
end
if isfinite(Opts.MinSpatialVE)
    keepMask = keepMask & isfinite(veEarly) & (veEarly >= Opts.MinSpatialVE);
end
if Opts.RequireSpatialSig
    keepMask = keepMask & isfinite(pEarly) & (pEarly < Opts.Palpha);
end

assert(any(keepMask), ...
    'No early Gaussian IT RFs passed the requested filters.');

x_rf = gxEarly(keepMask);
y_rf = gyEarly(keepMask);
RFrangeKeep = RFrange(keepMask);

fprintf('Building Gaussian-early IT geometry (%s)\n', char(monkeySuffix));
fprintf('Keeping %d / %d early Gaussian IT RF centers\n', nnz(keepMask), numel(keepMask));
fprintf('Filters: %s\n', char(describe_filters_local(Opts)));

Tall_IT_GaussianEarly = build_all_stim_tables(ALLCOORDS, RTAB384, x_rf, y_rf);

h = plot_projected_RFs_on_example_stim(Tall_IT_GaussianEarly, ALLCOORDS, RTAB384, ExampleStimulus, ...
    'MarkerSize', Opts.MarkerSize, 'Alpha', Opts.Alpha);
set(h.fig, 'Name', sprintf('All projected IT RFs on stim %d (Gaussian early, %s)', ...
    ExampleStimulus, char(monkeySuffix)), 'NumberTitle', 'off', ...
    'Tag', 'IT_projected_RFs_GaussianEarly');
title(h.ax, sprintf('All projected IT RFs on stim %d (Gaussian early, %s)', ...
    ExampleStimulus, char(monkeySuffix)));

OUTG = struct();
OUTG.Monkey = Monkey;
OUTG.monkeySuffix = monkeySuffix;
OUTG.ExampleStimulus = ExampleStimulus;
OUTG.filters = Opts;
OUTG.gaussPath = gaussPath;
OUTG.tallBasePath = tallBasePath;
OUTG.RFrange = RFrange;
OUTG.keepMask = keepMask;
OUTG.RFrangeKeep = RFrangeKeep;
OUTG.x_rf = x_rf;
OUTG.y_rf = y_rf;
OUTG.Tall_IT_GaussianEarly = Tall_IT_GaussianEarly;
OUTG.plotHandle = h;

if Opts.SaveTall
    outFile = fullfile(cfg.matDir, sprintf('Tall_IT_lines_GaussianEarly_%s.mat', char(monkeySuffix)));
    save(outFile, 'Tall_IT_GaussianEarly', 'ALLCOORDS', 'RTAB384', ...
        'Monkey', 'monkeySuffix', 'RFrangeKeep', 'x_rf', 'y_rf', 'ExampleStimulus', 'Opts', '-v7.3');
    fprintf('Saved Gaussian-early IT geometry table to %s\n', outFile);
    OUTG.outFile = outFile;
end
end

function Opts = normalize_opts_local(Opts)
defaults = struct();
defaults.RequireFiniteSigma = true;
defaults.RequireSpatialSig = false;
defaults.Palpha = 0.05;
defaults.MinSpatialVE = -inf;
defaults.MarkerSize = 4;
defaults.Alpha = 0.15;
defaults.SaveTall = false;

fn = fieldnames(defaults);
for i = 1:numel(fn)
    if ~isfield(Opts, fn{i}) || isempty(Opts.(fn{i}))
        Opts.(fn{i}) = defaults.(fn{i});
    end
end
end

function txt = describe_filters_local(Opts)
parts = {};
parts{end+1} = 'finite early Gaussian centers'; %#ok<AGROW>
if Opts.RequireFiniteSigma
    parts{end+1} = 'finite sigma'; %#ok<AGROW>
end
if isfinite(Opts.MinSpatialVE)
    parts{end+1} = sprintf('VE >= %.1f%%', Opts.MinSpatialVE); %#ok<AGROW>
end
if Opts.RequireSpatialSig
    parts{end+1} = sprintf('p < %.2f', Opts.Palpha); %#ok<AGROW>
end
txt = strjoin(parts, ', ');
end
