% Build stimulus-aligned RF geometry tables for IT sites in the line task.
% Mirrors the V4 setup in Line_Stimuli_V4.m, but keeps IT as a separate
% entry point so the existing V1/V4 workflows stay unchanged.

Monkey = 1; % 1 for Nilson, 2 for Figaro
TabFile = "ObjAtt_lines_monkeyN_20220201_B1";
ExampleStimulus = 68;
MakeExamplePlot = true;
ForceRebuild = false;

cfg = config();
repoRoot = cfg.repoRoot;
logDir = cfg.logsDir;

if Monkey == 1
    monkeySuffix = "N";
elseif Monkey == 2
    monkeySuffix = "F";
else
    error('Line_Stimuli_IT:InvalidMonkey', ...
        'Monkey must be 1 (Nilson) or 2 (Figaro).');
end

outFile = fullfile(cfg.matDir, sprintf('Tall_IT_lines_%s.mat', char(monkeySuffix)));
skipBuild = exist(outFile, 'file') == 2 && ~ForceRebuild;

cd(repoRoot);
load(fullfile(logDir, TabFile));      % loads ALLCOORDS
load(fullfile(logDir, 'RTAB384.mat'));

% Load RF centers/sizes and monkey-specific area ranges into the workspace.
RFs

RFrange = IT;
if isempty(RFrange)
    error('Line_Stimuli_IT:NoITSites', ...
        'No IT site indices are available for monkey %d.', Monkey);
end

x_rf = x(RFrange);
y_rf = y(RFrange);

fprintf('Using %d IT sites for monkey %s.\n', numel(RFrange), MonkeyName);

% Example single-stimulus geometry table for quick inspection.
stimNum = ExampleStimulus;
T = rf_table_target_distractor(ALLCOORDS, RTAB384, stimNum, x_rf, y_rf);
idxNaN = isnan(T.x_px) | isnan(T.y_px);
T.assignment(idxNaN) = "NaN";
T = add_arm_projection_metrics(T, ALLCOORDS, stimNum);
T = add_polar_about_s(T, ALLCOORDS, stimNum, 'OnlyBackground', true);
T = add_arc_about_s_edges(T, ALLCOORDS, RTAB384, stimNum);
[T, widthPx] = add_GC_normalization(T, RTAB384, stimNum); %#ok<NASGU>

fprintf(['Example stim %d: target=%d, distractor=%d, background=%d, ' ...
         'overlap=%d\n'], ...
    stimNum, ...
    nnz(T.assignment == "target"), ...
    nnz(T.assignment == "distractor"), ...
    nnz(T.assignment == "background"), ...
    nnz(T.overlap));

if MakeExamplePlot
    objSites = T.RF(T.assignment == "target" | T.assignment == "distractor");
    h = plot_stim_with_RFs(ALLCOORDS, RTAB384, stimNum, x_rf(objSites), y_rf(objSites)); %#ok<NASGU>
    title(sprintf('IT object-assigned RF centers on stim %d', stimNum));
end

% Build the full per-stimulus geometry table for all 384 line-task stimuli.
if skipBuild
    fprintf('Skipping full IT geometry build because output already exists:\n%s\n', outFile);
    fprintf('Set ForceRebuild = true to rebuild Tall_IT.\n');
    return;
end

Tall_IT = build_all_stim_tables(ALLCOORDS, RTAB384, x_rf, y_rf);

save(outFile, 'Tall_IT', 'ALLCOORDS', 'RTAB384', 'Monkey', 'MonkeyName', ...
    'RFrange', 'TabFile', '-v7.3');

fprintf('Saved IT geometry table to %s\n', outFile);
