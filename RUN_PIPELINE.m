% RUN_PIPELINE
% Entry point to initialize paths/config and manually run selected analyses.

clearvars;
close all;
clc;

% Repository root (this file is expected at the root).
repoRoot = fileparts(mfilename('fullpath'));

% Add only relevant code folders (avoid data/output directories).
addpath(repoRoot);
addpath(genpath(fullfile(repoRoot, 'analyses')));
addpath(genpath(fullfile(repoRoot, 'core')));
addpath(genpath(fullfile(repoRoot, 'utils')));

% Load configuration from config.m
cfg = config(); %#ok<NASGU>
repoRoot = cfg.repoRoot;

% -------------------------------------------------------------------------
% Run options by analysis category (uncomment what you want to execute)
% -------------------------------------------------------------------------

% LINE STIMULI
run(fullfile(repoRoot, 'analyses', 'line_stimuli', 'Line_Stimuli.m'));
run(fullfile(repoRoot, 'analyses', 'line_stimuli', 'Analyse_Line_Stimuli.m'));
% run(fullfile(repoRoot, 'analyses', 'line_stimuli', 'Attention_Line_Stimuli.m'));
% OUT_colorAtt = Attention_ColorPref_Line_Stimuli(); % interaction between color tuning and attenion

% LATENCY
% run(fullfile(repoRoot, 'analyses', 'latency', 'Latency_analysis_SVM.m'));
% run(fullfile(repoRoot, 'analyses', 'latency', 'Latency_analysis_color.m'));
% run(fullfile(repoRoot, 'analyses', 'latency', 'Latency_analysis_color_integrated.m'));
% run(fullfile(repoRoot, 'analyses', 'latency', 'Latency_analysis_execSVM_fast.m'));

% DECODING
% run(fullfile(repoRoot, 'analyses', 'decoding', 'Color_SVM.m'));
% run(fullfile(repoRoot, 'analyses', 'decoding', 'Color_decoding_V1.m'));
% run(fullfile(repoRoot, 'analyses', 'decoding', 'time_resolved_SVM.m'));

% PSTH
% run(fullfile(repoRoot, 'analyses', 'psth', 'PSTH_colorPref_V1.m'));
% run(fullfile(repoRoot, 'analyses', 'psth', 'PSTH_colorPref_V1_grayDist.m'));

% RF OVERLAY
% run(fullfile(repoRoot, 'analyses', 'rf_overlay', 'RF_on_stim_V1.m'));
% run(fullfile(repoRoot, 'analyses', 'rf_overlay', 'RF_on_stim.m'));
% run(fullfile(repoRoot, 'analyses', 'rf_overlay', 'AlignBitmapsRFs.m'));
% run(fullfile(repoRoot, 'analyses', 'rf_overlay', 'AllRFs_AlignBitmaps.m'));

% COLOR TUNING
% run(fullfile(repoRoot, 'analyses', 'color_tuning', 'ColorTuning_Capsules.m'));
% run(fullfile(repoRoot, 'analyses', 'color_tuning', 'color_tuning_cone.m'));

% STIMULUS GENERATION
% run(fullfile(repoRoot, 'analyses', 'stimulus_generation', 'Eval_RANDTAB.m'));
% run(fullfile(repoRoot, 'analyses', 'stimulus_generation', 'makle_randatab_shapes_compl.m'));

% EXPLORATORY
% run(fullfile(repoRoot, 'analyses', 'exploratory', 'Find_stimulus_pairs.m'));
% run(fullfile(repoRoot, 'analyses', 'exploratory', 'Compare_Paolos_SNR_Color_SNR.m'));
% run(fullfile(repoRoot, 'analyses', 'exploratory', 'test_incl_criteria.m'));

% MOVIE
% run(fullfile(repoRoot, 'analyses', 'movie', 'make_activity_movie.m'));
% run(fullfile(repoRoot, 'analyses', 'movie', 'make_activity_movie_wrapper_safe.m'));
