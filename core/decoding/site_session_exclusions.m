function T = site_session_exclusions(monkey)
% SITE_SESSION_EXCLUSIONS
% Central registry of site-by-session exclusions. This is the source of
% truth for cases where one recording session is unusable for a site, but
% the other session should be preserved.
%
% Columns:
%   monkeySuffix : "N" / "F"
%   siteGlobal   : raw/global channel index
%   day          : session/day id in ALLMAT(:,11)
%   reason       : short free-text note

T = table( ...
    strings(0,1), zeros(0,1), zeros(0,1), strings(0,1), ...
    'VariableNames', {'monkeySuffix','siteGlobal','day','reason'});

% Nilson IT session-quality exclusions identified from the session-split
% diagnostics in March 2026. Keep the good session for each site.
T = [T; { "N", 868, 2, "Session 2 lost stimulus-evoked signal" }]; %#ok<AGROW>
T = [T; { "N", 873, 1, "Session 1 strong baseline drift / instability" }]; %#ok<AGROW>

% Optional user-managed exclusions written by the interactive session
% quality reviewer.
cfg = config();
userPath = fullfile(cfg.resultsDir, 'site_session_exclusions_user.mat');
if exist(userPath, 'file') == 2
    S = load(userPath, 'Tuser');
    if isfield(S, 'Tuser') && istable(S.Tuser) && ~isempty(S.Tuser)
        needVars = {'monkeySuffix','siteGlobal','day','reason'};
        assert(all(ismember(needVars, S.Tuser.Properties.VariableNames)), ...
            '%s must contain table Tuser with variables monkeySuffix/siteGlobal/day/reason.', userPath);
        T = [T; S.Tuser(:, needVars)]; %#ok<AGROW>
    end
end

if ~isempty(T)
    key = strcat(string(T.monkeySuffix), "_", string(T.siteGlobal), "_", string(T.day));
    [~, ia] = unique(key, 'stable');
    T = T(sort(ia), :);
end

if nargin >= 1 && ~isempty(monkey)
    monkey = string(monkey);
    T = T(T.monkeySuffix == monkey, :);
end
end
