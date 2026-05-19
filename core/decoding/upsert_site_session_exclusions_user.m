function Tuser = upsert_site_session_exclusions_user(monkeySuffix, siteGlobal, day, reason)
% UPSERT_SITE_SESSION_EXCLUSIONS_USER
% Add or update one user-managed site-by-session exclusion entry.

cfg = config();
userPath = fullfile(cfg.resultsDir, 'site_session_exclusions_user.mat');

Tuser = table( ...
    strings(0,1), zeros(0,1), zeros(0,1), strings(0,1), ...
    'VariableNames', {'monkeySuffix','siteGlobal','day','reason'});

if exist(userPath, 'file') == 2
    S = load(userPath, 'Tuser');
    if isfield(S, 'Tuser') && istable(S.Tuser)
        Tuser = S.Tuser(:, Tuser.Properties.VariableNames);
    end
end

monkeySuffix = string(monkeySuffix);
siteGlobal = double(siteGlobal);
day = double(day);
reason = string(reason);

isSame = (string(Tuser.monkeySuffix) == monkeySuffix) & ...
         (double(Tuser.siteGlobal) == siteGlobal) & ...
         (double(Tuser.day) == day);
Tuser(isSame, :) = [];
Tuser = [Tuser; {monkeySuffix, siteGlobal, day, reason}]; %#ok<AGROW>

save(userPath, 'Tuser', '-v7.3');
fprintf('Saved user exclusion: monkey %s | site %d | day %d -> %s\n', ...
    char(monkeySuffix), siteGlobal, day, char(reason));
fprintf('User exclusion file: %s\n', userPath);
end
