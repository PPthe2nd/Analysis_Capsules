function tf = session_exclusion_cache_matches(S, monkeySuffix)
% SESSION_EXCLUSION_CACHE_MATCHES
% Return true when a cached analysis output was built with the same
% site-by-session exclusion list that is currently active.

current = normalize_table_local(site_session_exclusions(monkeySuffix));
cached = [];

if nargin >= 1 && isstruct(S)
    if isfield(S, 'siteSessionExclusions')
        cached = S.siteSessionExclusions;
    elseif isfield(S, 'OUT') && isstruct(S.OUT) && isfield(S.OUT, 'siteSessionExclusions')
        cached = S.OUT.siteSessionExclusions;
    end
end

if isempty(cached)
    tf = isempty(current);
    return;
end

cached = normalize_table_local(cached);
tf = isequaln(current, cached);
end

function T = normalize_table_local(T)
base = table( ...
    strings(0,1), zeros(0,1), zeros(0,1), strings(0,1), ...
    'VariableNames', {'monkeySuffix','siteGlobal','day','reason'});

if isempty(T)
    T = base;
    return;
end

assert(istable(T), 'siteSessionExclusions cache field must be a table.');
needVars = base.Properties.VariableNames;
assert(all(ismember(needVars, T.Properties.VariableNames)), ...
    'siteSessionExclusions table must contain monkeySuffix/siteGlobal/day/reason.');

T = T(:, needVars);
T.monkeySuffix = string(T.monkeySuffix);
T.siteGlobal = double(T.siteGlobal);
T.day = double(T.day);
T.reason = string(T.reason);

if ~isempty(T)
    T = sortrows(T, {'monkeySuffix','siteGlobal','day','reason'});
end
end
