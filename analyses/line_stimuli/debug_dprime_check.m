% debug_dprime_check
% Cross-check d' from attention_modulation_V1_3bin against an independent
% recomputation using DATA moments + Tall_V1 assignments.

clearvars;
close all;
clc;

cfg = config();

load(fullfile(cfg.matDir, 'Tall_V1_lines_N.mat'));   % Tall_V1, ALLCOORDS, RTAB384
R_resp = load_capsules_struct_exclusion_aware(fullfile(cfg.matDir, 'Resp_capsules_N_d12.mat'), 'N');
load(fullfile(cfg.matDir, 'SNR_V1_byColor_byWindow.mat')); % SNR

opts = struct();
opts.v1Sites = 1:512;
opts.timeIdx = 3;            % 300-500 ms by default
opts.excludeOverlap = true;
opts.epsDen = 1e-6;
opts.verbose = false;

OUT = attention_modulation_V1_3bin(R_resp, Tall_V1, SNR, opts);

% ---- Independent recomputation across all V1 sites ----
v1Sites = opts.v1Sites(:)';
nV1 = numel(v1Sites);
timeIdx = opts.timeIdx;
epsDen = opts.epsDen;

nTrials = R_resp.nTrials;
if isvector(nTrials)
    nTrials = nTrials(:)';
    perSiteTrials = false;
elseif ismatrix(nTrials) && size(nTrials,2) == 384
    perSiteTrials = true;
else
    error('R_resp.nTrials must be vector(384) or matrix(nSites x 384).');
end

b = SNR.muSpont(:);
topMat = [SNR.muYellowEarly(:), SNR.muYellowLate(:), SNR.muPurpleEarly(:), SNR.muPurpleLate(:)];
muTop = max(topMat, [], 2);

b = b(v1Sites);
muTop = muTop(v1Sites);
scale = muTop - b;
scale(~isfinite(scale) | scale <= 0) = NaN;

sumY_TY  = zeros(nV1,1); sumY2_TY = zeros(nV1,1); N_TY = zeros(nV1,1);
sumY_TP  = zeros(nV1,1); sumY2_TP = zeros(nV1,1); N_TP = zeros(nV1,1);
sumY_DY  = zeros(nV1,1); sumY2_DY = zeros(nV1,1); N_DY = zeros(nV1,1);
sumY_DP  = zeros(nV1,1); sumY2_DP = zeros(nV1,1); N_DP = zeros(nV1,1);

for stim = 1:384
    if ~isfield(Tall_V1(stim), 'T') || ~istable(Tall_V1(stim).T) || height(Tall_V1(stim).T) < max(v1Sites)
        continue;
    end
    Ttab = Tall_V1(stim).T;

    if ~ismember('assignment', Ttab.Properties.VariableNames) || ~ismember('center_color', Ttab.Properties.VariableNames)
        continue;
    end
    hasOverlap = ismember('overlap', Ttab.Properties.VariableNames);

    if ~perSiteTrials
        nTr_vec = repmat(nTrials(stim), nV1, 1);
    else
        nTr_vec = nTrials(v1Sites, stim);
    end

    asg = Ttab.assignment(v1Sites);
    col = Ttab.center_color(v1Sites);
    if hasOverlap
        ov = Ttab.overlap(v1Sites);
    else
        ov = zeros(nV1,1);
    end

    asgStr = toCellStrLocal(asg);
    colStr = toCellStrLocal(col);

    isBG = cellfun(@(x) strcmpi(x,'background'), asgStr);
    isT  = cellfun(@(x) strcmpi(x,'target'), asgStr);
    isD  = cellfun(@(x) strcmpi(x,'distractor'), asgStr);
    isY  = cellfun(@(x) contains(lower(x),'yellow'), colStr);
    isP  = cellfun(@(x) contains(lower(x),'purple'), colStr);

    if opts.excludeOverlap
        isOV = ov ~= 0;
    else
        isOV = false(nV1,1);
    end

    EX  = squeeze(R_resp.meanAct(v1Sites, stim, timeIdx));
    EX2 = squeeze(R_resp.meanSqAct(v1Sites, stim, timeIdx));

    EY  = (EX - b(:)) ./ scale(:);
    EY2 = (EX2 - 2.*b(:).*EX + b(:).^2) ./ (scale(:).^2);

    good = isfinite(EY) & isfinite(EY2) & isfinite(nTr_vec) & (nTr_vec > 0) & ~isBG & ~isOV & (isY | isP);

    m = good & isT & isY;
    if any(m)
        w = nTr_vec(m);
        sumY_TY(m)  = sumY_TY(m)  + w .* EY(m);
        sumY2_TY(m) = sumY2_TY(m) + w .* EY2(m);
        N_TY(m)     = N_TY(m)     + w;
    end

    m = good & isT & isP;
    if any(m)
        w = nTr_vec(m);
        sumY_TP(m)  = sumY_TP(m)  + w .* EY(m);
        sumY2_TP(m) = sumY2_TP(m) + w .* EY2(m);
        N_TP(m)     = N_TP(m)     + w;
    end

    m = good & isD & isY;
    if any(m)
        w = nTr_vec(m);
        sumY_DY(m)  = sumY_DY(m)  + w .* EY(m);
        sumY2_DY(m) = sumY2_DY(m) + w .* EY2(m);
        N_DY(m)     = N_DY(m)     + w;
    end

    m = good & isD & isP;
    if any(m)
        w = nTr_vec(m);
        sumY_DP(m)  = sumY_DP(m)  + w .* EY(m);
        sumY2_DP(m) = sumY2_DP(m) + w .* EY2(m);
        N_DP(m)     = N_DP(m)     + w;
    end
end

mu_TY = sumY_TY ./ N_TY; m2_TY = sumY2_TY ./ N_TY;
mu_TP = sumY_TP ./ N_TP; m2_TP = sumY2_TP ./ N_TP;
mu_DY = sumY_DY ./ N_DY; m2_DY = sumY2_DY ./ N_DY;
mu_DP = sumY_DP ./ N_DP; m2_DP = sumY2_DP ./ N_DP;

mu_TY(N_TY<=0)=NaN; m2_TY(N_TY<=0)=NaN;
mu_TP(N_TP<=0)=NaN; m2_TP(N_TP<=0)=NaN;
mu_DY(N_DY<=0)=NaN; m2_DY(N_DY<=0)=NaN;
mu_DP(N_DP<=0)=NaN; m2_DP(N_DP<=0)=NaN;

var_TY = max(0, m2_TY - mu_TY.^2);
var_TP = max(0, m2_TP - mu_TP.^2);
var_DY = max(0, m2_DY - mu_DY.^2);
var_DP = max(0, m2_DP - mu_DP.^2);

wY = min(N_TY, N_DY);
wP = min(N_TP, N_DP);
w = wY + wP;

muT = (wY .* mu_TY + wP .* mu_TP) ./ w;
muD = (wY .* mu_DY + wP .* mu_DP) ./ w;
m2T = (wY .* m2_TY + wP .* m2_TP) ./ w;
m2D = (wY .* m2_DY + wP .* m2_DP) ./ w;

muT(w<=0)=NaN; muD(w<=0)=NaN; m2T(w<=0)=NaN; m2D(w<=0)=NaN;

varT = max(0, m2T - muT.^2);
varD = max(0, m2D - muD.^2);

den = abs((muT + muD)/2);
den(den < epsDen) = epsDen;
idx_manual = (muT - muD) ./ den;

sdP = sqrt(0.5*(varT + varD));
sdP(sdP < epsDen) = epsDen;
dprime_manual = (muT - muD) ./ sdP;

validSite = (wY > 0) | (wP > 0);
idx_manual(~validSite) = NaN;
dprime_manual(~validSite) = NaN;

% ---- Compare ----
dDiff = dprime_manual - OUT.dprime;
iDiff = idx_manual - OUT.index;

absDDiff = abs(dDiff);
absIDiff = abs(iDiff);

fprintf('\n=== d'' CHECK SUMMARY ===\n');
fprintf('Valid sites (manual): %d / %d\n', nnz(validSite), nV1);
fprintf('max |dprime diff|: %.6g\n', max(absDDiff, [], 'omitnan'));
fprintf('median |dprime diff|: %.6g\n', median(absDDiff, 'omitnan'));
fprintf('max |index diff|: %.6g\n', max(absIDiff, [], 'omitnan'));
fprintf('median |index diff|: %.6g\n', median(absIDiff, 'omitnan'));

thr = 1e-8;
bad = find(absDDiff > thr & isfinite(absDDiff));
fprintf('Sites with |dprime diff| > %.1e: %d\n', thr, numel(bad));

if ~isempty(bad)
    [~,ord] = sort(absDDiff(bad), 'descend');
    top = bad(ord(1:min(15,end)));
    T = table(top(:), OUT.dprime(top), dprime_manual(top), dDiff(top), ...
        OUT.index(top), idx_manual(top), iDiff(top), ...
        OUT.wY(top), OUT.wP(top), ...
        'VariableNames', {'site','dprime_fun','dprime_manual','dprime_diff', ...
                          'index_fun','index_manual','index_diff','wY','wP'});
    disp(T);
else
    fprintf('No mismatches above threshold.\n');
end

% Optional focused diagnostics
sitesToInspect = [100 200 300];
sitesToInspect = sitesToInspect(sitesToInspect >= 1 & sitesToInspect <= nV1);
if ~isempty(sitesToInspect)
    fprintf('\n=== Site Diagnostics ===\n');
    for s = sitesToInspect
        fprintf(['Site %d | d'' fun=%.5f manual=%.5f | idx fun=%.5f manual=%.5f | ' ...
                 'N_TY=%.1f N_TP=%.1f N_DY=%.1f N_DP=%.1f\n'], ...
            s, OUT.dprime(s), dprime_manual(s), OUT.index(s), idx_manual(s), ...
            N_TY(s), N_TP(s), N_DY(s), N_DP(s));
    end
end

function C = toCellStrLocal(x)
if iscell(x)
    C = x;
elseif isstring(x)
    C = cellstr(x);
elseif iscategorical(x)
    C = cellstr(x);
elseif ischar(x)
    C = cellstr(x);
else
    C = cellstr(string(x));
end
C = C(:);
for i = 1:numel(C)
    if iscell(C{i}) && numel(C{i}) == 1
        C{i} = C{i}{1};
    end
end
end
