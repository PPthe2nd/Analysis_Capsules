function OUT = Compare_FrontBackCentroid_Tuning_IT()
% COMPARE_FRONTBACKCENTROID_TUNING_IT
% Compare IT variance explained for front/attended and back/occluded
% capsule centroid tuning.

%% Settings
P = struct();
P.minUsefulVE = 10;      % minimum max(VE_front, VE_back) to count a site as useful
P.minVEDiff = 5;         % minimum VE difference to call one centroid better
P.makeFigures = true;

OUTfront = AttendedCentroid_Tuning_IT;
OUTback = OccludedCentroid_Tuning_IT;

assert(isequal(OUTfront.RFrange(:), OUTback.RFrange(:)), ...
    'Front/back centroid analyses must use the same IT site set.');
assert(strcmp(char(OUTfront.monkeySuffix), char(OUTback.monkeySuffix)), ...
    'Front/back centroid analyses must use the same monkey.');

veFrontEarly = [OUTfront.FitEarly.r2TrainPct].';
veFrontLate = [OUTfront.FitLate.r2TrainPct].';
veBackEarly = [OUTback.FitEarly.r2TrainPct].';
veBackLate = [OUTback.FitLate.r2TrainPct].';

diffEarly = veFrontEarly - veBackEarly;
diffLate = veFrontLate - veBackLate;

useEarly = isfinite(veFrontEarly) & isfinite(veBackEarly) & ...
    (max(veFrontEarly, veBackEarly) > P.minUsefulVE);
useLate = isfinite(veFrontLate) & isfinite(veBackLate) & ...
    (max(veFrontLate, veBackLate) > P.minUsefulVE);

frontBetterEarly = useEarly & (diffEarly > P.minVEDiff);
backBetterEarly = useEarly & (diffEarly < -P.minVEDiff);
similarEarly = useEarly & (abs(diffEarly) <= P.minVEDiff);

frontBetterLate = useLate & (diffLate > P.minVEDiff);
backBetterLate = useLate & (diffLate < -P.minVEDiff);
similarLate = useLate & (abs(diffLate) <= P.minVEDiff);

fprintf('Front/back centroid comparison for monkey %s\n', char(OUTfront.monkeySuffix));
fprintf('Useful sites (early): %d | front better: %d | back better: %d | similar: %d\n', ...
    nnz(useEarly), nnz(frontBetterEarly), nnz(backBetterEarly), nnz(similarEarly));
fprintf('Useful sites (late):  %d | front better: %d | back better: %d | similar: %d\n', ...
    nnz(useLate), nnz(frontBetterLate), nnz(backBetterLate), nnz(similarLate));

if P.makeFigures
    limEarly = compute_square_lim([veFrontEarly(useEarly); veBackEarly(useEarly)]);
    limLate = compute_square_lim([veFrontLate(useLate); veBackLate(useLate)]);

    figure('Color', 'w');
    subplot(1,2,1);
    scatter(veFrontEarly(useEarly), veBackEarly(useEarly), 30, 'filled');
    hold on;
    plot([0 limEarly], [0 limEarly], 'k-');
    xlim([0 limEarly]);
    ylim([0 limEarly]);
    axis square;
    xlabel('Front centroid VE (%)');
    ylabel('Back centroid VE (%)');
    title(sprintf('Early | useful=%d', nnz(useEarly)));
    grid on;

    subplot(1,2,2);
    scatter(veFrontLate(useLate), veBackLate(useLate), 30, 'filled');
    hold on;
    plot([0 limLate], [0 limLate], 'k-');
    xlim([0 limLate]);
    ylim([0 limLate]);
    axis square;
    xlabel('Front centroid VE (%)');
    ylabel('Back centroid VE (%)');
    title(sprintf('Late | useful=%d', nnz(useLate)));
    grid on;

    figure('Color', 'w');
    histogram(diffEarly(useEarly), 30, 'FaceColor', [0.25 0.45 0.85], 'EdgeColor', 'none');
    hold on;
    histogram(diffLate(useLate), 30, 'FaceColor', [0.85 0.35 0.25], 'EdgeColor', 'none');
    xline(0, 'k-');
    xline(P.minVEDiff, 'k--');
    xline(-P.minVEDiff, 'k--');
    xlabel('VE_{front} - VE_{back} (%)');
    ylabel('N sites');
    title(sprintf('Front vs back centroid VE difference | useful>%g%%', P.minUsefulVE));
    legend('Early', 'Late');
    grid on;
end

Tlate = table((1:numel(OUTfront.RFrange)).', OUTfront.RFrange(:), ...
    veFrontLate, OUTfront.pValueApproxLate(:), veBackLate, OUTback.pValueApproxLate(:), ...
    'VariableNames', {'localSite','globalSiteInR','veFrontLate','pFrontLate','veBackLate','pBackLate'});
Tlate.veDiffLate = Tlate.veFrontLate - Tlate.veBackLate;
Tlate = sortrows(Tlate, 'veDiffLate', 'descend');

disp('Top IT sites better explained by front centroid in the late window:');
disp(Tlate(1:min(10,height(Tlate)), :));
disp('Top IT sites better explained by back centroid in the late window:');
disp(flipud(Tlate(max(1,height(Tlate)-9):height(Tlate), :)));

OUT.P = P;
OUT.monkeySuffix = OUTfront.monkeySuffix;
OUT.RFrange = OUTfront.RFrange;
OUT.veFrontEarly = veFrontEarly;
OUT.veFrontLate = veFrontLate;
OUT.veBackEarly = veBackEarly;
OUT.veBackLate = veBackLate;
OUT.diffEarly = diffEarly;
OUT.diffLate = diffLate;
OUT.useEarly = useEarly;
OUT.useLate = useLate;
OUT.frontBetterEarly = frontBetterEarly;
OUT.backBetterEarly = backBetterEarly;
OUT.similarEarly = similarEarly;
OUT.frontBetterLate = frontBetterLate;
OUT.backBetterLate = backBetterLate;
OUT.similarLate = similarLate;
OUT.TableLateCompare = Tlate;
end

function lim = compute_square_lim(v)
if isempty(v) || ~any(isfinite(v))
    lim = 1;
    return;
end
lim = max(v(isfinite(v)));
if ~isfinite(lim) || lim <= 0
    lim = 1;
end
lim = 1.05 * lim;
end
