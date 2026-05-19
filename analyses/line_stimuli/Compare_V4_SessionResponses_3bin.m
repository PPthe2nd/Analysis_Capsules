function OUT = Compare_V4_SessionResponses_3bin(Puser)
% COMPARE_V4_SESSIONRESPONSES_3BIN
% Compare day-1 vs day-2 V4 responses using the standard coarse 3 windows
% from the existing SNR_capsules_*_d12.mat pipeline.

if nargin < 1
    Puser = struct();
end

P = struct();
P.Monkey = 1;          % 1 = Nilson, 2 = Figaro
P.SNRthr = 0.7;
P.MinObjectStim = 1;
P.onlyCorrect = true;
P.correctCol = 9;
P.correctVal = 1;
P.dayCol = 11;
P.sessions = [1 2];
P.expectedMaxStim = 384;
P.chunkTrials = 200;
P.useCache = true;
P.saveResult = true;
P.plotFigure = true;

if ~isempty(Puser)
    fn = fieldnames(Puser);
    for i = 1:numel(fn)
        P.(fn{i}) = Puser.(fn{i});
    end
end

OUT = compare_area_session_responses_3bin('V4', P.Monkey, P);
end
