function PrintResults(patientSet, save)
% Prints (and optionally saves) results of simulation.
% INPUTS:
%   patientSet - existing table of results
%   save - logical, saves plots and table if true

global CONFIG

T = table;
for ii = 1:length(patientSet)
    P = patientSet{ii};
    T = TabulateResults(T, P);
end

if (save)
    SaveOpenFigures(P.source);
    writetable(T, fullfile(CONFIG.RESULTPATH, P.source+"table.csv"));
end

PanelDebugPlots();
disp(T);