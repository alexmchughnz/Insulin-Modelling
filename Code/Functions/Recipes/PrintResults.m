function PrintResults(patientSet, save)
% Recipe to print (and optionally save) results of simulation.
% INPUTS:
%   patientSet - existing table of results
%   save - logical, saves plots and table if true

global CONFIG

T = table;
for ii = 1:length(patientSet)
    P = patientSet{ii};
    
    PlotModelResults(P);
    T = TabulateResults(T, P);
end

if (save)
    SaveOpenFigures(P.source);
    writetable(T, fullfile(CONFIG.RESULTPATH, P.source+"table.csv"));
end

PanelDebugPlots(1);
disp(T);