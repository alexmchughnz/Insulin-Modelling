function PrintResults(patientSet, save)
% Prints (and optionally saves) results of simulation.
% INPUTS:
%   patientSet - existing table of results
%   save - logical, saves plots and table if true

global CONFIG

source = patientSet{end}.source;
tables = TabulateResults(patientSet);

if (save)
    SaveOpenFigures(source);
    
    for tt = 1:length(tables)
        T = tables{tt};
        
        title = T.Properties.Description;
        filename = fullfile(CONFIG.RESULTPATH, source + title + ".csv");
        
        writetable(T, filename);
    end
end

PanelDebugPlots();

for tt = 1:length(tables)
    T = tables{tt};
    
    disp(T.Properties.Description);
    disp(T);
    disp(newline);
end