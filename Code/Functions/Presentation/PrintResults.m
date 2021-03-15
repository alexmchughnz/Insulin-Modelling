function PrintResults(patientSet, recipeFunction)
% Prints (and optionally saves) results of simulation.
% INPUTS:
%   patientSet - existing table of results

global CONFIG

source = patientSet{end}.source;

recipeStruct = functions(recipeFunction);
recipe = recipeStruct.function;

tables = TabulateResults(patientSet);

if (CONFIG.SAVERESULTS)
    SaveOpenFigures(source, recipe);
    
    for tt = 1:length(tables)
        T = tables{tt};
        
        title = T.Properties.Description;
        filename = fullfile(CONFIG.RESULTPATH, source + recipe + title + ".csv");
        
        writetable(T, filename, ...
                "WriteRowNames", true);
    end
end

PanelDebugPlots();

for tt = 1:length(tables)
    T = tables{tt};
    
    disp(T.Properties.Description);
    disp(T);
    disp(newline);
end