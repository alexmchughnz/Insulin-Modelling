function PrintResults(T)
% Prints (and optionally saves) results of simulation.
% INPUTS:
%   patientSet - existing table of results

global CONFIG


recipeStruct = functions(T.recipe);
recipeName = string(recipeStruct.function);

%% Tables
tables = TabulateResults(T.resultSet);
for tt = 1:length(tables)
    table = tables{tt};
    title = string(table.Properties.Description);
    
    disp(title);
    disp(table);
    disp(newline);
    
    if (CONFIG.SAVERESULTS)
        trialPath = fullfile(T.source, recipeName);
        
        % Append label if present.
        if ~isempty(T.label)
            trialPath = fullfile(trialPath, T.label);
        end
        
        try
            filepath = fullfile(CONFIG.RESULTPATH, trialPath, title+".csv");
            writetable(table, filepath, "WriteRowNames", true);
        catch
            disp("Results file open - cannot save!")
        end
    end
end

%% Plots
if (CONFIG.SAVERESULTS)
    SaveOpenFigures(T, trialPath);
end

monitor = 3;
PanelDebugPlots(monitor);
