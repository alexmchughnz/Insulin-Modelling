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
        fileLabel = "";
        if T.label ~= ""
            trialPath = fullfile(trialPath, T.label);
            fileLabel = "-"+T.label;
        end
        
        try
            filePath = fullfile(CONFIG.RESULTPATH, trialPath, title+"_"+T.source+recipeName+fileLabel+".csv");
            writetable(table, filePath, "WriteRowNames", true);
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
