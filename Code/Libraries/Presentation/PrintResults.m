function PrintResults(Trial)
% Prints (and optionally saves) results of simulation.
% INPUTS:
%   patientSet - existing table of results

global CONFIG


recipeStruct = functions(Trial.recipe);
recipeName = string(recipeStruct.function);

%% Tables
tables = TabulateResults(Trial.resultSet);

for tt = 1:length(tables)
    table = tables{tt};
    title = string(table.Properties.Description);
    
    disp(title);
    disp(table);
    disp(newline);
    
    if (CONFIG.SAVERESULTS)
        trialPath = fullfile(Trial.source, recipeName);
        
        % Append label if present.
        fileLabel = "";
        if Trial.label ~= ""
            trialPath = fullfile(trialPath, Trial.label);
            fileLabel = "-"+Trial.label;
        end
        
        try
            filePath = fullfile(CONFIG.RESULTPATH, trialPath, title+"_"+Trial.source+recipeName+fileLabel+".csv");
            writetable(table, filePath, "WriteRowNames", true);
        catch
            disp("Results file open - cannot save!")
        end
    end
end

%% Plots
if (CONFIG.SAVERESULTS)
    SaveOpenFigures(Trial, trialPath);
end

monitor = 3;
PanelDebugPlots(monitor);
