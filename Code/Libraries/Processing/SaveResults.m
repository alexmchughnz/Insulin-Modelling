function SaveResults(Trial)
% Prints and saves results of simulation.
% INPUTS:
%   patientSet - existing table of results

recipeStruct = functions(Trial.recipe);
recipeName = string(recipeStruct.function);

%% Tables
tables = TabulateResults(Trial.patientSet);

for tt = 1:length(tables)
    table = tables{tt};
    title = string(table.Properties.Description);
    
    disp(title);
    disp(table);
    disp(newline);
    
    if (Trial.Config.SAVERESULTS)
        try
            filePath = fullfile(Trial.Config.RESULTPATH, Trial.outputPath, title+"_"+Trial.source+recipeName+Trial.label+".csv");
            writetable(table, filePath, "WriteRowNames", true);
            disp("Results file saved.")
        catch err
            if err.identifier == "MATLAB:table:write:FileOpenError"
                disp("Results file open - cannot save!")
            else
                rethrow(err)
            end
        end
    end
    
end
