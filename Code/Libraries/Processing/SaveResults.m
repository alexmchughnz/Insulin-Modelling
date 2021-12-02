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
            filePath = fullfile(Trial.Config.RESULTPATH, Trial.outputPath, title+"_"+Trial.source+recipeName+fileLabel+".csv");
            writetable(table, filePath, "WriteRowNames", true);
        catch err
            
            disp("Results file open - cannot save!")
        end
    end
end

end
