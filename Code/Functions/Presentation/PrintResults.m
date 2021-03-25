function PrintResults(patientSet, recipeFunction, tag)
% Prints (and optionally saves) results of simulation.
% INPUTS:
%   patientSet - existing table of results

global CONFIG

source = patientSet{end}.source;

recipeStruct = functions(recipeFunction);
recipe = string(recipeStruct.function);

tables = TabulateResults(patientSet);

if (CONFIG.SAVERESULTS)
    SaveOpenFigures(tag, source, recipe);
end

PanelDebugPlots();


%% Tables
for tt = 1:length(tables)
    T = tables{tt};
    title = string(T.Properties.Description);
    
    disp(title);
    disp(T);
    disp(newline);
    
    if (CONFIG.SAVERESULTS)
        T = tables{tt};
        
        sourcedir = fullfile(CONFIG.RESULTPATH, source);
        if ~exist(sourcedir, 'dir')
            mkdir(sourcedir);
        end
        
        recipedir = fullfile(sourcedir, recipe);
        if ~exist(recipedir, 'dir')
            mkdir(recipedir);
        end
        
        filepath = fullfile(recipedir, tag+recipe+title+".csv");
        writetable(T, filepath, "WriteRowNames", true);
    end
end