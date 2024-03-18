function PrintResults(patientSet, recipeFunction, tag)
% Prints (and optionally saves) results of simulation.
% INPUTS:
%   patientSet - existing table of results

global CONFIG

source = patientSet{end}.source;
recipeStruct = functions(recipeFunction);
recipe = string(recipeStruct.function);


%% Tables
tables = TabulateResults(patientSet);
for tt = 1:length(tables)
    T = tables{tt};
    title = string(T.Properties.Description);
    
    disp(title);
    disp(T);
    disp(newline);
    
end

%% Plots
monitor = 3;
PanelDebugPlots(monitor);
