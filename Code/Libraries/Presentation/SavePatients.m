function SavePatients(T, patientSet, location)

global CONFIG

if ~exist("location", "var")
    location = CONFIG.RESULTPATH;
end

if length(patientSet) == 1
    patientSet = {patientSet};
end


% Extract recipe name.
recipeStruct = functions(T.recipe);
recipeName = string(recipeStruct.function);

% Save each patient struct in the correct Trial > Recipe > Label directory.
for ii = 1:length(patientSet)
    P = patientSet{ii};
    name = MakeValidName(CONFIG.PATIENTFILEFORMAT(T, P));
    
    trialPath = fullfile(T.source, recipeName);
    
    % Append label if present.
    fileLabel = "";
    if T.label ~= ""
        trialPath = fullfile(trialPath, T.label);
        fileLabel = "-"+T.label;
    end
    
    trialDir = fullfile(location, trialPath);
    if ~isfolder(trialDir)
        mkdir(trialDir);
    end
    save(fullfile(trialDir, name+"_"+recipeName+fileLabel), '-struct', 'P');
    
    message = sprintf("Saved patient %s.", name);
    PrintStatusUpdate(P, message);
end

end

