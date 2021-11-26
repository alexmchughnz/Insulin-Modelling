function SavePatients(Trial, patientSet, location)

if ~exist("location", "var")
    location = Trial.Config.RESULTPATH;
end

if length(patientSet) == 1
    patientSet = {patientSet};
end


% Extract recipe name.
recipeStruct = functions(Trial.recipe);
recipeName = string(recipeStruct.function);

% Save each patient struct in the correct Trial > Recipe > Label directory.
for ii = 1:length(patientSet)
    P = patientSet{ii};
    name = MakeValidName(Trial.Config.PATIENTFILEFORMAT(Trial, P));
    
    trialPath = fullfile(Trial.source, recipeName);
    
    % Append label if present.
    fileLabel = "";
    if Trial.label ~= ""
        trialPath = fullfile(trialPath, Trial.label);
        fileLabel = "-"+Trial.label;
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

