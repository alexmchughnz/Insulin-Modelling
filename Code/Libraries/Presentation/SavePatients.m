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
    code = MakeValidName(P.patientCode);    
    
    if isempty(T.label)
        trialDir = fullfile(location, T.source, recipeName);
    else        
        trialDir = fullfile(location, T.source, recipeName, T.label);
    end    
    
    if ~isfolder(trialDir)
        mkdir(trialDir);
    end
    save(fullfile(trialDir, code), '-struct', 'P');
    
    message = sprintf("Saved patient %s.", code);
    PrintStatusUpdate(P, message);
end

end

