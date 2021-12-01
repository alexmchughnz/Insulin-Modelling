function SavePatients(Trial, patientSet)

if length(patientSet) == 1
    patientSet = {patientSet};
end

% Save each patient struct in the correct Trial > Recipe > Label directory.
for ii = 1:length(patientSet)
    P = patientSet{ii};
    name = MakeValidName(Trial.Config.PATIENTFILEFORMAT(Trial, P));
    
    resultsDir = fullfile(Trial.Config.RESULTPATH, Trial.outputPath);
    if ~isfolder(resultsDir)
        mkdir(resultsDir);
    end

    % Remove big figure handles from saved structs.
    saveP = P;
    saveP = rmfield(saveP, "figures");

    save(fullfile(resultsDir, name), '-struct', 'saveP');
    
    message = sprintf("Saved patient %s.", name);
    PrintStatusUpdate(P, message);
end

end

