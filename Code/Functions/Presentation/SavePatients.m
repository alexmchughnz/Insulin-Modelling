function SavePatients(patientSet, location)

global CONFIG
if ~exist("location", "var")
    location = CONFIG.RESULTPATH;
end

for ii = 1:length(patientSet)
    if length(patientSet) == 1
        P = patientSet;
    else
        P = patientSet{ii};
    end
    
    code = MakeValidName(P.patientCode);
    
    patientDir = fullfile(location, P.source);
    if ~isfolder(patientDir)
        mkdir(patientDir);
    end
    save(fullfile(patientDir, code), '-struct', 'P');
    
    PrintStatusUpdate(P, "Saved patient.")
end

end

