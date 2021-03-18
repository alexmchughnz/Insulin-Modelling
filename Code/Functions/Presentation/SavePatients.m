function SavePatients(patientSet)

global CONFIG

for ii = 1:length(patientSet)
    P = patientSet{ii};
    code = matlab.lang.makeValidName(P.patientCode);     
    
    patientDir = fullfile(CONFIG.RESULTPATH, P.source);
    if ~isfolder(patientDir)
        mkdir(patientDir);
    end    
    save(fullfile(patientDir, code), '-struct', 'P');
    
    PrintStatusUpdate(P, "Saved patient.")
end

end

