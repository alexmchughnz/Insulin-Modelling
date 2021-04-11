function SavePatients(patientSet, location)

global CONFIG
if ~exist("location", "var")
    location = CONFIG.RESULTPATH;
end


for ii = 1:length(patientSet)
    P = patientSet{ii};
    code = matlab.lang.makeValidName(P.patientCode);     
    
    patientDir = fullfile(location, P.source);
    if ~isfolder(patientDir)
        mkdir(patientDir);
    end    
    save(fullfile(patientDir, code), '-struct', 'P');
    
    PrintStatusUpdate(P, "Saved patient.")
end

end

 