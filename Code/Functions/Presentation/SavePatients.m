function SavePatients(patientSet)

global CONFIG

for ii = 1:length(patientSet)
    P = patientSet{ii};
    code = matlab.lang.makeValidName(P.patientCode);     
    
    save(fullfile(CONFIG.DATAPATH, P.source, P.code), '-struct', 'P');
    PrintStatusUpdate(P, "Saved patient.")
end

end

