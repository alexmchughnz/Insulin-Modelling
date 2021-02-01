function SavePatients(patientSet)

global CONFIG

for ii = 1:length(patientSet)
    P = patientSet{ii};
    save(fullfile(CONFIG.DATAPATH, P.source, P.patientCode), '-struct', 'P');
    PrintStatusUpdate(P, "Saved patient.")
end

end

