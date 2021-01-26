function patientSet = LoadData(source, patientNums)

global CONFIG

%% Load data.
if source == "Detemir"
    patientSet = MakeDetemir(patientNums);
end



%% Save patient structs.
for ii = 1:length(patientSet)
    P = patientSet{ii};
    filename = CONFIG.PATIENTFORMAT(P.source, P.patientNum);
    save(fullfile(CONFIG.DATAPATH, P.source, filename), '-struct', 'P');
    fprintf('%s: Saved patient data.\n', P.patientCode);
end

end

