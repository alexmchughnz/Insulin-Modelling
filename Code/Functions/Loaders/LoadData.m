function patientSet = LoadData(source, patientNums)

global CONFIG


%% Set up patient set.
for ii = 1:length(patientNums)    
    P.source = source;
    P.patientNum = patientNums(ii);
    P.patientCode = CONFIG.PATIENTFORMAT(P);
    patientSet{ii} = P;
end

%% Load data.
if source == "Detemir"
    patientSet = MakeDetemir(patientSet);
elseif source == "DISST"
    patientSet = MakeDISST(patientSet);
elseif source == "CREBRF"
    patientSet = MakeCREBRF(patientSet);
elseif source == "OGTTLui"
    patientSet = MakeOGTTLui(patientSet);
end



%% Save patient structs.
for ii = 1:length(patientSet)
    P = patientSet{ii};
    save(fullfile(CONFIG.DATAPATH, P.source, P.patientCode), '-struct', 'P');
    fprintf('%s: Saved patient data.\n', P.patientCode);
end

end

