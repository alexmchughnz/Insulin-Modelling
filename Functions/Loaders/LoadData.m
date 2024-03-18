function patientSet = LoadData(source, patientNums)

global CONFIG

isSource = @(dataset) contains(source, dataset, 'IgnoreCase', true);


%% Set up patient set.
for ii = 1:length(patientNums)
    P.source = source;
    P.patientNum = patientNums(ii);
    P.patientCode = CONFIG.PATIENTFORMAT(P);
    P.patientSuffix = "";
    patientSet{ii} = P;
end

%% Load data.
patientSet = MakeTemplate(patientSet);

%% Add remaining elements to patient structs.
for ii = 1:length(patientSet)
    patientSet{ii} = AddParameters(patientSet{ii});
    patientSet{ii} = AddPersistents(patientSet{ii});
    
    patientSet{ii} = AddBolusArrays(patientSet{ii});
    patientSet{ii} = AddPlasmaInsulinInputArray(patientSet{ii});
end

end

