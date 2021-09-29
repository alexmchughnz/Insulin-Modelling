function patientSet = LoadData(source, patientNums)

global CONFIG

isSource = @(dataset) contains(source, dataset, 'IgnoreCase', true);

%% Select function and numbers.

if isSource("Detemir")
    MakeDataFunction = @MakeDetemir;
    allNums = [1 3 4];
elseif isSource("DISST")
    MakeDataFunction = @MakeDISST;
    allNums = [1:50];
    bestNums = [3 5 7 8 9 13 14 16 24 25];
elseif isSource("CREBRF")
    MakeDataFunction = @MakeCREBRF;
    allNums = [146 95 68 140 12 19 147 154 33 85 126 46 156 104 72 79 ...
        73 65 78 105 138 158 87 198 128 169 186 153 115 209 196 160 145 ...
        216 166 171 220 259 240 253 235 194 263 251];
    bestNums = [12 128 146 160 166 169 171 196 198 216];
elseif isSource("OGTTLui")
    MakeDataFunction = @(pSet) MakeOGTTLui(pSet, CONFIG.ENABLELOADERPLOTS);
    allNums = [1 2 4 5 14 16 22 23 25 30];
    bestNums = [5 22 23 25 30];
end

% Replace function if mock.
if isSource("Mock")
    MakeDataFunction = @MakeMock;
end

% Replace numbers if selecting all/best.
if isequal(patientNums, "all")
    patientNums = allNums;
elseif isequal(patientNums, "best")
    patientNums = bestNums;
end


%% Set up patient set.
for ii = 1:length(patientNums)
    P.source = source;
    P.patientNum = patientNums(ii);
    P.patientCode = CONFIG.PATIENTFORMAT(P);
    patientSet{ii} = P;
end

%% Load data.
patientSet = MakeDataFunction(patientSet);

%% Add remaining elements to patient structs.
for ii = 1:length(patientSet)
    patientSet{ii} = AddParameters(patientSet{ii});
    patientSet{ii} = AddPersistents(patientSet{ii});
    
    patientSet{ii} = AddTrialInputs(patientSet{ii});
end

end

