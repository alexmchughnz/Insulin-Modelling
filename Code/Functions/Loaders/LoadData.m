function patientSet = LoadData(source, patientNums, showPlots)

global CONFIG
DP = DebugPlots().LoadData;

isSource = @(dataset) contains(source, dataset, 'IgnoreCase', true);

%% Select function and numbers.

if isSource("Detemir")
    MakeDataFunction = @MakeDetemir;
    allNums = [1 3 4];
elseif isSource("DISST")
    MakeDataFunction = @MakeDISST;
    allNums = [1:50];
elseif isSource("CREBRF")
    MakeDataFunction = @MakeCREBRF;
    allNums = [146 95 68 140 12 19 147 154 33 85 126 46 156 104 72 79 ...
        73 65 78 105 138 158 87 198 128 169 186 153 115 209 196 160 145 ...
        216 166 171 220 259 240 253 235 194 263 251];
elseif isSource("OGTTLui")
    MakeDataFunction = @(pSet) MakeOGTTLui(pSet, showPlots);
    allNums = [1 2 4 5 14 16 22 23 25 30];
end

% Replace function if mock.
if isSource("Mock")
    MakeDataFunction = @MakeMock;
end

% Replace numbers if selecting all.
if isequal(patientNums, "all")
    patientNums = allNums;
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
    
    patientSet{ii} = AddBolusArrays(patientSet{ii});
    patientSet{ii} = AddPlasmaInsulinInputArray(patientSet{ii});
end

%% Produce data arrays from loaded information.
for ii = 1:length(patientSet)
end

%% Debug Plots
if DP.GlucoseInput
    for ii = 1:length(patientSet)
        P = patientSet{ii};
        MakeDebugPlot("Glucose Input", P, DP);
        
        subplot(2,1,1)
        plot(P.results.tArray, P.results.GBolus)
        ylabel("G Bolus [mmol/min]")
        
        subplot(2,1,2)
        plot(P.results.tArray, P.data.GInfusion)
        ylabel("G Infusion [mmol/min]")
    end
end

end

