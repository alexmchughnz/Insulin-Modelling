function patientSet = LoadData(source, patientNums)

global CONFIG DEBUGPLOTS

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

%% Apply parameters and persistents to patient structs.
for ii = 1:length(patientSet)
    patientSet{ii} = LoadParameters(patientSet{ii});
    patientSet{ii} = LoadPersistents(patientSet{ii});
end

%% Debug Plots
DP = DEBUGPLOTS.LoadData;
if DP.GlucoseInput
    for ii = 1:length(patientSet)
        P = patientSet{ii};
        MakeDebugPlot("Glucose Input", P, DP);

        subplot(2,1,1)
        plot(P.results.tArray, P.data.GBolus(P.results.tArray))
        ylabel("G Bolus [mmol/min]")

        subplot(2,1,2)
        plot(P.results.tArray, P.data.GInfusion)
        ylabel("G Infusion [mmol/min]")
    end
end

end

