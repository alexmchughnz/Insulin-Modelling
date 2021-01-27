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

%% Save patient structs.
for ii = 1:length(patientSet)
    P = patientSet{ii};
    save(fullfile(CONFIG.DATAPATH, P.source, P.patientCode), '-struct', 'P');
    fprintf('%s: Saved patient data.\n', P.patientCode);
end

%% Debug Plots
DP = DEBUGPLOTS.LoadData;
if DP.GlucoseInput
    for ii = 1:length(patientSet)
        P = patientSet{ii};
        MakeDebugPlot(P, DP);

        GI = zeros(size(P.results.tArray));
        GB = zeros(size(P.results.tArray));
        for tt = 1:length(P.results.tArray)
            time = P.results.tArray(tt);
            GI(tt) = GetGlucoseDelivery(time, P);
            GB(tt) = P.data.GBolus(time);
        end

        subplot(3,1,1)
        plot(P.results.tArray, GI)
        title(sprintf("P%s: G Input", P.patientCode))
        ylabel("[mmol/min]")

        subplot(3,1,2)
        plot(P.results.tArray, GB)
        title(sprintf("P%s: G Bolus", P.patientCode))
        ylabel("[mmol/min]")

        subplot(3,1,3)
        plot(P.results.tArray, P.data.GInfusion)
        title(sprintf("P%s: G Infusion", P.patientCode))
        ylabel("[mmol/min]")
    end
end

end

