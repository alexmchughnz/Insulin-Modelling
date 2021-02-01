function P = LoadPersistents(P)

global CONFIG

filename = fullfile(CONFIG.DATAPATH, P.source, P.patientCode);

if isfile(filename)
    loadP = load(fullfile(CONFIG.DATAPATH, P.source, P.patientCode));
    
    if isfield(loadP, "persistents")
        P.persistents = loadP.persistents;
    end
    
end

end
