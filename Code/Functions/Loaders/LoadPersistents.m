function P = LoadPersistents(P)

global CONFIG

filename = fullfile(CONFIG.RESULTPATH, P.source, P.patientCode+'.mat');

if isfile(filename)
    loadP = load(fullfile(CONFIG.RESULTPATH, P.source, P.patientCode));
    
    if isfield(loadP, "persistents")
        P.persistents = loadP.persistents;
    end
    
end

end
