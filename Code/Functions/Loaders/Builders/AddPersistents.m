function P = AddPersistents(P)

global CONFIG

code = MakeValidName(P.patientCode);

filename = fullfile(CONFIG.RESULTPATH, P.source, code+'.mat');

if isfile(filename)
    loadP = load(fullfile(CONFIG.RESULTPATH, P.source, code));
    
    if isfield(loadP, "persistents")
        P.persistents = loadP.persistents;
    end
    
end

end
