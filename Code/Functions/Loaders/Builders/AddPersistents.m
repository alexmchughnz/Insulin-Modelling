function P = AddPersistents(P)

global CONFIG

code = matlab.lang.makeValidName(P.patientCode);

filename = fullfile(CONFIG.RESULTPATH, P.source, code+'.mat');

if isfile(filename)
    loadP = load(fullfile(CONFIG.RESULTPATH, P.source, code));
    
    if isfield(loadP, "persistents")
        fields = fieldnames(loadP.persistents);
        for ii = 1:numel(fields)
            P.persistents.(fields{ii}) = loadP.persistents.(fields{ii});
        end
    end
    
end

end
