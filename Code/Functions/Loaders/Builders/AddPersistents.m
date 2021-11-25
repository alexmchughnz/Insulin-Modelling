function P = AddPersistents(Trial, P)

global CONFIG

code = MakeValidName(P.patientCode);

filename = fullfile(CONFIG.RESULTPATH, Trial.source, code+'.mat');

if isfile(filename)
    loadP = load(fullfile(CONFIG.RESULTPATH, Trial.source, code));
    
    if isfield(loadP, "persistents")
        fields = fieldnames(loadP.persistents);
        for ii = 1:numel(fields)
            P.persistents.(fields{ii}) = loadP.persistents.(fields{ii});
        end
    end
    
end

end
