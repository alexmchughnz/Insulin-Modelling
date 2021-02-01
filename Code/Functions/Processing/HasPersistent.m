function hasPersistent = HasPersistent(P, name)
    
    hasPersistent = false;
    
    if isfield(P, "persistents")
        if isfield(P.persistents, name)
            hasPersistent = true;
        end        
    end
end