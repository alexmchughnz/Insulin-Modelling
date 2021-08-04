function [P, hasPersistent] = GetPersistent(P, name)
    
    hasPersistent = false;  
    
    % P may be a modified version of some original P.
    % Reload persistents based on patient code and try again.
%     P = AddPersistents(P);  
    
    % Check if this struct has persistent saved.
    if isfield(P, "persistents")
        if isfield(P.persistents, name)
            hasPersistent = true;
        end        
    end
    
end