function [tI, vI] = GetIFromITotal(P)

if P.data.IType == 'detemir'
    % Time of reading in sim [min]
    % Plasma insulin + Detemir [mU/L]
    [tITotal, vITotal] = GetSimTime(P, P.data.ITotal);
    
    if ~isfield(P.results, 'IDF')
        % Forward simulate ID model for IDF.
        P = IDModel(P);
    end
    
    % Get IDF values at ITotal measure times.
    IDF = P.results.IDF; % [mU/L]
    [~, iiI] = ismember(tITotal, P.results.tArray);
    
    % Seperate plasma I from total measurement.
    vI =  vITotal - IDF(iiI);
    tI = tITotal;
    
else    
    % Non-detemir trial.
    [tI, vI] = GetSimTime(P, P.data.I);
end

end

