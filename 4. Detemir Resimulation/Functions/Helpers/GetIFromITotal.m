function [tI, vI] = GetIFromITotal(P)

if isfield(P.data, 'ITotal')    
    % Time of reading in sim [min]
    % Plasma insulin + Detemir [mU/L]
    [tITotal, vITotal] = GetSimTime(P, P.data.ITotal);
    
    % Forward simulate ID model for IDF.
    P = IDModel(P);
    IDF = P.results.IDF; % [mU/L]
    
    % Get IDF values at ITotal measure times.
    [~, iiI] = ismember(tITotal, P.results.tArray);
    
    % Seperate plasma I from total measurement.
    vI =  vITotal - IDF(iiI);
    tI = tITotal;
else
    % Non-detemir trial.
    [tI, vI] = GetSimTime(P, P.data.I);
end

end

