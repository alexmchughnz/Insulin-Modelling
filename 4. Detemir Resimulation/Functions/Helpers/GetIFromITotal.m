function [tI, vI] = GetIFromITotal(P)
% NOTE: takes data in pmol, and returns in mU!!
global C

if isfield(P.data, 'ITotal')    
    % Time of reading in sim [min]
    % Plasma insulin + Detemir [pmol/L]
    [tITotal, vITotal] = GetSimTime(P, P.data.ITotal);
    vITotal = C.pmol2mU(vITotal);
    
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
    vI = C.pmol2mU(vI);
end

end

