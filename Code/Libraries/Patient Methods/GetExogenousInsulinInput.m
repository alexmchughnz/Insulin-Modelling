function [Uex, P] = GetExogenousInsulinInput(P, forceResim)
% Returns array of input rate to I compartment over time.
% 
% INPUT:
%   P - patient struct
%   forceResim - true if SCModel() needs to be run again to update QLocal.
% OUTPUT:
%   P - patient struct updated with insulin input array


if ~exist("forceResim", "var")
    forceResim = false;
end

if P.data.IType == "human"
    if P.data.IDelivery == "intravenous"
        Uex = P.results.IBolus;  % [mU/min]
        
    elseif P.data.IDelivery == "subcutaneous"
        if forceResim || ~isfield(P.results, "QLocal")
            % Simulate SC model if not already done.
            % This may need to be re-run if any SC model parameters (e.g.
            % IInput or rate parameters) are changed.
            P = SCModel(P);
        end
        
        Uex = P.parameters.SC.ks3*P.results.QLocal;  % [mU/min]
    end
    
else
    Uex = zeros(size(P.results.tArray));
end

end