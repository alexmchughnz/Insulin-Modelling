function P = AddPlasmaInsulinInputArray(P, forceReSim)
% Returns array of insulin input rate over time.
% INPUT:
%   P - patient struct
% OUTPUT:
%   P - patient struct updated with insulin input array


if ~exist("forceReSim", "var")
    forceReSim = false;
end

if P.data.IType == "human"
    if P.data.IDelivery == "intravenous"
        IInput = P.results.IBolus;  % [mU/min]
        
    elseif P.data.IDelivery == "subcutaneous"
        SC = P.parameters.SC;
        
        if forceReSim || ~isfield(P.results, "QLocal")
            % Simulate SC model if not already done.
            P = SCModel(P);
        end
        
        IInput = SC.ks3*P.results.QLocal;  % [mU/min]
        
    else
        IInput = 0;
    end
    
end

P.results.IInput = IInput; % [mU/min]

end