function P = MakePlasmaInsulinInputArray(P)
% Returns array of insulin input rate over time.
% INPUT:
%   P - patient struct
% OUTPUT:
%   P - patient struct updated with insulin input array

SC = P.parameters.SC;

if P.data.IType == "human"
    if P.data.IDelivery == "intravenous"
        IInput = P.data.IBolus(P.results.tArray);  % [mU/min]
        
    elseif P.data.IDelivery == "subcutaneous"
        if ~isfield(P.results, "QLocal")
            % Simulate SC model if not already done.
            P = SCModel(P);
        end
        
        IInput = SC.ks3*P.results.QLocal;  % [mU/min]
        
    else
        IInput = 0;
    end
    
end

P.data.IInput = IInput; % [mU/min]

end