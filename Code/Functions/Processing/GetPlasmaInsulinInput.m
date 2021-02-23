function IInput = GetPlasmaInsulinInput(t, P)
% Fits data using MLR to find nL and xL.
% INPUT:
%   t   - examined time point(s) [min]
%   P   - patient struct
% OUTPUT:
%   IInput - input to plasma insulin compartment at time t [mU/min]

SC = P.parameters.SC;

if P.data.IType == "human"
    if P.data.IDelivery == "intravenous"
        IInput = P.data.IBolus(t);  % [mU/min]

    elseif P.data.IDelivery == "subcutaneous"    
        if ~isfield(P.results, "QLocal")
            % Simulate SC model if not already done.
            P = SCModel(P);
        end
        
        n = GetTimeIndex(t, P.results.tArray);  % Index of current timestep.
        IInput = SC.ks3*P.results.QLocal(n);  % [mU/min]

    else
        IInput = 0;
    end

end

