function P = FixAndFitParameters(P, integralSystemFunc, fixedParams)
% Fits data using MLR to fit for parameters.
% INPUT:
%   P   - patient struct
%   forcenLxL - [nL, xL] to force (for plots)
% OUTPUT:
%   P   - modified patient struct with nL and xL


if ~exist("fixedParams", "var") || all(isnan(fixedParams))
    % No fixed parameters.
    P = IntegralFitParameters(P, integralSystemFunc);
    
else
    [A, b, ~, ~, paramNames] = integralSystemFunc(P);
    
    if all(~isnan(fixedParams))
        % Two fixed parameters.
        param1 = fixedParams(1);
        param2 = fixedParams(2);
        
    else
        % One fixed parameter.
        if ~isnan(fixedParams(1))
            % First parameter fixed.
            param1 = fixedParams(1);
            param2 = A(:,2) \ (b - A(:,1)*param1);
            
        elseif ~isnan(fixedParams(2))
            % Second parameter fixed.
            param2 = fixedParams(2);
            param1 = A(:,1) \ (b - A(:,2)*param2);
        end
    end
    
    % Save parameters.
    P.results.(paramNames(1)) = param1;
    P.results.(paramNames(2)) = param2;
end

end