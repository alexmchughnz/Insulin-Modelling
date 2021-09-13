function P = FindOptimalValue(P, fieldPath, searchGrid, objFunc, applyFuncs)
% Performs a grid search to find an optimal value of a parameter.
% INPUTS:
%   P      - patient struct
% OUTPUT:
%   P   - modified patient struct with optimal parameter value.   

if ~exist("applyFuncs", "var")
    applyFuncs = {};
elseif length(applyFuncs) == 1
    applyFuncs = {applyFuncs};
end


%% Setup
fieldSteps = split(fieldPath, '.');
N = numel(searchGrid);
objGrid = nan(size(searchGrid));

%% Search
for ii = 1:N   
    testValue = searchGrid(ii);
    
    % Apply test value.
    copyP = setfield(P, fieldSteps{:}, testValue);
    
    % Apply designated functions.
    for ff = 1:length(applyFuncs)
        func = applyFuncs{ff};
        copyP = func(copyP);
    end
    
    % Run objective function to get error for this test value.
    objGrid(ii) = objFunc(copyP);
end

%% Saving
% Get optimal value of field.
[~, iiOpt] = min(objGrid);
optValue = searchGrid(iiOpt);

% Write to struct.
P = setfield(P, fieldSteps{:}, optValue);
        
message = sprintf("Optimal %s = %g (index %d/%d)", fieldSteps{end}, optValue, iiOpt, N);
PrintStatusUpdate(P, message);     


end




