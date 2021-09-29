function P = LineSearchOptimum(P, fieldPath, searchRange, errorCriteria, funcsToApply, varargin)
% Performs a line search to find an optimal value of a parameter.
% INPUTS:
%   P             - patient struct
%   fieldPath     - string path of field to alter, e.g. "P.field.subfield"
%   searchLine    - array of field values to search over
%   errorCriteria - function(P) yielding error at each field value
%   funcsToApply  - cell array of function handles to apply to each struct
%                   BEFORE running error criteria
%   varargin      - varargin{n} contains any and all arguments for
%                   funcsToApply{n}
%
% OUTPUT:
%   P   - modified patient struct with optimal field value with respect to
%         errorCriteria


% Pre-formatting
if ~exist("funcsToApply", "var")
    funcsToApply = {};
elseif length(funcsToApply) == 1
    funcsToApply = {funcsToApply};
end

if isempty(varargin)
    argumentsToApply = cell(size(funcsToApply));
else
    argumentsToApply = varargin;
end


%% Setup
N = numel(searchRange);
fieldSteps = split(fieldPath, '.');
residualsArray = nan(size(searchRange));


message = sprintf("Line searching %s from %.3g to %.3g (N = %d)...", ...
                    fieldPath, searchRange(1), searchRange(end), N);
PrintStatusUpdate(P, message);

%% Search
for ii = 1:N
    testValue = searchRange(ii);
    
    % Apply test value.
    copyP = setfield(P, fieldSteps{:}, testValue);
    
    % Apply designated functions.
    for ff = 1:length(funcsToApply)
        func = funcsToApply{ff};
        args = argumentsToApply{ff};
        
        if ~isempty(args)
            copyP = func(copyP, args(:));
        else
            copyP = func(copyP);
        end        
    end
    
    % Run objective function to get error for this test value.
    residualsArray(ii) = errorCriteria(copyP);
end

%% Saving
% Get optimal value of field.
[~, iiOpt] = min(residualsArray);
optValue = searchRange(iiOpt);

% Write to struct.
P = setfield(P, fieldSteps{:}, optValue);

message = sprintf("Optimal %s = %g (index %d/%d)", fieldSteps{end}, optValue, iiOpt, N);
PrintStatusUpdate(P, message);

%% Plotting
plotvars.residualsArray = residualsArray;
plotvars.searchRange = searchRange;
plotvars.optValue = optValue;
MakePlots(P, fieldSteps{end}, plotvars);

end

function MakePlots(P, fieldName, plotvars)
DP = DebugPlots().LineSearchOptimum;

if DP.ErrorFunction
   MakeDebugPlot("Line Search "+fieldName, P, DP);
   
   plot(plotvars.searchRange, plotvars.residualsArray, 'r');
   line([plotvars.optValue plotvars.optValue], ylim());
   
   legend off
   
   xlabel(fieldName)
   ylabel("Error")
   
   legend
end

end

