function P = LineSearchOptimum(P, fieldName, searchRange, errorCriteria, funcsToApply, varargin)
% Performs a line search to find an optimal value of a parameter.
% INPUTS:
%   P             - patient struct
%   fieldName     - string name of field to alter, e.g. "P.field.subfield"
%   searchRange   - array of field values to search over
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
elseif numel(funcsToApply) == 1
    funcsToApply = {funcsToApply};
end

if isempty(varargin)
    argumentsToApply = cell(size(funcsToApply));
else
    argumentsToApply = varargin;
end


%% Setup
N = numel(searchRange);
path = split(fieldName, '.');
residualsArray = nan(size(searchRange));


message = sprintf("Line searching %s from %.3g to %.3g (N = %d)...", ...
                    fieldName, searchRange(1), searchRange(end), N);
PrintStatusUpdate(P, message);

%% Search
for ii = 1:N
    testValue = searchRange(ii);
    
    % Apply test value.
    copyP = setfield(P, path{:}, testValue);
    
    % Apply designated functions.
    for ff = 1:numel(funcsToApply)
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
P = setfield(P, path{:}, optValue);

message = sprintf("Optimal %s = %g (index %d/%d)", path{end}, optValue, iiOpt, N);
PrintStatusUpdate(P, message);

%% Plotting
plotvars.residualsArray = residualsArray;
plotvars.searchRange = searchRange;
plotvars.optValue = optValue;

P = MakePlots(P, path{end}, plotvars);

end

function P = MakePlots(P, fieldName, plotvars)
   tag = "LineSearchOptimum";

   %% Error Function
   P = AddFigure(P, tag, "ErrorFunction");
   
   plot(plotvars.searchRange, plotvars.residualsArray, 'r');
   line([plotvars.optValue plotvars.optValue], ylim());
   
   legend off
   
   xlabel(fieldName)
   ylabel("Error")
   
   legend
end

