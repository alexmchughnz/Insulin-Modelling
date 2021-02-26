function PArray = AdjustSimulationUen(P)
% Recipe for adjusting Uen and counter balancing by changing inputs.
% INPUTS:
%   P  - patient struct
% OUTPUT:
%   PArray  - updated patient structs

%% Plots
plots = DebugPlots();

plots.SolveSystem.Glucose = 1;
plots.SolveSystem.Insulin = 1;

DebugPlots(plots);


%% Variables
UenProportion = 1.15;

%% Setup
% Simulate patient as is.
P = EstimateInsulinSecretion(P);
P = FitHepaticClearance(P);
P = FindGutEmptyingRate(P);
P = FitInsulinSensitivity(P, false);
targetSI = P.results.SI;

% Create a copy of the original P, with the minimum req. variables.
adjP = rmfield(P, 'results');
adjP.patientCode = P.patientCode+"-adj";
fieldsToKeep = ["tArray", "Uen", "nL", "xL", "d2", "SI"];
for ii = 1:length(fieldsToKeep)
    field = fieldsToKeep(ii);
    adjP.results.(field) = P.results.(field);
end

%% Functions
% Simulate adjP with adjusted Uen to get result SI.
adjP.results.Uen = P.results.Uen .* UenProportion;

% Find IInput value that minimises error in SI.
function error = GetSIError(IProportion)
    % Make new copy of patient and adjust IInput and bolus functions.
    copyP = adjP;
    copyP.data.vIBolus = adjP.data.vIBolus .* IProportion;
    copyP = MakeBolusFunctions(copyP);
    
    % Fit SI and store results.
    copyP = FitInsulinSensitivity(copyP, false);
    newSI = copyP.results.SI;    
    error = abs(newSI-targetSI);
end
IProportion = fminbnd(@GetSIError, 0.5, 1.5);


% Save IProportions.
P.results.IInputProportion = 1.0;
adjP.results.IInputProportion = IProportion;

% Just finishing simulation for plots!
P = SolveSystem(P, true);
adjP = SolveSystem(adjP, true);

PArray = {P adjP};

end