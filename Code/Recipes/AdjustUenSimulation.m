function PArray = AdjustUenSimulation(P)
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
UenProportion = 0.85;
deltaIProportion = 0.1;
IProportion = 1.00 + deltaIProportion;

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
resultP = FitInsulinSensitivity(adjP, false);
resultSI = resultP.results.SI;

SIGap = abs(targetSI - resultSI);
newGap = Inf;

% Iterate to find the IInputProportion that makes adjP's SI equal to P's.
while newGap > 1e-7   
    % Make new copy of patient and adjust IInput and bolus functions.
    copyP = adjP;
    copyP.data.vIBolus = P.data.vIBolus .* IProportion;
    copyP = MakeBolusFunctions(copyP);
    
    % Fit SI and store results.    
    copyP = FitInsulinSensitivity(copyP, false);
    newSI = copyP.results.SI;
    
    % IInput needs to go down if SI is too low, or up if SI is too high.
    newGap = abs(targetSI-newSI);
    isSIImproving =  newGap < SIGap;
    relativeGapSize = newGap / SIGap;
    
    if isSIImproving
        deltaIProportion = deltaIProportion * relativeGapSize;
    else        
        deltaIProportion = -deltaIProportion * relativeGapSize;
    end    
    IProportion = 1 + deltaIProportion;
end

% Save IProportions.
P.results.IInputProportion = 1.0;
adjP.results.IInputProportion = IProportion;

% Just finishing simulation for plots!
P = SolveSystem(P);
adjP = SolveSystem(adjP);

PArray = {P adjP};

%% Debug Plots
plots.AdjustUenSimulation = struct();
DP = plots.AdjustUenSimulation;

% MakeDebugPlot("SI Adjustments", P, DP);
% 
% plt = plot(basalP.results.tArray, arrayify(basalP, basalSI), ':');
% plt.DisplayName = 'Basal SI';
% 
% for ii = 1:length(IProportion)
%     plt = plot(P.results.tArray, arrayify(P, fullSI), ':');
%     plt.DisplayName = sprintf("Full SI (%d%%)", IProportion(ii));    
% end
% 
% if isfield(P.data, 'tIBolus')
%     plt = line([P.data.tIBolus'; P.data.tIBolus'], ylim, ...
%         'Color', 'r', ...
%         'LineStyle', '--');
%     plt(1).DisplayName = 'Bolus Input';
%     for pp = 2:length(plt)
%         plt(pp).HandleVisibility = 'off';
%     end
% end
% 
% xlabel('Time')
% ylabel('$S_I$ [L/mU/min]')
% 
% legend()

end
