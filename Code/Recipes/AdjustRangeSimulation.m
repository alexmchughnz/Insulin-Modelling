function PArray = AdjustRangeSimulation(P)
% Recipe for basic fitting and forward simulating a patient.
% INPUTS:
%   P  - patient struct
% OUTPUT:
%   PArray  - updated patient structs

%% Plots
plots = DebugPlots();

plots.SolveSystem.Glucose = 1;
plots.SolveSystem.Insulin = 1;

DebugPlots(plots);


%% Setup
adjP = P;
adjP.patientCode = P.patientCode+"-adj";

% Simulate patient as is.
P = EstimateInsulinSecretion(P);
P = FitHepaticClearance(P);
P = FindGutEmptyingRate(P);
P = FitInsulinSensitivity(P, false);
targetSI = P.results.SI


%% Functions
% Simulate patient with adjusted Uen.
UenProportion = 0.01;

adjP = EstimateInsulinSecretion(adjP);
adjP.results.Uen = adjP.results.Uen .* UenProportion;
adjP = FitHepaticClearance(adjP);
adjP = FindGutEmptyingRate(adjP);
adjP = FitInsulinSensitivity(adjP, false);
resultSI = adjP.results.SI

SolveSystem(P);
SolveSystem(adjP);


PArray = {P adjP};

%% Debug Plots
% plots.SimInsulinUptake = struct();
% DP = plots.SimInsulinUptake;
% 
% arrayify = @(P, value) value .* ones(size(P.results.tArray));
% 
% 
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
