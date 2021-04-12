function [P, mockP] = MakeMockPatientSim(P)
% Recipe for basic fitting and forward simulating a patient.
% INPUTS:
%   P  - patient struct
% OUTPUT:
%   P  - updated patient struct

global CONFIG

%% Plots
plots = DebugPlots();

plots.EstimateInsulinSecretion.Uen = 1;
plots.EstimateInsulinSecretion.CPep = 1;
plots.FitHepaticClearance.GraphicalID = 1;
plots.FitInsulinSensitivity.SI = 1;
plots.SolveSystem.Glucose = 1;
plots.SolveSystem.Insulin = 1;

DebugPlots(plots);

%% Functions
% First run a basic sim.
P = SimpleSim(P);

% Now sample model results at times of measurements.
[tI, vI] = GetResultsSample(P, P.data.I.time, P.results.I);
[tG, vG] = GetResultsSample(P, P.data.G.time, P.results.G);

% Produce a basic patient struct with simulated data.
mockP = rmfield(P, "results");
mockP.results.tArray = P.results.tArray;

mockP.source = "Mock" + P.source;
mockP.patientCode = CONFIG.PATIENTFORMAT(mockP);
mockP.patientNum = 1000 + P.patientNum;

mockP.data = P.data;
mockP.data.I.value = vI;
mockP.data.G.value = vG;

SavePatients({mockP}, CONFIG.DATAPATH);

% Now run that basic sim again for plots!
mockP = SimpleSim(mockP);





end

