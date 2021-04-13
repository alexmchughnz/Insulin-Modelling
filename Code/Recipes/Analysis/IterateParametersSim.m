function PArray = IterateParametersSim(P)
% INPUTS:
%   P  - patient struct
% OUTPUT:
%   PArray  - updated patient structs

%% Plots
plots = DebugPlots();

plots.EstimateInsulinSecretion.Uen = false;
plots.EstimateInsulinSecretion.CPep = false;
plots.FitHepaticClearance.GraphicalID = false;
plots.FitHepaticClearance.Convergence = false;
plots.SolveSystem.Glucose = false;

DebugPlots(plots);

%% Setup
PArray = {};

%% Functions
JLKArray = [0.1 : 0.1 : 1.0];  % "Justified Loss from injection in Knee"
ks3RateArray = [1.0 : 0.2 : 2.0];

for kk = 1:numel(ks3RateArray)
    % Apply change to ks3 parameter.
    kkP = P;
    kkP = ScalePatientField(kkP, ks3RateArray(kk), "parameters", "SC", "ks3");

    for jj = 1:numel(JLKArray)
        % Apply IInput proportion.
        jjP = ScalePatientField(kkP, JLKArray(jj), "data", "vIBolus");
        jjP = AddBolusArrays(jjP);
        
        forceReSim = true;
        jjP = AddPlasmaInsulinInputArray(jjP, forceReSim); % Re-simulate SC model to change QLocal.
        
        % Now run full simulation for patient.
        jjP = SimpleSim(jjP);        
        
        % Find measure of variance due to insulin error.
        [jjP, hasSSE] = GetPersistent(jjP, "stddevSSE");
        if ~hasSSE
            stdDevPc = 5/100;
            N = 1000;
            jjP = AnalyseInsulinVariance(jjP, stdDevPc, N);
        end
        
        % Save.
        PArray{end+1} = jjP;
    end
end
