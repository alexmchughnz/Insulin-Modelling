function PArray = MatchnLSim(P)
% Recipe for adjusting Uen and counter balancing by changing nL rate.
% INPUTS:
%   P  - patient struct
% OUTPUT:
%   PArray  - updated patient structs

%% Plots
plots = DebugPlots();

DebugPlots(plots);


%% Setup
P = EstimateInsulinSecretion(P);
P = FitHepaticClearance(P);
P = FindGutEmptyingRate(P);
P = FitInsulinSensitivity(P);

%% Functions

% Find nL value that minimises error in insulin.
IInputScales = 1.00 + [0 -0.15 -0.30];
numI = length(IInputScales);

UenScales = 1.00 + [-0.15 0 +0.15];
numU = length(UenScales);

numRuns = numI*numU;

runtime = tic;
for ii = 1:numI
    iP = P;
    
    % Adjust Iinput values.
    iP = ScalePatientField(IInputScales(ii), iP, "data", "vIBolus");
    iP = MakeBolusFunctions(iP);
        
    for uu = 1:numU        
        % Edit patient.
        uP = iP;
        uP.patientNum = 1000*uP.patientNum + 10*uu + ii;
        
        % Adjust Uen values.
        uP = ScalePatientField(UenScales(uu), uP, "results", "Uen");
        
        % Simulate for optimal insulin error.       
        GetInsulinError = MakeInsulinErrorFunc(uP);   
        lowerBound = 0.00;
        upperBound = 2.00;
        [IInputScale, IError] = fminbnd(GetInsulinError, lowerBound, upperBound); 
        IScales(uu, ii) = IInputScale;    
        IErrors(uu, ii) = IError;
        
        count = (ii-1)*numI + uu;  
        runtime = PrintTimeRemaining("MatchnLSim", runtime, count, numRuns, P);
    end
end

% Plot patient and save results.
SolveSystem(P, true);
P.results.MatchnL.IScales = IScales;
P.results.MatchnL.IErrors = IErrors;
P.results.MatchnL.IInputScales = IInputScales;
P.results.MatchnL.UenScales = UenScales;

% Plot optimal fit for unchanged parameters.
uuOpt = find(UenScales == 1.0);
iiOpt = find(IInputScales == 1.0);
optimalIScale = IScales(uuOpt, iiOpt);
optimalP = ScalePatientField(optimalIScale, P, "data", "vIBolus");
optimalP = MakeBolusFunctions(optimalP);
SolveSystem(optimalP, true);

PArray = {P optimalP};

end


function insulinErrorFunc = MakeInsulinErrorFunc(P)

    function error = GetInsulinError(IInputScale)
        % Adjust insulin input.
        P = ScalePatientField(IInputScale, P, "data", "vIBolus");
        P = MakeBolusFunctions(P);
        
        % Retrieve insulin fit error.
        P = SolveSystem(P, false);
        error = P.results.fits.insulinSSE;
    end

insulinErrorFunc = @GetInsulinError;
end
