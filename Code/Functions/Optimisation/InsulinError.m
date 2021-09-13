function error = InsulinError(P)

% Simulate I(t).
P = SolveSystem(P);  

% Get RMS error of simulated G to measured data.
[tI, vI] = GetData(P.data.I);
[~, simI] = GetResultsSample(P, tI, P.results.I);

error = rms(simI - vI);

end