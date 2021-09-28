function error = GlucoseError(P)

% Simulate G(t).
P = GCModel(P);

% Get RMS error of simulated G to measured data.
[tG, vG] = GetData(P.data.G);
[~, simG] = GetResultsSample(P, tG, P.results.G);

error = rms(simG - vG);

end