function P = ApplyInsulinLossFactor(P, JLK)

% JLK is an optional argument so that this function can be iterated over
% e.g. with LineSearchOptimum.
if exist("JLK", "var")
    P.results.JLK = JLK;
end

% Apply IInput proportion.
P = ScalePatientField(P, P.results.JLK, "data", "vIBolus");
P = AddTrialInputs(P);

end

