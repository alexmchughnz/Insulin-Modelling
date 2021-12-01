function P = ApplyInsulinLossFactor(P, JLK)

% JLK is an optional argument so that this function can be iterated over
% e.g. with LineSearchOptimum, by setting the JLK field in advance.
if exist("JLK", "var")
    P.results.JLK = JLK;
end

% Apply IInput proportion.
P = ScalePatientField(P, "data.vIBolus", P.results.JLK);
P = AddTrialInputs(P);

end

