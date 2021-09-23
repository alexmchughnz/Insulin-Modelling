function P = ApplyInsulinLossFactor(P, JLK)

P.results.JLK = JLK;

% Apply IInput proportion.
P = ScalePatientField(P, P.results.JLK, "data", "vIBolus");
P = AddTrialInputs(P);

end

