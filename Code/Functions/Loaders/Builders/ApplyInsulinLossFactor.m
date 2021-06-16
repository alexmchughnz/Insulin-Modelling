function P = ApplyInsulinLossFactor(P, JLK)

P.results.JLK = JLK;

% Apply IInput proportion.
P = ScalePatientField(P, P.results.JLK, "data", "vIBolus");
P = AddBolusArrays(P);

% Re-simulate SC model to change QLocal (if applicable).
forceReSim = true;
P = AddPlasmaInsulinInputArray(P, forceReSim); 

end

