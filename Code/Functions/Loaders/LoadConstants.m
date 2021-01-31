function C = LoadConstants()

C.MGlucose = 180.156;  % Glucose molar mass [g/mol]
C.MCPeptide = 3020.29;  % C-peptide molar mass [g/mol]
C.mU2pmol = @(mU) mU * 6.0;    % Insulin [mU]  -> Insulin [pmol]
C.pmol2mU = @(pmol) pmol / 6.0;  % Insulin [pmol] -> Insulin [mU]
C.IU18Factor = 18;  % TODO: Not sure what this yet! Something involving dL?

end