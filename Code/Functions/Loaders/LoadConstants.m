function CONST = LoadConstants()

CONST.MGlucose = 180.156;  % Glucose molar mass [g/mol]
CONST.MCPeptide = 3020.29;  % C-peptide molar mass [g/mol]
CONST.mU2pmol = @(mU) mU * 6.0;    % Insulin [mU]  -> Insulin [pmol]
CONST.pmol2mU = @(pmol) pmol / 6.0;  % Insulin [pmol] -> Insulin [mU]
CONST.IU18Factor = 18;  % TODO: Not sure what this yet! Something involving dL?

CONST.COLUMNWISE = 1;
CONST.ROWWISE = 2;
CONST.PAGEWISE = 3;

end