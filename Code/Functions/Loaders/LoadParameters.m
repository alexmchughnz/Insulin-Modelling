function P = LoadParameters(P)
P.params.CP = CPParameters(P);
P.params.GC = GCParameters(P);
P.params.GI = GIParameters(P);
if P.data.IType == "detemir"
    P.params.ID = IDParameters(P);
elseif P.data.IDelivery == "subcutaneous"
    P.params.SC = SCParameters(P);
end
end


function CP = CPParameters(P)
if P.data.BMI > 30
    CHalfLife1 = 4.55;       % Half-life of C-peptide in compartment 1 [min]
    F = 0.78;
else
    CHalfLife1 = 4.95;       % Half-life of C-peptide in compartment 1 [min]
    F = 0.76;
end

CHalfLife2 = 0.14*P.data.age + 29.2;

a = log(2)/CHalfLife1;
b = log(2)/CHalfLife2;

% Rate constants.
k2 = F*(b-a) + a;
k3 = a*b/(k2); % NOTE: Original code had no factor of 1/2; PDD's thesis does.
k1 = a + b - k2 - k3;

CP.k1 = k1;
CP.k2 = k2;
CP.k3 = k3;
end
