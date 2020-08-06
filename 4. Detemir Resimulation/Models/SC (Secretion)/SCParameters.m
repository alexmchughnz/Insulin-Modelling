global SC


SC.k = @Getk1k2k3;

function [k1, k2, k3] = Getk1k2k3(P)

if P.data.BMI > 30
    CHalfLife1 = 4.55;       % Half-life of C-peptide in compartment 1 [min]
    F = 0.78;
else
    CHalfLife1 = 4.55;       % Half-life of C-peptide in compartment 1 [min]
    F = 0.78;
end

CHalfLife2 = 0.14*P.data.age + 29.2;

a = log(2)/CHalfLife1;
b = log(2)/CHalfLife2;


% Rate constants.
    k2 = F*(b-a) + a;
    k3 = a*b/(k2); % NOTE: Original code had no factor of 1/2; PDD's thesis does.
    k1 = a + b - k2 - k3;
end