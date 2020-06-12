% Parameters
F = 0.76;
CHalfLife1 = 4.95;       % Half-life of C-peptide in compartment 1 [min]
CHalflife2 = 32.4;         % Half-life of C-peptide in compartment 2 [min]
a = log(2)/CHalfLife1;
b = log(2)/CHalflife2;

SC.k2 = F*(b-a) + a;     % Rate constants
SC.k3 = a*b/(2*SC.k2); % NOTE: Original code had no factor of 1/2; PDD's thesis does.
SC.k1 = a + b - SC.k2 - SC.k3;